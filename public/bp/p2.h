#pragma once

#include "./details.h"
#include "./types.h"

// protocol 2
// for a given P, the prover proves that it has vectors a, b for which
// q = g*a + h*b + u*<a,b>

namespace bp {
template <typename GET_G, typename GET_F>
P2Proof P2Prove(P1Committment const& p1_committment, GET_G const& get_g,
                GET_F const& get_f, Challenge const& challenge) {
  Tick tick(__FUNCTION__);
  using namespace details;
  P2Proof p2_proof;
  uint64_t count = challenge.count();
  bool zero_u = challenge.zero_u();
  auto g_count = PackGCount(count);
  std::vector<Fr> a, b;
  a.reserve(g_count);
  b.reserve(g_count);
  for (size_t i = 0; i < g_count; ++i) {
    a.push_back(get_f(i));
    b.push_back(get_f(i + g_count));
  }

  std::vector<G1> g, h;
  g.reserve(g_count);
  h.reserve(g_count);
  for (size_t i = 0; i < g_count; ++i) {
    g.push_back(get_g(i));
    h.push_back(get_g(i + g_count));
  }

  p2_proof.q = p1_committment.q(challenge.u());

  G1 p = p2_proof.q;

  p2_proof.left.resize(challenge.x().size());
  p2_proof.right.resize(challenge.x().size());
  for (size_t loop = 0; loop < challenge.x().size(); ++loop) {
    auto const& x = challenge.x()[loop];
    auto const& x_inverse = challenge.x_inverse()[loop];
    auto const& x_square = challenge.x_square()[loop];
    auto const& x_square_inverse = challenge.x_square_inverse()[loop];

    auto nn = g.size() / 2;

    G1 L, R;

    // L = MultiExpGH(&g[nn], &a[0], &h[0], &b[nn], nn);
    // R = MultiExpGH(&g[0], &a[nn], &h[nn], &b[0], nn);
    std::array<parallel::Task, 2> tasks;
    tasks[0] = [&L, &g, &a, &h, &b, nn]() {
      L = details::MultiExpGH(&g[nn], &a[0], &h[0], &b[nn], nn);
    };
    tasks[1] = [&R, &g, &a, &h, &b, nn]() {
      R = details::MultiExpGH(&g[0], &a[nn], &h[nn], &b[0], nn);
    };
    parallel::Invoke(tasks);

    if (!zero_u) {
      auto CL = InnerProduct(&a[0], &b[nn], nn);
      auto CR = InnerProduct(&a[nn], &b[0], nn);
      L += challenge.u() * CL;
      R += challenge.u() * CR;
    }

    std::vector<G1> gg;
    gg.resize(nn);
    auto parallel_f = [&x_inverse, &gg, &g, &x, nn](int64_t i) {
      gg[i] = MultiExp(g[i], x_inverse, g[nn + i], x);
    };
    parallel::For((int64_t)nn, parallel_f);

    std::vector<G1> hh;
    hh.resize(nn);
    auto parallel_f2 = [&x_inverse, &hh, &h, &x, nn](int64_t i) {
      hh[i] = MultiExp(h[i], x, h[nn + i], x_inverse);
    };
    parallel::For((int64_t)nn, parallel_f2);

    auto pp = p + MultiExp(L, x_square, R, x_square_inverse);

    std::vector<Fr> aa;
    aa.resize(nn);
    auto parallel_f3 = [&aa, &a, &x_inverse, &x, nn](int64_t i) {
      aa[i] = a[i] * x + a[nn + i] * x_inverse;
    };
    parallel::For((int64_t)nn, parallel_f3);

    std::vector<Fr> bb;
    bb.resize(nn);
    auto parallel_f4 = [&bb, &b, &x_inverse, &x, nn](int64_t i) {
      bb[i] = b[i] * x_inverse + b[nn + i] * x;
    };
    parallel::For((int64_t)nn, parallel_f4);

    g.swap(gg);
    h.swap(hh);
    p = pp;
    a.swap(aa);
    b.swap(bb);

    p2_proof.left[loop] = L;
    p2_proof.right[loop] = R;
  }

  assert(g.size() == 1 && h.size() == 1);
  assert(a.size() == 1 && b.size() == 1);

  p2_proof.a = a[0];
  p2_proof.b = b[0];

  return p2_proof;
}

template <typename GET_G>
bool P2Verify(P1Committment const& p1_committment, GET_G const& get_g,
              P2Proof const& p2_proof, Challenge const& challenge) {
  Tick tick(__FUNCTION__);
  using namespace details;
  uint64_t count = challenge.count();
  auto g_count = PackGCount(count);
  auto x_count = PackXCount(count);
  // std::cout << "proof.p: " << proof.p << std::endl;
  // std::cout << "proof.a: " << proof.a << std::endl;
  // std::cout << "proof.b: " << proof.b << std::endl;
  // std::cout << "challenge: " << challenge_seed << std::endl;

  auto const& u = challenge.u();
  G1 p = p1_committment.q(u);

  assert(p2_proof.left.size() == x_count && p2_proof.right.size() == x_count);
  if (p2_proof.left.size() != x_count || p2_proof.right.size() != x_count)
    return false;

  auto const& x = challenge.x();
  auto const& x_inverse = challenge.x_inverse();
  auto const& x_square = challenge.x_square();
  auto const& x_square_inverse = challenge.x_square_inverse();

  std::vector<G1> g, h;
  g.resize(g_count);
  h.resize(g_count);
  for (size_t i = 0; i < g_count; ++i) {
    g[i] = get_g(i);
    h[i] = get_g(i + g_count);
  }

  // true: 1; false: -1
  auto get_b = [x_count](size_t i, size_t j) -> bool {
    auto pow_j = (size_t)1 << (x_count - 1 - j);
    return (i & pow_j) ? true : false;
  };

  std::vector<Fr> ss;
  std::vector<Fr> ss_inverse;
  ss.resize(g_count);
  ss_inverse.resize(g_count);

  auto parallel_f = [&ss, x_count, &get_b, &x, &x_inverse,
                     &ss_inverse](int64_t i) {
    ss[i] = FrOne();
    for (size_t j = 0; j < x_count; ++j) {
      auto b = get_b(i, j);
      assert(!x[j].isZero());
      ss[i] = ss[i] * (b ? x[j] : x_inverse[j]);
    }
    Fr::inv(ss_inverse[i], ss[i]);
  };
  parallel::For((int64_t)g_count, parallel_f);

  G1 last_g, last_h;

  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&last_g, &g, &ss, g_count]() {
    last_g = MultiExpBdlo12(&g[0], &ss[0], g_count);
  };
  tasks[1] = [&last_h, &h, &ss_inverse, g_count]() {
    last_h = MultiExpBdlo12(&h[0], &ss_inverse[0], g_count);
  };
  parallel::Invoke(tasks);
  
  G1 out = MultiExp(last_g, p2_proof.a, last_h, p2_proof.b);

  out += u * (p2_proof.a * p2_proof.b);

  if (x_count) {    
    p += details::MultiExpGH(&p2_proof.left[0], &x_square[0], &p2_proof.right[0],
                    &x_square_inverse[0], x_count);
  }

  // std::cout << "out: " << out << std::endl;
  // std::cout << "p: " << p << std::endl;
  bool ret = out == p;
  assert(ret);
  return ret;
}

}  // namespace bp