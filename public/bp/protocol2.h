#pragma once

#include "./utils/fst.h"
#include "ecc/ecc.h"
#include "misc/misc.h"

// protocol 2
// for a given P, the prover proves that it has vectors a, b for which
// P = g*a + h*b + u*<a,b>

namespace bp {

struct RoundX {
  Fr x;
  Fr inv;
  Fr square;
  Fr square_inv;
};

void GenerateRoundX(h256_t const& seed, RoundX& round_x) {
  round_x.x = H256ToFr(seed);
  Fr::inv(round_x.inv, round_x.x);
  Fr::sqr(round_x.square, round_x.x);
  Fr::inv(round_x.square_inv, round_x.square);
}

void UpdateSeed(h256_t& seed, G1 const& p, G1 const& u, size_t size) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, p);
  HashUpdate(hash, u);
  HashUpdate(hash, size);
  hash.Final(seed.data());
}

void UpdateSeed(h256_t& seed, G1 const& L, G1 const& R) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, L);
  HashUpdate(hash, R);
  hash.Final(seed.data());
}

inline G1 MultiExpGH(G1 const* g, Fr const* a, G1 const* h, Fr const* b,
                     size_t n) {
  auto get_g = [g, h, n](size_t i) -> G1 const& {
    return i < n ? g[i] : h[i - n];
  };
  auto get_f = [a, b, n](size_t i) -> Fr const& {
    return i < n ? a[i] : b[i - n];
  };
  return MultiExpBdlo12<G1>(get_g, get_f, n * 2, true);
}

struct Protocol2Proof {
  // G1 p;
  std::vector<G1> left;  // size = log(g.size())
  std::vector<G1> right;
  Fr a;
  Fr b;
};

inline bool operator==(Protocol2Proof const& left,
                       Protocol2Proof const& right) {
  return left.left == right.left && left.right == right.right &&
         left.a == right.a && left.b == right.b;
}

inline bool operator!=(Protocol2Proof const& left,
                       Protocol2Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Protocol2Proof const& t) {
  ar& YAS_OBJECT_NVP("bp.p2proof", ("l", t.left), ("r", t.right), ("a", t.a),
                     ("b", t.b));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Protocol2Proof& t) {
  ar& YAS_OBJECT_NVP("bp.p2proof", ("l", t.left), ("r", t.right), ("a", t.a),
                     ("b", t.b));
}

Protocol2Proof Protocol2Prove(h256_t seed, G1 p, G1 const& u,
                              std::vector<G1>&& g, std::vector<G1>&& h,
                              std::vector<Fr>&& a, std::vector<Fr>&& b,
                              Fr const& c) {
  Tick tick(__FUNCTION__);

  assert(g.size() == h.size());
  assert(g.size() == a.size());
  assert(g.size() == b.size());
  assert(InnerProduct(a, b) == c);
  assert(MultiExpBdlo12(g, a) + MultiExpBdlo12(h, b) + u * c == p);
  UpdateSeed(seed, p, u, a.size());

  auto n = g.size();
  auto align_n = misc::Pow2UB(n);
  if (align_n > n) {
    g.resize(align_n);
    h.resize(align_n);
    a.resize(align_n);
    b.resize(align_n);
    std::fill(g.begin() + n, g.end(), G1Zero());
    std::fill(h.begin() + n, h.end(), G1Zero());
    std::fill(a.begin() + n, a.end(), FrZero());
    std::fill(b.begin() + n, b.end(), FrZero());
  }

  Protocol2Proof proof;
  auto rounds = misc::Log2UB(align_n);
  proof.left.resize(rounds);
  proof.right.resize(rounds);
  RoundX round_x;
  for (size_t loop = 0; loop < rounds; ++loop) {
    GenerateRoundX(seed, round_x);
    auto nn = g.size() / 2;

    G1 L, R;
    std::array<parallel::Task, 2> tasks;
    tasks[0] = [&L, &g, &a, &h, &b, &u, nn]() {
      L = MultiExpGH(&g[nn], &a[0], &h[0], &b[nn], nn);
      auto CL = InnerProduct(&a[0], &b[nn], nn);
      L += u * CL;
    };
    tasks[1] = [&R, &g, &a, &h, &b, &u, nn]() {
      R = MultiExpGH(&g[0], &a[nn], &h[nn], &b[0], nn);
      auto CR = InnerProduct(&a[nn], &b[0], nn);
      R += u * CR;
    };
    parallel::Invoke(tasks);

    std::vector<G1> gg(nn);
    auto parallel_f = [&round_x, &gg, &g, nn](int64_t i) {
      gg[i] = MultiExp(g[i], round_x.inv, g[nn + i], round_x.x);
    };
    parallel::For((int64_t)nn, parallel_f);

    std::vector<G1> hh(nn);
    auto parallel_f2 = [&round_x, &hh, &h, nn](int64_t i) {
      hh[i] = MultiExp(h[i], round_x.x, h[nn + i], round_x.inv);
    };
    parallel::For((int64_t)nn, parallel_f2);

    auto pp = p + MultiExp(L, round_x.square, R, round_x.square_inv);

    std::vector<Fr> aa(nn);
    auto parallel_f3 = [&aa, &a, &round_x, nn](int64_t i) {
      aa[i] = a[i] * round_x.x + a[nn + i] * round_x.inv;
    };
    parallel::For((int64_t)nn, parallel_f3);

    std::vector<Fr> bb(nn);
    auto parallel_f4 = [&bb, &b, &round_x, nn](int64_t i) {
      bb[i] = b[i] * round_x.inv + b[nn + i] * round_x.x;
    };
    parallel::For((int64_t)nn, parallel_f4);

    g.swap(gg);
    h.swap(hh);
    p = pp;
    a.swap(aa);
    b.swap(bb);

    proof.left[loop] = L;
    proof.right[loop] = R;

    UpdateSeed(seed, L, R);
  }

  assert(g.size() == 1 && h.size() == 1);
  assert(a.size() == 1 && b.size() == 1);

  proof.a = a[0];
  proof.b = b[0];

  return proof;
}

bool Protocol2Verify(h256_t seed, G1 p, G1 const& u, std::vector<G1>&& g,
                     std::vector<G1>&& h, Protocol2Proof const& proof) {
  Tick tick(__FUNCTION__);
  assert(g.size() == h.size());
  UpdateSeed(seed, p, u, g.size());

  auto n = g.size();
  auto align_n = misc::Pow2UB(n);
  if (align_n > n) {
    g.resize(align_n);
    h.resize(align_n);
    std::fill(g.begin() + n, g.end(), G1Zero());
    std::fill(h.begin() + n, h.end(), G1Zero());
  }
  auto rounds = misc::Log2UB(align_n);
  if (proof.left.size() != rounds || proof.right.size() != rounds) {
    assert(false);
    return false;
  }

  std::vector<RoundX> rounds_x(rounds);
  for (size_t loop = 0; loop < rounds; ++loop) {
    GenerateRoundX(seed, rounds_x[loop]);
    UpdateSeed(seed, proof.left[loop], proof.right[loop]);
  }

  // true: 1; false: -1
  auto get_b = [rounds](size_t i, size_t j) -> bool {
    auto pow_j = (size_t)1 << (rounds - 1 - j);
    return (i & pow_j) ? true : false;
  };

  std::vector<Fr> ss(align_n);
  std::vector<Fr> ss_inverse(align_n);
  auto parallel_f = [&ss, rounds, &get_b, &rounds_x, &ss_inverse](size_t i) {
    ss[i] = FrOne();
    for (size_t j = 0; j < rounds_x.size(); ++j) {
      auto b = get_b(i, j);
      assert(!rounds_x[j].x.isZero());
      ss[i] = ss[i] * (b ? rounds_x[j].x : rounds_x[j].inv);
    }
  };
  parallel::For(align_n, parallel_f);
  ss_inverse = ss;
  FrInv(ss_inverse.data(), ss_inverse.size());

  G1 last_g, last_h;

  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&last_g, &g, &ss]() {
    last_g = MultiExpBdlo12(&g[0], &ss[0], g.size());
  };
  tasks[1] = [&last_h, &h, &ss_inverse]() {
    last_h = MultiExpBdlo12(&h[0], &ss_inverse[0], h.size());
  };
  parallel::Invoke(tasks);

  G1 out = MultiExp(last_g, proof.a, last_h, proof.b);

  out += u * (proof.a * proof.b);

  if (rounds) {
    auto get_g = [&proof](size_t i) -> G1 const& {
      auto size = proof.left.size();
      return i < size ? proof.left[i] : proof.right[i - size];
    };
    auto get_f = [&rounds_x](size_t i) -> Fr const& {
      auto size = rounds_x.size();
      return i < size ? rounds_x[i].square : rounds_x[i - size].square_inv;
    };
    p += MultiExpBdlo12<G1>(get_g, get_f, rounds * 2);
  }

  bool ret = out == p;
  assert(ret);
  return ret;
}

inline bool TestProtocol2(int64_t n) {
  Tick tick(__FUNCTION__);
  h256_t seed = misc::RandH256();
  std::vector<G1> g(n);
  G1Rand(g.data(), n);

  std::vector<G1> h(n);
  G1Rand(h.data(), n);

  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  std::vector<Fr> b(n);
  FrRand(b.data(), n);

  Fr c = InnerProduct(a, b);
  G1 u = G1Rand();
  G1 p = MultiExpBdlo12(g, a) + MultiExpBdlo12(h, b) + u * c;

  auto g2 = g;
  auto h2 = h;
  auto proof = Protocol2Prove(seed, p, u, std::move(g), std::move(h),
                              std::move(a), std::move(b), c);

  bool success =
      Protocol2Verify(seed, p, u, std::move(g2), std::move(h2), proof);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace bp