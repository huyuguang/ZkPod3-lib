#pragma once

#include <map>
#include <memory>

#include "./details.h"
#include "./sec51b.h"
#include "utils/fst.h"

// t is a public vector which size = n
// x, y are secret matric which size = m*n
// commit {x1}, {y1}...., z
// open: com(gx, x1), com(gx, x2)...com(gy, y1)... com(gz, z)
// prove z = <x1,y1 o t> + <x2,y2 o t>...
namespace groth09 {

struct Sec52b {
  struct CommitmentPub {
    std::vector<G1> a;  // a.size = m
    std::vector<G1> b;  // b.size = m
    G1 c;
    int64_t m() const { return a.size(); }
    bool CheckFormat() const { return !a.empty() && (a.size() == b.size()); }
  };

  struct CommitmentSec {
    std::vector<Fr> r;  // r.size = m
    std::vector<Fr> s;  // s.size = m
    Fr t;
  };

  struct VerifyInput {
    VerifyInput(std::vector<Fr> const& t, CommitmentPub const& com_pub,
                GetRefG1 const& get_gx, GetRefG1 const& get_gy, G1 const& gz)
        : t(t), com_pub(com_pub), get_gx(get_gx), get_gy(get_gy), gz(gz) {}
    std::vector<Fr> const& t;
    CommitmentPub const& com_pub;
    int64_t m() const { return com_pub.m(); }
    bool CheckFormat(int64_t check_n) const {
      if (t.size() != (uint64_t)check_n) return false;
      return com_pub.CheckFormat();
    }
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    G1 const& gz;
  };

  struct ProveInput {
    std::vector<std::vector<Fr>> const& x;   // m*n
    std::vector<std::vector<Fr>> const& y;   // m*n
    std::vector<Fr> const& t;                // n
    std::vector<std::vector<Fr>> const& yt;  // m*n
    Fr const z;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    G1 const& gz;
    std::vector<Fr> sum_xy;  // m * 2 - 1

    int64_t m() const { return x.size(); }
    int64_t n() const { return x[0].size(); }

    ProveInput(std::vector<std::vector<Fr>> const& x,
               std::vector<std::vector<Fr>> const& y, std::vector<Fr> const& t,
               std::vector<std::vector<Fr>> const& yt, Fr const& z,
               GetRefG1 const& get_gx, GetRefG1 const& get_gy, G1 const& gz)
        : x(x),
          y(y),
          t(t),
          yt(yt),
          z(z),
          get_gx(get_gx),
          get_gy(get_gy),
          gz(gz),
          sum_xy(BuildSumXY()) {
#ifdef _DEBUG
      assert(x.size() == y.size() && !x.empty());
      assert(x.size() == yt.size());
      assert(t.size() == (size_t)n());
      Fr check_z = FrZero();
      for (int64_t i = 0; i < m(); ++i) {
        assert(x[i].size() == (size_t)n());
        assert(y[i].size() == (size_t)n());
        assert(yt[i].size() == (size_t)n());
        assert(yt[i] == HadamardProduct(y[i], t));
        check_z += InnerProduct(x[i], yt[i]);
      }
      assert(z == check_z);
#endif
    }

   private:
    Fr ComputeSumOfXY(int64_t l) {
      Fr ret = FrZero();
      // l = m + i - j - 1
      int64_t min_j = m() - l - 1;
      int64_t max_j = m() * 2 - l - 1;
      if (min_j < 0) min_j = 0;
      if (max_j > m()) max_j = m();
      for (int64_t j = min_j; j < max_j; ++j) {
        int64_t i = l + j + 1LL - m();
        assert(i >= 0 && i < m());
        ret += InnerProduct(x[i], yt[j]);
      }
      return ret;
    }

    std::vector<Fr> BuildSumXY() {
      Tick tick(__FN__);
      std::vector<Fr> ret(m() * 2 - 1);

      auto parallel_f = [this, &ret](int64_t l) { ret[l] = ComputeSumOfXY(l); };
      parallel::For(m() * 2 - 1, parallel_f);
      return ret;
    }
  };

  struct CommitmentExtPub {
    std::vector<G1> cl;  // cl.size = 2m-1
    bool CheckFormat(int64_t check_m) const {
      return (int64_t)cl.size() == check_m * 2 - 1;
    }
    int64_t m() const { return (cl.size() + 1) / 2; }
  };

  struct CommitmentExtSec {
    std::vector<Fr> tl;  // tl.size = 2m
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;
    Sec51b::Proof proof_51;
    int64_t n() const { return proof_51.n(); }
    int64_t m() const { return com_ext_pub.m(); }

    bool CheckFormat(int64_t check_m) const {
      return com_ext_pub.CheckFormat(check_m) && proof_51.CheckFormat();
    }
  };

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                         ProveInput const& input) {
    Tick tick(__FN__);
    auto const m = input.m();
    // auto const n = input.n();
    com_sec.r.resize(m);
    FrRand(com_sec.r.data(), m);

    com_sec.s.resize(m);
    FrRand(com_sec.s.data(), m);

    com_sec.t = FrRand();

    com_pub.a.resize(m);
    com_pub.b.resize(m);

    auto parallel_f = [&input, &com_pub, &com_sec](int64_t i) {
      com_pub.a[i] = pc::ComputeCom(input.get_gx, input.x[i], com_sec.r[i]);
      com_pub.b[i] = pc::ComputeCom(input.get_gy, input.y[i], com_sec.s[i]);
    };
    parallel::For(m, parallel_f);

    com_pub.c = pc::ComputeCom(input.gz, input.z, com_sec.t);
  }

  static void ComputeComExt(CommitmentExtPub& com_ext_pub,
                            CommitmentExtSec& com_ext_sec,
                            ProveInput const& input,
                            CommitmentPub const& com_pub,
                            CommitmentSec const& com_sec) {
    Tick tick(__FN__);

    auto const m = input.m();
    com_ext_sec.tl.resize(m * 2);
    FrRand(com_ext_sec.tl.data(), m * 2);
    com_ext_sec.tl[m - 1] = com_sec.t;

    com_ext_pub.cl.resize(m * 2 - 1);

    auto parallel_f = [&input, &com_ext_pub, &com_ext_sec](int64_t l) {
      com_ext_pub.cl[l] =
          pc::ComputeCom(input.gz, input.sum_xy[l], com_ext_sec.tl[l]);
    };
    parallel::For(m * 2 - 1, parallel_f);

    CHECK(com_ext_pub.cl[m - 1] == com_pub.c, "");
  }

  // compute challenge1 by commitment1 and commitment2
  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                         CommitmentExtPub const& com_ext_pub) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_pub.a);
    HashUpdate(hash, com_pub.b);
    HashUpdate(hash, com_pub.c);
    HashUpdate(hash, com_ext_pub.cl);
    hash.Final(seed.data());
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input,
                    CommitmentPub const& com_pub,
                    CommitmentSec const& com_sec) {
    Tick tick(__FN__);
    auto m = input.m();
    auto n = input.n();

    CommitmentExtSec com_ext_sec;
    ComputeComExt(proof.com_ext_pub, com_ext_sec, input, com_pub, com_sec);

    UpdateSeed(seed, com_pub, proof.com_ext_pub);
    Fr e = H256ToFr(seed);

    std::vector<Fr> e_pow;
    std::vector<Fr> e_pow_reverse;
    ComputePowOfE(e, m, e_pow, e_pow_reverse);

    std::vector<Fr> x(n, FrZero());
    // std::fill(x.begin(), x.end(), FrZero());

    auto parallel_fx = [m, &x, &e_pow, &input](int64_t j) {
      for (int64_t i = 0; i < m; ++i) {
        x[j] += e_pow[i] * input.x[i][j];
      }
    };
    parallel::For(n, parallel_fx);

    std::vector<Fr> y(n, FrZero());
    // std::fill(y.begin(), y.end(), FrZero());

    auto parallel_fy = [m, &y, &e_pow_reverse, &input](int64_t j) {
      for (int64_t i = 0; i < m; ++i) {
        y[j] += e_pow_reverse[i] * input.y[i][j];
      }
    };
    parallel::For(n, parallel_fy);

    auto const& t = input.t;
    std::vector<Fr> yt = HadamardProduct(y, t);
    auto z = InnerProduct(x, yt);
    Sec51b::ProveInput input_51(x, y, t, yt, z, input.get_gx, input.get_gy,
                                input.gz);

    assert(z == InnerProduct(input.sum_xy, e_pow));

    Sec51b::CommitmentPub com_pub_51;
    Sec51b::CommitmentSec com_sec_51;
    com_pub_51.a = MultiExpBdlo12(com_pub.a, e_pow);
    com_pub_51.b = MultiExpBdlo12(com_pub.b, e_pow_reverse);
    com_pub_51.c = MultiExpBdlo12(proof.com_ext_pub.cl, e_pow);

    com_sec_51.r = FrZero();
    for (int64_t i = 0; i < m; ++i) {
      com_sec_51.r += e_pow[i] * com_sec.r[i];
    }

    com_sec_51.s = FrZero();
    for (int64_t i = 0; i < m; ++i) {
      com_sec_51.s += e_pow_reverse[i] * com_sec.s[i];
    }

    com_sec_51.t = FrZero();
    for (int64_t l = 0; l < m * 2 - 1; ++l) {
      com_sec_51.t += com_ext_sec.tl[l] * e_pow[l];
    }

#ifdef _DEBUG
    auto check_a = pc::ComputeCom(input_51.get_gx, input_51.x, com_sec_51.r);
    assert(check_a == com_pub_51.a);
    auto check_b = pc::ComputeCom(input_51.get_gy, input_51.y, com_sec_51.s);
    assert(check_b == com_pub_51.b);
    auto check_c = pc::ComputeCom(input_51.gz, input_51.z, com_sec_51.t);
    assert(check_c == com_pub_51.c);
#endif

    Sec51b::Prove(proof.proof_51, seed, input_51, com_pub_51, com_sec_51);
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    auto m = proof.m();

    auto const& com_pub = input.com_pub;
    auto const& com_ext_pub = proof.com_ext_pub;

    UpdateSeed(seed, com_pub, com_ext_pub);
    Fr e = H256ToFr(seed);

    std::vector<Fr> e_pow;
    std::vector<Fr> e_pow_reverse;
    ComputePowOfE(e, m, e_pow, e_pow_reverse);

    std::cout << Tick::GetIndentString() << "2 multiexp(" << m << ")\n";
    std::cout << Tick::GetIndentString() << "multiexp(" << 2 * m - 1 << ")\n";

    Sec51b::CommitmentPub com_pub_51;
    std::array<parallel::VoidTask, 3> tasks;
    tasks[0] = [&com_pub_51, &com_pub, e_pow]() {
      com_pub_51.a = MultiExpBdlo12(com_pub.a, e_pow);
    };
    tasks[1] = [&com_pub_51, &com_pub, e_pow_reverse]() {
      com_pub_51.b = MultiExpBdlo12(com_pub.b, e_pow_reverse);
    };
    tasks[2] = [&com_pub_51, &com_ext_pub, e_pow]() {
      com_pub_51.c = MultiExpBdlo12(com_ext_pub.cl, e_pow);
    };
    parallel::Invoke(tasks);

    Sec51b::VerifyInput input_51(input.t, com_pub_51, input.get_gx,
                                 input.get_gy, input.gz);
    return Sec51b::Verify(proof.proof_51, seed, input_51);
  }

  static bool Test(int64_t m, int64_t n) {
    Tick tick(__FN__);
    std::cout << "m=" << m << ", n=" << n << "\n";

    std::vector<std::vector<Fr>> x(m);
    for (auto& i : x) {
      i.resize(n);
      FrRand(i.data(), n);
    }

    std::vector<std::vector<Fr>> y(m);
    for (auto& i : y) {
      i.resize(n);
      FrRand(i.data(), n);
    }

    std::vector<Fr> t(n);
    FrRand(t.data(), t.size());

    h256_t seed = misc::RandH256();

    int64_t x_g_offset = 10;
    int64_t y_g_offset = 550;
    GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
      return pc::PcG()[x_g_offset + i];
    };
    GetRefG1 get_gy = [y_g_offset](int64_t i) -> G1 const& {
      return pc::PcG()[y_g_offset + i];
    };

    std::vector<std::vector<Fr>> yt(m);
    Fr z = FrZero();
    for (int64_t i = 0; i < m; ++i) {
      yt[i] = HadamardProduct(y[i], t);
      z += InnerProduct(x[i], yt[i]);
    }

    ProveInput prove_input(x, y, t, yt, z, get_gx, get_gy, pc::PcU());
    CommitmentPub com_pub;
    CommitmentSec com_sec;
    ComputeCom(com_pub, com_sec, prove_input);

    Proof proof;
    Prove(proof, seed, prove_input, com_pub, com_sec);

    VerifyInput verify_input(t, com_pub, get_gx, get_gy, pc::PcU());
    bool success = Verify(proof, seed, verify_input);
    std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
    return success;
  }

 private:
  static void ComputePowOfE(Fr const& e, int64_t m, std::vector<Fr>& vec,
                            std::vector<Fr>& rev) {
    // Tick tick(__FUNCTION__, std::to_string(m));
    vec.resize(m * 2 - 1);
    rev.resize(m);

    vec[0] = FrOne();
    for (int64_t i = 1; i < m * 2 - 1; ++i) {
      vec[i] = e * vec[i - 1];
    }
    for (int64_t i = 0; i < m; ++i) {
      rev[i] = vec[m - i - 1];
    }
  }
};

}  // namespace groth09