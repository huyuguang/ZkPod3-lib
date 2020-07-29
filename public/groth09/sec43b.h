#pragma once

#include "./details.h"
#include "./sec53b.h"
#include "hyrax/a2.h"
#include "hyrax/a3.h"

// x, y, z: secret matrix<Fr>, size =m*n
// open: a1=com(gx, x1)...am=com(gx, xm)
// open: b1=com(gy, y1)...bm=com(gy, ym)
// open: c1=com(gz, z1)...cm=com(gz, zm)
// prove: z=x o y (o is hadamard product)
// proof size: 2*log(m)+6 G1, 3n+5 Fr
// prove cost: 2*log(m)*mulexp(n)
// verify cost: 2*mulexp(n)
namespace groth09 {

template <typename Sec53, typename HyraxA>
struct Sec43b {
  // input_53's gz can be any value, here just use pc::PcU()
  static G1 const& SelectSec53Gz() { return pc::PcU(); }

  struct CommitmentPub {
    std::vector<G1> a;  // a.size = m
    std::vector<G1> b;  // b.size = m
    std::vector<G1> c;  // c.size = m
    void Align() {
      int64_t old_m = a.size();
      int64_t new_m = (int64_t)misc::Pow2UB(old_m);
      if (new_m > old_m) {
        // G1 const& g0 = G1Zero();
        a.resize(new_m, G1Zero());
        // std::fill(a.begin() + old_m, a.end(), g0);
        b.resize(new_m, G1Zero());
        // std::fill(b.begin() + old_m, b.end(), g0);
        c.resize(new_m, G1Zero());
        // std::fill(c.begin() + old_m, c.end(), g0);
      }
    }
  };

  struct CommitmentSec {
    std::vector<Fr> r;  // r.size = m
    std::vector<Fr> s;  // s.size = m
    std::vector<Fr> t;  // t.size = m
    void Align() {
      int64_t old_m = r.size();
      int64_t new_m = (int64_t)misc::Pow2UB(old_m);
      if (new_m > old_m) {
        // Fr const& f0 = FrZero();
        r.resize(new_m, FrZero());
        // std::fill(r.begin() + old_m, r.end(), f0);
        s.resize(new_m, FrZero());
        // std::fill(s.begin() + old_m, s.end(), f0);
        t.resize(new_m, FrZero());
        // std::fill(t.begin() + old_m, t.end(), f0);
      }
    }
  };

  struct Proof {
    G1 c;
    typename Sec53::Proof proof_53;  // 2*log(m)+4 G1, 2n+3 Fr
    typename HyraxA::Proof proof_a;  // 2 G1, n+2 Fr

    bool operator==(Proof const& right) const {
      return c == right.c && proof_53 == right.proof_53 &&
             proof_a == right.proof_a;
    }
    bool operator!=(Proof const& right) const { return !(*this == right); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("43b.pf", ("c", c), ("53p", proof_53),
                         ("ap", proof_a));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("43b.pf", ("c", c), ("53p", proof_53),
                         ("ap", proof_a));
    }
  };

  struct ProveInput {
    int64_t const original_m;
    std::vector<std::vector<Fr>> x;  // m*n
    std::vector<std::vector<Fr>> y;
    std::vector<std::vector<Fr>> z;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 const& get_gz;

    int64_t m() const { return x.size(); }
    int64_t n() const { return x[0].size(); }
    void Take(std::vector<std::vector<Fr>>& ox,
              std::vector<std::vector<Fr>>& oy,
              std::vector<std::vector<Fr>>& oz) {
      ox = std::move(x);
      oy = std::move(y);
      oz = std::move(z);
    }

    ProveInput(int64_t original_m, std::vector<std::vector<Fr>>&& ix,
               std::vector<std::vector<Fr>>&& iy,
               std::vector<std::vector<Fr>>&& iz, GetRefG1 const& get_gx,
               GetRefG1 const& get_gy, GetRefG1 const& get_gz)
        : original_m(original_m),
          x(std::move(ix)),
          y(std::move(iy)),
          z(std::move(iz)),
          get_gx(get_gx),
          get_gy(get_gy),
          get_gz(get_gz) {
      // Tick tick(__FN__);
      assert(!x.empty());
      assert(x.size() == y.size());
      assert(x.size() == z.size());
      for (auto i = 0LL; i < m(); ++i) {
        assert(x[i].size() == (size_t)n());
        assert(y[i].size() == (size_t)n());
        assert(z[i].size() == (size_t)n());
      }
    }

    // pad some trivial value
    void Align() {
      // Tick tick(__FN__);
      int64_t old_m = m();
      int64_t new_m = (int64_t)misc::Pow2UB(old_m);
      if (old_m == new_m) return;

      auto const& f0 = FrZero();
      x.resize(new_m);
      y.resize(new_m);
      z.resize(new_m);

      std::vector<Fr> vf0(n(), f0);

      for (int64_t i = old_m; i < new_m; ++i) {
        x[i] = vf0;
        y[i] = vf0;
        z[i] = vf0;
        // auto& xi = x[i];
        // xi.resize(n());
        // std::fill(xi.begin(), xi.end(), f0);
        // auto& yi = y[i];
        // yi.resize(n());
        // std::fill(yi.begin(), yi.end(), f0);
        // auto& zi = z[i];
        // zi.resize(n());
        // std::fill(zi.begin(), zi.end(), f0);
      }
    }
  };

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                         ProveInput const& input) {
    Tick tick(__FN__);
    auto const m = input.m();
    com_sec.r.resize(m);
    FrRand(com_sec.r.data(), m);

    com_sec.s.resize(m);
    FrRand(com_sec.s.data(), m);

    com_sec.t.resize(m);
    FrRand(com_sec.t.data(), m);

    com_pub.a.resize(m);
    com_pub.b.resize(m);
    com_pub.c.resize(m);

    auto parallel_f = [&com_sec, &com_pub, &input](int64_t i) {
      com_pub.a[i] = pc::ComputeCom(input.get_gx, input.x[i], com_sec.r[i]);
      com_pub.b[i] = pc::ComputeCom(input.get_gy, input.y[i], com_sec.s[i]);
      com_pub.c[i] = pc::ComputeCom(input.get_gz, input.z[i], com_sec.t[i]);
    };
    parallel::For(m, parallel_f);
  }

  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub, int64_t m,
                         int64_t n) {
    // Tick tick(__FN__);
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_pub.a);
    HashUpdate(hash, com_pub.b);
    HashUpdate(hash, com_pub.c);
    HashUpdate(hash, m);
    HashUpdate(hash, n);
    hash.Final(seed.data());
  }

  static void ComputeChallengeKT(h256_t const& seed, std::vector<Fr>& k,
                                 std::vector<Fr>& t) {
    ComputeFst(seed, "gro09::sec43b::k", k);
    ComputeFst(seed, "gro09::sec43b::t", t);
  }

  // pad some trivial value
  static void AlignData(ProveInput& input, CommitmentPub& com_pub,
                        CommitmentSec& com_sec) {
    // Tick tick(__FN__);
    input.Align();
    com_pub.Align();
    com_sec.Align();
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput&& input,
                    CommitmentPub com_pub, CommitmentSec com_sec) {
    Tick tick(__FN__);
    auto m = input.m();
    auto n = input.n();
    auto original_m = input.original_m;

    std::cout << "m: " << m << ", n:" << n << "\n";

    UpdateSeed(seed, com_pub, m, n);
    std::vector<Fr> k(m);
    std::vector<Fr> t(n);
    ComputeChallengeKT(seed, k, t);

    Fr input_53_z;
    G1 com_pub_53_c;
    Fr com_sec_53_t;

    std::vector<std::vector<Fr>> input_x;
    std::vector<std::vector<Fr>> input_y;
    std::vector<std::vector<Fr>> input_z;
    input.Take(input_x, input_y, input_z);

    struct Para53 {
      Para53(int64_t m) : input_yt(m) {}
      typename Sec53::CommitmentSec com_sec_53;
      typename Sec53::CommitmentPub com_pub_53;
      std::vector<std::vector<Fr>> input_yt;
      std::unique_ptr<typename Sec53::ProveInput> input_53;
    } para53(m);

    {
      // Tick tick53(" sec43b->Sec53");

      auto& com_sec_53 = para53.com_sec_53;
      auto& com_pub_53 = para53.com_pub_53;

      auto parallel_f = [&input_x, &k](int64_t i) { input_x[i] *= k[i]; };
      parallel::For(original_m, parallel_f, original_m < 1024);

      auto& input_yt = para53.input_yt;

      Fr z = FrZero();

      {
        // Tick tickz("Sec53 compute z");
        // 2*m*n fr mul
        for (int64_t i = 0; i < original_m; ++i) {
          input_yt[i] = HadamardProduct(input_y[i], t);
          z += InnerProduct(input_x[i], input_yt[i]);
        }
        for (int64_t i = original_m; i < m; ++i) {
          input_yt[i].resize(n, FrZero());
          // std::fill(input_yt[i].begin(), input_yt[i].end(), FrZero());
        }
      }

      para53.input_53.reset(new typename Sec53::ProveInput(
          std::move(input_x), std::move(input_y), t, std::move(input_yt), z,
          input.get_gx, input.get_gy, SelectSec53Gz()));

      auto& input_53 = *para53.input_53;

      {
        // Tick tickz("Sec53 compute com_sec_53 com_pub_53");
        input_53_z = input_53.z;
        com_sec_53.r.resize(m, FrZero());
        com_pub_53.a.resize(m, G1Zero());
        auto parallel_f2 = [&com_sec, &com_pub, &com_sec_53, &com_pub_53,
                            &k](int64_t i) {
          com_sec_53.r[i] = com_sec.r[i] * k[i];
          com_pub_53.a[i] = com_pub.a[i] * k[i];
        };
        parallel::For(original_m, parallel_f2, original_m < 16 * 1024);
        // std::fill(com_sec_53.r.begin() + original_m, com_sec_53.r.end(),
        //          FrZero());
        // std::fill(com_pub_53.a.begin() + original_m, com_pub_53.a.end(),
        //          G1Zero());

        com_sec_53.s = com_sec.s;
        com_sec_53.t = FrRand();
        com_sec_53_t = com_sec_53.t;

        com_pub_53.b = com_pub.b;
        com_pub_53.c = pc::ComputeCom(input_53.gz, input_53_z, com_sec_53.t);
        proof.c = com_pub_53.c;  // verifier can not compute c by com_pub.c
        com_pub_53_c = com_pub_53.c;
      }
    }

    struct ParaHy {
      ParaHy(int64_t n) : zk(n, FrZero()) {}
      typename HyraxA::CommitmentPub com_pub_hy;
      typename HyraxA::CommitmentSec com_sec_hy;
      std::vector<Fr> zk;
      std::unique_ptr<typename HyraxA::ProveInput> input_hy;
    } parahy(n);

    {
      // Tick tick53("sec43b->hyraxa");
      auto& com_pub_hy = parahy.com_pub_hy;
      auto& com_sec_hy = parahy.com_sec_hy;
      auto& zk = parahy.zk;

      auto parallel_f = [&input_z, &zk, &k, m](int64_t j) {
        for (int64_t i = 0; i < m; ++i) {
          zk[j] += input_z[i][j] * k[i];
        }
      };
      parallel::For(n, parallel_f, n < 16 * 1024);

      parahy.input_hy.reset(new typename HyraxA::ProveInput(
          zk, t, input_53_z, input.get_gz, SelectSec53Gz()));
      auto& input_hy = *parahy.input_hy;
      (void)input_hy;

      com_sec_hy.r_xi = InnerProduct(com_sec.t, k);
      com_sec_hy.r_tau = com_sec_53_t;

      com_pub_hy.tau = com_pub_53_c;

      // do not need to compute the com1_pub_hy.xi in release build
      assert(input_hy.y == InnerProduct(input_hy.x, input_hy.a));

      com_pub_hy.xi = MultiExpBdlo12(com_pub.c, k);

#ifdef _DEBUG
      auto check_xi = pc::ComputeCom(input.get_gz, input_hy.x, com_sec_hy.r_xi);
      assert(check_xi == com_pub_hy.xi);
#endif
    }

    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&para53, &proof, &seed]() {
      Sec53::Prove(proof.proof_53, seed, std::move(*para53.input_53),
                   std::move(para53.com_pub_53), std::move(para53.com_sec_53));
    };

    tasks[1] = [&parahy, &proof, &seed]() {
      HyraxA::Prove(proof.proof_a, seed, std::move(*parahy.input_hy),
                    std::move(parahy.com_pub_hy), std::move(parahy.com_sec_hy));
    };

    parallel::Invoke(tasks);
  }

  struct VerifyInput {
    VerifyInput(int64_t m, int64_t n, CommitmentPub const& com_pub,
                GetRefG1 const& get_gx, GetRefG1 const& get_gy,
                GetRefG1 const& get_gz)
        : m(misc::Pow2UB(m)),
          n(n),
          com_pub(com_pub),
          get_gx(get_gx),
          get_gy(get_gy),
          get_gz(get_gz) {}
    int64_t const m;
    int64_t const n;
    CommitmentPub const& com_pub;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 const& get_gz;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    // Tick tick(__FN__);
    auto m = input.m;
    auto n = input.n;

    auto const& com_pub = input.com_pub;
    UpdateSeed(seed, com_pub, m, n);
    std::vector<Fr> k(m);
    std::vector<Fr> t(n);
    ComputeChallengeKT(seed, k, t);

    std::array<parallel::VoidTask, 2> tasks;
    bool ret_53 = false;
    tasks[0] = [&ret_53, &proof, &input, m, &com_pub, &k, &t, &seed]() {
      typename Sec53::CommitmentPub com_pub_53;
      com_pub_53.c = proof.c;
      com_pub_53.b = input.com_pub.b;
      com_pub_53.a.resize(m);
      auto parallel_f = [&com_pub_53, &com_pub, &k](int64_t i) {
        com_pub_53.a[i] = com_pub.a[i] * k[i];
      };
      parallel::For(m, parallel_f, m < 1024);

      typename Sec53::VerifyInput input_53(t, com_pub_53, input.get_gx,
                                           input.get_gy, SelectSec53Gz());
      ret_53 = Sec53::Verify(proof.proof_53, seed, input_53);
      assert(ret_53);
    };

    bool ret_a2 = false;
    tasks[1] = [&ret_a2, &com_pub, &proof, &t, &k, &seed, &input]() {
      typename HyraxA::CommitmentPub com_pub_hy(MultiExpBdlo12(com_pub.c, k),
                                                proof.c);
      typename HyraxA::VerifyInput input_hy(t, com_pub_hy, input.get_gz,
                                            SelectSec53Gz());
      ret_a2 = HyraxA::Verify(proof.proof_a, seed, input_hy);
      assert(ret_a2);
    };

    parallel::Invoke(tasks);

    if (!ret_53 || !ret_a2) {
      std::cout << "ret_53: " << ret_53 << ", ret_a2: " << ret_a2 << "\n";
      assert(false);
    }
    return ret_53 && ret_a2;
  }

  static bool Test(int64_t m, int64_t n);
};

template <typename Sec53, typename HyraxA>
bool Sec43b<Sec53, HyraxA>::Test(int64_t m, int64_t n) {
  Tick tick(__FN__);
  std::cout << "old_m=" << m << ", n=" << n << "\n";

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

  std::vector<std::vector<Fr>> z(m);
  for (int64_t i = 0; i < m; ++i) {
    z[i] = HadamardProduct(x[i], y[i]);
  }

  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 0;
  int64_t y_g_offset = 0;
  int64_t z_g_offset = 0;
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  GetRefG1 get_gy = [y_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[y_g_offset + i];
  };
  GetRefG1 get_gz = [z_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[z_g_offset + i];
  };

  ProveInput prove_input(m, std::move(x), std::move(y), std::move(z), get_gx,
                         get_gy, get_gz);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prove_input);

  AlignData(prove_input, com_pub, com_sec);

  Proof proof;
  Prove(proof, seed, std::move(prove_input), com_pub, com_sec);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  if (proof != proof2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  VerifyInput verify_input(m, n, com_pub, get_gx, get_gy, get_gz);
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace groth09
