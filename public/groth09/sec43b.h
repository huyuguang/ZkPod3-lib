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
  };

  struct CommitmentSec {
    std::vector<Fr> r;  // r.size = m
    std::vector<Fr> s;  // s.size = m
    std::vector<Fr> t;  // t.size = m
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
    std::vector<std::vector<Fr>> x;  // m*n
    std::vector<std::vector<Fr>> y;
    std::vector<std::vector<Fr>> z;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 const& get_gz;
    size_t n_ = 0;

    int64_t m() const { return x.size(); }
    int64_t n() const { return n_; }
    std::string to_string() const {
      return std::to_string(m()) + "*" + std::to_string(n());
    }

    void Take(std::vector<std::vector<Fr>>& ox,
              std::vector<std::vector<Fr>>& oy,
              std::vector<std::vector<Fr>>& oz) {
      ox = std::move(x);
      oy = std::move(y);
      oz = std::move(z);
    }

    ProveInput(std::vector<std::vector<Fr>>&& x,
               std::vector<std::vector<Fr>>&& y,
               std::vector<std::vector<Fr>>&& z, GetRefG1 const& get_gx,
               GetRefG1 const& get_gy, GetRefG1 const& get_gz)
        : x(std::move(x)),
          y(std::move(y)),
          z(std::move(z)),
          get_gx(get_gx),
          get_gy(get_gy),
          get_gz(get_gz) {
      Check();
    }

   private:
    void Check() {
      CHECK(!x.empty(), "");
      CHECK(x.size() == y.size() && x.size() == z.size(), "");
      for (auto i = 0LL; i < m(); ++i) {
        CHECK(x[i].size() == y[i].size() && x[i].size() == z[i].size(), "");
        n_ = std::max(n_, x[i].size());
      }
    }
  };

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                         ProveInput const& input) {
    Tick tick(__FN__, input.to_string());
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
      std::array<parallel::VoidTask, 3> tasks;
      tasks[0] = [&com_pub, &input, &com_sec, i]() {
        com_pub.a[i] = pc::ComputeCom(input.get_gx, input.x[i], com_sec.r[i]);
      };
      tasks[1] = [&com_pub, &input, &com_sec, i]() {
        com_pub.b[i] = pc::ComputeCom(input.get_gy, input.y[i], com_sec.s[i]);
      };
      tasks[2] = [&com_pub, &input, &com_sec, i]() {
        com_pub.c[i] = pc::ComputeCom(input.get_gz, input.z[i], com_sec.t[i]);
      };
      parallel::Invoke(tasks);
    };
    parallel::For(m, parallel_f);
  }

  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub, int64_t m,
                         int64_t n) {
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

  static void Prove(Proof& proof, h256_t seed, ProveInput&& input,
                    CommitmentPub const& com_pub,
                    CommitmentSec const& com_sec) {
    Tick tick(__FN__, input.to_string());

    auto m = input.m();
    auto n = input.n();

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
      Tick tick53(" sec43b->Sec53");

      auto& com_sec_53 = para53.com_sec_53;
      auto& com_pub_53 = para53.com_pub_53;

      auto parallel_f = [&input_x, &k](int64_t i) { input_x[i] *= k[i]; };
      parallel::For(m, parallel_f);

      auto& input_yt = para53.input_yt;

      Fr z = FrZero();

      {
        Tick tickz("Sec53 compute z");
        // 2*m*n fr mul
        std::vector<Fr> ip(m);
        auto pf = [&input_x, &input_yt, &input_y, &t, &ip](int64_t i) {
          input_yt[i] = HadamardProduct(input_y[i], t);
          ip[i] = InnerProduct(input_x[i], input_yt[i]);
        };
        parallel::For(m, pf);
        z = std::accumulate(ip.begin(), ip.end(), FrZero());
        // for (int64_t i = 0; i < m; ++i) {
        //  input_yt[i] = HadamardProduct(input_y[i], t);
        //  z += InnerProduct(input_x[i], input_yt[i]);
        //}
      }

      para53.input_53.reset(new typename Sec53::ProveInput(
          std::move(input_x), std::move(input_y), t, std::move(input_yt), z,
          input.get_gx, input.get_gy, SelectSec53Gz()));

      auto& input_53 = *para53.input_53;

      {
        Tick tickz("Sec53 compute com_sec_53 com_pub_53");
        input_53_z = input_53.z;
        com_sec_53.r.resize(m);
        com_pub_53.a.resize(m);
        auto parallel_f2 = [&com_sec, &com_pub, &com_sec_53, &com_pub_53,
                            &k](int64_t i) {
          com_sec_53.r[i] = com_sec.r[i] * k[i];
          com_pub_53.a[i] = com_pub.a[i] * k[i];
        };
        parallel::For(m, parallel_f2);

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
      Tick tick53("sec43b->hyraxa");
      auto& com_pub_hy = parahy.com_pub_hy;
      auto& com_sec_hy = parahy.com_sec_hy;
      auto& zk = parahy.zk;

      auto parallel_f = [&input_z, &zk, &k, m](int64_t j) {
        for (int64_t i = 0; i < m; ++i) {
          auto const& input_zi = input_z[i];
          if (j < (int64_t)input_zi.size()) {
            zk[j] += input_zi[j] * k[i];
          }
        }
      };
      parallel::For(n, parallel_f, n < 16 * 1024);

      parahy.input_hy.reset(new typename HyraxA::ProveInput(
          "43b", zk, t, input_53_z, input.get_gz, SelectSec53Gz()));
      auto& input_hy = *parahy.input_hy;
      (void)input_hy;

      com_sec_hy.r_xi = InnerProduct(com_sec.t, k);
      com_sec_hy.r_tau = com_sec_53_t;

      com_pub_hy.tau = com_pub_53_c;

      // do not need to compute the com1_pub_hy.xi in release build
      DCHECK(input_hy.y == InnerProduct(input_hy.x, input_hy.a), "");

      com_pub_hy.xi = MultiExpBdlo12(com_pub.c, k);

      DCHECK(pc::ComputeCom(input.get_gz, input_hy.x, com_sec_hy.r_xi) ==
                 com_pub_hy.xi,
             "");
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
    VerifyInput(std::vector<size_t> const& mn, CommitmentPub const& com_pub,
                GetRefG1 const& get_gx, GetRefG1 const& get_gy,
                GetRefG1 const& get_gz)
        : mn(std::move(mn)),
          com_pub(com_pub),
          get_gx(get_gx),
          get_gy(get_gy),
          get_gz(get_gz) {}
    std::vector<size_t> const& mn;
    CommitmentPub const& com_pub;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 const& get_gz;
    size_t m() const { return mn.size(); }
    size_t n() const { return *std::max_element(mn.begin(), mn.end()); }
    std::string to_string() const {
      return std::to_string(m()) + "*" + std::to_string(n());
    }
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FN__, input.to_string());
    auto m = input.m();
    auto n = input.n();

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

      typename Sec53::VerifyInput input_53(input.mn, t, std::move(com_pub_53),
                                           input.get_gx, input.get_gy,
                                           SelectSec53Gz());
      ret_53 = Sec53::Verify(proof.proof_53, seed, std::move(input_53));
      assert(ret_53);
    };

    bool ret_a2 = false;
    tasks[1] = [&ret_a2, &com_pub, &proof, &t, &k, &seed, &input]() {
      typename HyraxA::CommitmentPub com_pub_hy(MultiExpBdlo12(com_pub.c, k),
                                                proof.c);
      typename HyraxA::VerifyInput input_hy("43b", t, com_pub_hy, input.get_gz,
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
  std::cout << Tick::GetIndentString() << "old_m=" << m << ", n=" << n << "\n";

  std::vector<std::vector<Fr>> x(m);
  for (auto& i : x) {
    i.resize(n);
    FrRand(i);
  }

  std::vector<std::vector<Fr>> y(m);
  for (auto& i : y) {
    i.resize(n);
    FrRand(i);
  }

  x.front().resize(n / 3);
  y.front().resize(n / 3);

  x.resize(x.size() + 1);
  x.back().resize(n / 2);
  FrRand(x.back());

  y.resize(y.size() + 1);
  y.back().resize(n / 2);
  FrRand(y.back());
  ++m;

  std::vector<size_t> mn(m);
  for (int64_t i = 0; i < m; ++i) {
    mn[i] = x[i].size();
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

  ProveInput prove_input(std::move(x), std::move(y), std::move(z), get_gx,
                         get_gy, get_gz);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prove_input);

  Proof proof;
  Prove(proof, seed, std::move(prove_input), com_pub, com_sec);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << Tick::GetIndentString()
            << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  CHECK(proof == proof2, "");
#endif

  VerifyInput verify_input(mn, com_pub, get_gx, get_gy, get_gz);
  bool success = Verify(proof, seed, verify_input);
  std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace groth09
