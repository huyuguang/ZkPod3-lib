#pragma once

#include "./a3.h"
#include "./details.h"

// batch a3
// a_i: public vector<Fr>, i\in[0,m-1], a_i.size() maybe neq a_j.size()
// x_i: secret vector<Fr>, i\in[0,m-1], x_i.size() maybe neq x_j.size()
// a_i.size() must eq x_i.size()
// z: secret Fr
// open: com(gx,x), com(gz,z)
// prove: z = \sum_{i=0}^{m}<x_i,a_i>
// proof size: 2*log(m) G1 + (2+2log(n))G1 + 2Fr

namespace hyrax {
struct A4 {
  struct CommitmentPub {
    std::vector<G1> cx;
    G1 cz;
    int64_t m() const { return cx.size(); }
  };

  struct CommitmentSec {
    std::vector<Fr> r;  // r.size = m
    Fr t;
  };

  struct VerifyInput {
    VerifyInput(std::string const& tag, CommitmentPub&& com_pub,
                GetRefG1 const& get_gx, std::vector<std::vector<Fr>>&& a,
                G1 const& gz)
        : tag(tag),
          com_pub(std::move(com_pub)),
          get_gx(get_gx),
          a(std::move(a)),
          gz(gz) {
      Check();
    }
    std::string tag;
    CommitmentPub com_pub;
    GetRefG1 const& get_gx;
    std::vector<std::vector<Fr>> a;
    G1 gz;
    size_t max_n = 0;

    int64_t m() const { return com_pub.m(); }
    int64_t n() const { return max_n; }
    std::string to_string() const {
      return tag + ": " + std::to_string(m()) + "*" + std::to_string(n());
    }

    void SortAndAlign() {
      auto order = GetSortOrder(a);
      PermuteAndAlign(order, com_pub);
      PermuteAndAlign(order, a);
    }

    void Update(Fr const& e) {
      auto m2 = m() / 2;
      std::vector<std::vector<Fr>> a2(m2);
      auto pf = [this, &a2, &e](int64_t i) {
        a2[i] = a[2 * i] * e + a[2 * i + 1];
      };
      parallel::For(m2, pf);
      a2.swap(a);
    }

   private:
    void Check() {
      CHECK(com_pub.m() == (int64_t)a.size(), "");
      for (auto const& i : a) {
        max_n = std::max(max_n, i.size());
      }
    }
  };

  struct ProveInput {
    std::string tag;
    std::vector<std::vector<Fr>> x;
    std::vector<std::vector<Fr>> a;
    Fr z;
    GetRefG1 const& get_gx;
    G1 const& gz;
    size_t max_n = 0;

    int64_t m() const { return (int64_t)x.size(); }
    int64_t n() const { return (int64_t)max_n; }
    std::string to_string() const {
      return tag + ": " + std::to_string(m()) + "*" + std::to_string(n());
    }
    ProveInput(std::string const& tag, std::vector<std::vector<Fr>>&& x,
               std::vector<std::vector<Fr>>&& a, Fr const& z,
               GetRefG1 const& get_gx, G1 const& gz)
        : tag(tag),
          x(std::move(x)),
          a(std::move(a)),
          z(z),
          get_gx(get_gx),
          gz(gz) {
      Check();
    }

    void SortAndAlign(CommitmentPub& com_pub, CommitmentSec& com_sec) {
      auto order = GetSortOrder(a);
      PermuteAndAlign(order, com_pub);
      PermuteAndAlign(order, com_sec);
      PermuteAndAlign(order, x);
      PermuteAndAlign(order, a);
    }

    void Update(Fr const& alpha, Fr const& beta, Fr const& e, Fr const& ee) {
      Tick tick(__FN__, to_string());
      auto m2 = m() / 2;
      std::vector<std::vector<Fr>> x2(m2);
      std::vector<std::vector<Fr>> a2(m2);
      auto pf = [this, &x2, &a2, &e](int64_t i) {
        x2[i] = x[2 * i + 1] * e + x[2 * i];
        a2[i] = a[2 * i] * e + a[2 * i + 1];
      };
      parallel::For(m2, pf, n() < 1024);
      x2.swap(x);
      a2.swap(a);

      z = alpha * ee + z * e + beta;
    }

   private:
    void Check() {
      CHECK(!x.empty() && x.size() == a.size(), "");
      for (size_t i = 0; i < x.size(); ++i) {
        CHECK(x[i].size() == a[i].size(), "");
        max_n = std::max(max_n, x[i].size());
      }

      if (DEBUG_CHECK) {
        Fr check_z = FrZero();
        for (int64_t i = 0; i < m(); ++i) {
          check_z += InnerProduct(x[i], a[i]);
        }
        CHECK(z == check_z, "");
      }
    }
  };

  struct CommitmentExtPub {
    // 5.3 Recursive
    std::vector<G1> cl;  // size = log(m)
    std::vector<G1> cu;  // size = log(m)
    bool CheckFormat(int64_t check_m) const {
      if (cl.size() != cu.size()) return false;
      return m() == check_m;
    }
    int64_t m() const { return 1LL << cl.size(); }
    bool operator==(CommitmentExtPub const& right) const {
      return cl == right.cl && cu == right.cu;
    }

    bool operator!=(CommitmentExtPub const& right) const {
      return !(*this == right);
    }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("a4.cep", ("cl", cl), ("cu", cu));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("a4.cep", ("cl", cl), ("cu", cu));
    }
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;  // 2*log(m) G1
    A3::Proof proof_a3;            // (2+2log(n))G1 + 2Fr

    int64_t m() const { return com_ext_pub.m(); }

    bool CheckFormat(int64_t check_m) const {
      return com_ext_pub.CheckFormat(check_m) && proof_a3.CheckFormat();
    }
    bool operator==(Proof const& right) const {
      return com_ext_pub == right.com_ext_pub && proof_a3 == right.proof_a3;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("a4.pf", ("c", com_ext_pub), ("r", proof_a3));
    }

    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("a4.pf", ("c", com_ext_pub), ("r", proof_a3));
    }
  };

  static void ComputeCom(ProveInput const& input, CommitmentPub* com_pub,
                         CommitmentSec const& com_sec) {
    Tick tick(__FN__, input.to_string());
    auto const m = input.m();
    // auto const n = input.n();

    com_pub->cx.resize(m);

    auto parallel_f = [&input, &com_pub, &com_sec](int64_t i) {
      com_pub->cx[i] = pc::ComputeCom(input.get_gx, input.x[i], com_sec.r[i]);
    };
    parallel::For(m, parallel_f);

    com_pub->cz = pc::ComputeCom(input.gz, input.z, com_sec.t);
  }

  static void ComputeCom(ProveInput const& input, CommitmentPub* com_pub,
                         CommitmentSec* com_sec) {
    // Tick tick(__FN__);
    auto const m = input.m();
    com_sec->r.resize(m);
    FrRand(com_sec->r.data(), m);

    com_sec->t = FrRand();

    ComputeCom(input, com_pub, *com_sec);
  }

  static void ProveFinal(Proof& proof, h256_t const& seed,
                         ProveInput const& input, CommitmentPub const& com_pub,
                         CommitmentSec const& com_sec) {
    Tick tick(__FN__, input.to_string());
    CHECK(input.m() == 1, "");

    A3::ProveInput input_a3(input.tag, input.x[0], input.a[0], input.z,
                            input.get_gx, input.gz);
    A3::CommitmentPub com_pub_a3(com_pub.cx[0], com_pub.cz);
    A3::CommitmentSec com_sec_a3(com_sec.r[0], com_sec.t);
    A3::Prove(proof.proof_a3, seed, input_a3, com_pub_a3, com_sec_a3);
  }

  static void ComputeSigmaXA(ProveInput const& input, Fr* alpha, Fr* beta) {
    Tick tick(__FN__, input.to_string());
    int64_t m = input.m();
    auto m2 = m / 2;
    std::vector<Fr> xa1(m2, FrZero());
    std::vector<Fr> xa2(m2, FrZero());
    auto parallel_f = [&input, &xa1, &xa2](int64_t i) {
      auto const& x1 = input.x[2 * i + 1];
      auto const& a1 = input.a[2 * i];
      xa1[i] = InnerProduct(x1, a1);

      auto const& x2 = input.x[2 * i];
      auto const& a2 = input.a[2 * i + 1];
      xa2[i] = InnerProduct(x2, a2);
    };
    parallel::For(m2, parallel_f, m2 < 1024);

    *alpha = parallel::Accumulate(xa1.begin(), xa1.end(), FrZero());
    *beta = parallel::Accumulate(xa2.begin(), xa2.end(), FrZero());
  }

  static void UpdateCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                        Fr const& tl, Fr const& tu, G1 const& cl, G1 const& cu,
                        Fr const& e, Fr const& ee) {
    Tick tick(__FN__);
    CommitmentPub com_pub2;
    CommitmentSec com_sec2;
    auto m2 = com_pub.cx.size() / 2;
    com_pub2.cx.resize(m2);
    com_sec2.r.resize(m2);

    auto parallel_f = [&com_pub, &com_sec, &com_pub2, &com_sec2,
                       &e](int64_t i) {
      auto& cx2 = com_pub2.cx;
      auto const& cx = com_pub.cx;
      cx2[i] = cx[2 * i] + cx[2 * i + 1] * e;

      auto& r2 = com_sec2.r;
      auto const& r = com_sec.r;
      r2[i] = r[2 * i] + r[2 * i + 1] * e;
    };
    parallel::For((int64_t)m2, parallel_f, m2 < 1024);

    com_pub2.cz = cl * ee + com_pub.cz * e + cu;
    com_sec2.t = tl * ee + com_sec.t * e + tu;
    com_pub = std::move(com_pub2);
    com_sec = std::move(com_sec2);
  }

  static Fr ComputeChallenge(h256_t const& seed, CommitmentPub const& com_pub,
                             G1 const& cl, G1 const& cu) {
    // Tick tick(__FN__);
    CryptoPP::Keccak_256 hash;
    h256_t digest;
    HashUpdate(hash, seed);
    HashUpdate(hash, cl);
    HashUpdate(hash, cu);
    HashUpdate(hash, com_pub.cx);
    HashUpdate(hash, com_pub.cz);
    hash.Final(digest.data());
    return H256ToFr(digest);
  }

  static void ProveRecursive(Proof& proof, h256_t& seed, ProveInput& input,
                             CommitmentPub& com_pub, CommitmentSec& com_sec) {
    Tick tick(__FN__, input.to_string());
    assert(input.m() > 1);

    Fr alpha, beta;
    ComputeSigmaXA(input, &alpha, &beta);

    // compute cl, cu
    Fr tl = FrRand();
    Fr tu = FrRand();
    G1 cl = pc::ComputeCom(input.gz, alpha, tl);
    G1 cu = pc::ComputeCom(input.gz, beta, tu);
    proof.com_ext_pub.cl.push_back(cl);
    proof.com_ext_pub.cu.push_back(cu);

    // challenge
    Fr e = ComputeChallenge(seed, com_pub, cl, cu);
    Fr ee = e * e;
    seed = FrToBin(e);

    input.Update(alpha, beta, e, ee);

    UpdateCom(com_pub, com_sec, tl, tu, cl, cu, e, ee);

    // debug check com_pub2 and com_sec2
    if (DEBUG_CHECK) {
      CommitmentPub check_com_pub;
      ComputeCom(input, &check_com_pub, com_sec);
      CHECK(check_com_pub.cx == com_pub.cx, "");
      CHECK(check_com_pub.cz == com_pub.cz, "");
    }
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput&& input,
                    CommitmentPub&& com_pub, CommitmentSec&& com_sec) {
    Tick tick(__FN__, input.to_string());

    input.SortAndAlign(com_pub, com_sec);

    while (input.m() > 1) {
      ProveRecursive(proof, seed, input, com_pub, com_sec);
    }
    return ProveFinal(proof, seed, input, com_pub, com_sec);
  }

  static bool Verify(Proof const& proof, h256_t seed, VerifyInput input) {
    Tick tick(__FN__, input.to_string());

    input.SortAndAlign();

    if (!proof.CheckFormat(input.m())) {
      assert(false);
      return false;
    }

    CommitmentPub& com_pub = input.com_pub;

    for (size_t loop = 0; loop < proof.com_ext_pub.cl.size(); ++loop) {
      // challenge
      auto const& cl = proof.com_ext_pub.cl[loop];
      auto const& cu = proof.com_ext_pub.cu[loop];
      Fr e = ComputeChallenge(seed, com_pub, cl, cu);
      Fr ee = e * e;
      seed = FrToBin(e);

      input.Update(e);

      std::vector<G1> cx2(com_pub.m() / 2);
      G1 cz2;

      auto m2 = com_pub.m() / 2;
      auto parallel_f = [&com_pub, &cx2, &e](int64_t i) {
        auto const& cx = com_pub.cx;
        cx2[i] = cx[2 * i] + cx[2 * i + 1] * e;
      };
      parallel::For(m2, parallel_f, m2 < 1024);

      cz2 = cl * ee + com_pub.cz * e + cu;

      com_pub.cx = std::move(cx2);
      com_pub.cz = std::move(cz2);
    }

    assert(com_pub.m() == 1);
    assert(input.a.size() == 1);

    A3::CommitmentPub com_pub_a3(com_pub.cx[0], com_pub.cz);
    A3::VerifyInput verifier_input_a3(input.tag, input.a[0], com_pub_a3,
                                      input.get_gx, input.gz);
    return A3::Verify(proof.proof_a3, seed, verifier_input_a3);
  }

 private:
  static std::vector<size_t> GetSortOrder(std::vector<size_t> const& mn) {
    std::vector<size_t> order(mn.size());
    for (size_t i = 0; i < order.size(); ++i) {
      order[i] = i;
    }

    std::stable_sort(order.begin(), order.end(),
                     [&mn](size_t a, size_t b) { return mn[a] > mn[b]; });

    return order;
  }

  static std::vector<size_t> GetSortOrder(
      std::vector<std::vector<Fr>> const& data) {
    std::vector<size_t> mn(data.size());
    for (size_t i = 0; i < mn.size(); ++i) {
      mn[i] = data[i].size();
    }
    return GetSortOrder(mn);
  }

  template <typename T>
  static void Permute(std::vector<size_t> const& order, std::vector<T>& v) {
    CHECK(order.size() == v.size(), "");

    std::vector<T> v2(v.size());
    for (size_t i = 0; i < order.size(); ++i) {
      v2[i] = std::move(v[order[i]]);
    }
    v.swap(v2);
  }

  static void PermuteAndAlign(std::vector<size_t> const& order,
                              CommitmentPub& v) {
    Permute(order, v.cx);
    int64_t old_m = v.cx.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      v.cx.resize(new_m, G1Zero());
    }
  }

  static void PermuteAndAlign(std::vector<size_t> const& order,
                              CommitmentSec& v) {
    Permute(order, v.r);
    int64_t old_m = v.r.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      v.r.resize(new_m, FrZero());
    }
  }

  static void PermuteAndAlign(std::vector<size_t> const& order,
                              std::vector<std::vector<Fr>>& v) {
    Permute(order, v);
    int64_t old_m = v.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      v.resize(new_m);
    }
  }

 public:
  static bool Test(int64_t m, int64_t n) {
    Tick tick(__FN__);
    std::cout << "m=" << m << ", n=" << n << "\n";

    std::vector<std::vector<Fr>> x(m);
    for (auto& i : x) {
      i.resize(n);
      FrRand(i.data(), n);
    }
    x.back().resize(n / 2);

    std::vector<std::vector<Fr>> a(m);
    for (auto& i : a) {
      i.resize(n);
      FrRand(i.data(), n);
    }
    a.back().resize(n / 2);

    Fr z = FrZero();
    for (int64_t i = 0; i < m; ++i) {
      z += InnerProduct(x[i], a[i]);
    }

    h256_t seed = misc::RandH256();

    int64_t x_g_offset = 20;
    GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
      return pc::PcG()[x_g_offset + i];
    };

    auto a_copy = a;
    ProveInput prove_input("test", std::move(x), std::move(a), z, get_gx,
                           pc::PcU());
    CommitmentPub com_pub;
    CommitmentSec com_sec;
    ComputeCom(prove_input, &com_pub, &com_sec);

    Proof proof;
    auto copy_com_pub = com_pub;
    Prove(proof, seed, std::move(prove_input), std::move(com_pub),
          std::move(com_sec));

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

    VerifyInput verify_input("test", std::move(copy_com_pub), get_gx,
                             std::move(a_copy), pc::PcU());
    bool success = Verify(proof, seed, verify_input);
    std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
    return success;
  }
};
}  // namespace hyrax
