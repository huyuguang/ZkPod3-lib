#pragma once

#include <libsnark/gadgetlib1/protoboard.hpp>

#include "./details.h"
#include "groth09/groth09.h"

// pb: r1cs of one circuit. m=pb.num_constraints(),s=pb.num_variables()
// w: matrix<Fr, s, n>, all var(witness) of the circuits.
// com_w_r: array<Fr, s>, random fr.
// com_w: array<Fr, s>, com_w[i] = com(w[i], com_w_r[i]).
// for open n instance of circuit, open com_w, prove consistency of the com_w

namespace clink {

template <typename Policy>
struct ParallelR1cs {       
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using Sec43 = typename Policy::Sec43;
  using Proof = typename Sec43::Proof;

  struct ProveInput {
    ProveInput(libsnark::protoboard<Fr> const& pb,
                std::vector<std::vector<Fr>>&& w, std::vector<G1> const& com_w,
                std::vector<Fr> const& com_w_r, int64_t g_offset)
        : pb(pb),
          com_w(com_w),
          com_w_r(com_w_r),
          g_offset(g_offset),
          m((int64_t)pb.num_constraints()),
          s((int64_t)pb.num_variables()),
          n((int64_t)w[0].size()),
          constraint_system(pb.get_constraint_system()),
          constraints(constraint_system.constraints) {
#ifdef _DEBUG
      assert((int64_t)w.size() == s);
      assert((int64_t)com_w.size() == s);
      assert((int64_t)com_w_r.size() == s);
      for (int64_t i = 0; i < s; ++i) {
        assert(com_w[i] == PcComputeCommitmentG(g_offset, w[i], com_w_r[i]));
      }
      for (size_t i = 0; i < constraint_system.primary_input_size; ++i) {
        assert(com_w_r[i] == 0);
      }
#endif

      x.resize(m);
      for (auto& i : x) i.resize(n);

      y.resize(m);
      for (auto& i : y) i.resize(n);

      z.resize(m);
      for (auto& i : z) i.resize(n);

      auto parallel_f = [this, &w](int64_t j) {
        std::vector<Fr> witness(s);
        for (int64_t i = 0; i < s; ++i) {
          witness[i] = w[i][j];
        }
        for (int64_t i = 0; i < m; ++i) {
          auto const& constraint = constraints[i];
          x[i][j] = constraint.a.evaluate(witness);
          y[i][j] = constraint.b.evaluate(witness);
          z[i][j] = constraint.c.evaluate(witness);
          assert(z[i][j] == x[i][j] * y[i][j]);
        }
      };
      parallel::For(n, parallel_f);
    }

    libsnark::protoboard<Fr> const& pb;
    std::vector<G1> const& com_w;
    std::vector<Fr> const& com_w_r;
    int64_t const g_offset;
    int64_t const m;
    int64_t const s;
    int64_t const n;
    libsnark::r1cs_constraint_system<Fr> constraint_system;
    std::vector<libsnark::r1cs_constraint<Fr>> const& constraints;
    std::vector<std::vector<Fr>> x;
    std::vector<std::vector<Fr>> y;
    std::vector<std::vector<Fr>> z;
  };

  // w: s*n
  static void Prove(Proof& proof, h256_t seed, ProveInput&& input) {
    Tick tick(__FUNCTION__);
    int64_t m = input.m;
    int64_t n = input.n;

    UpdateSeed(seed, input.com_w);

    // prove hadamard product
    typename Sec43::ProveInput input_43(m, std::move(input.x), std::move(input.y),
                                 std::move(input.z), input.g_offset,
                                 input.g_offset, input.g_offset);

    typename Sec43::CommitmentPub com_pub;
    typename Sec43::CommitmentSec com_sec;
    BuildHpCom(m, n, input.com_w, input.com_w_r, input.constraints,
                        input.g_offset, com_pub, com_sec);
    DebugCheckHpCom(m, input_43, com_pub, com_sec);

    Sec43::AlignData(input_43, com_pub, com_sec);
    Sec43::Prove(proof, seed, std::move(input_43), std::move(com_pub),
                  std::move(com_sec));
  }

  struct VerifyInput {
    VerifyInput(int64_t n, libsnark::protoboard<Fr> const& pb,
                  std::vector<G1> const& com_w,
                  std::vector<std::vector<Fr>> const& public_w,
                  int64_t g_offset)
        : n(n),
          pb(pb),
          com_w(com_w),
          public_w(public_w),
          g_offset(g_offset),
          m((int64_t)pb.num_constraints()),
          s((int64_t)pb.num_variables()),
          constraint_system(pb.get_constraint_system()),
          constraints(constraint_system.constraints) {}

    bool Check() const {
      if ((int64_t)com_w.size() != s) {
        assert(false);
        return false;
      }

      for (auto const& i : public_w) {
        if ((int64_t)i.size() != n) {
          assert(false);
          return false;
        }
      }
      // check public_w
      if (public_w.size() != constraint_system.primary_input_size) {
        assert(false);
        return false;
      }

      bool all_success = false;
      auto parallel_f = [this](int64_t i) {
        return com_w[i] ==
               PcComputeCommitmentG(g_offset, public_w[i], FrZero());
      };
      parallel::For(&all_success, (int64_t)public_w.size(), parallel_f);
      if (!all_success) {
        std::cerr << "ASSERT: " << __FUNCTION__ << ": " << __LINE__ << "\n";
        return false;
      }
      return true;
    }

    int64_t const n;
    libsnark::protoboard<Fr> const& pb;
    std::vector<G1> const& com_w;
    std::vector<std::vector<Fr>> const& public_w;
    int64_t const g_offset;
    int64_t const m;
    int64_t const s;
    libsnark::r1cs_constraint_system<Fr> constraint_system;
    std::vector<libsnark::r1cs_constraint<Fr>> const& constraints;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FUNCTION__);

    if (!input.Check()) {
      assert(false);
      return false;
    }

    UpdateSeed(seed, input.com_w);

    typename Sec43::CommitmentPub com_pub;
    BuildHpCom(input.m, input.n, input.com_w, input.constraints,
                        input.g_offset, com_pub);
    com_pub.Align();
    typename Sec43::VerifyInput input_43(input.m, input.n, com_pub, input.g_offset,
                                   input.g_offset, input.g_offset);
    return Sec43::Verify(proof, seed, input_43);
  }

  // see match.h
  static bool Test() { return true; }

 private:
  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  static void BuildHpCom(std::vector<G1> const& com_w,
                         std::vector<Fr> const& com_w_r,
                         libsnark::linear_combination<Fr> const& lc,
                         G1& com_pub, Fr& com_sec, G1 const& sigma_g) {
    // Tick tick(__FUNCTION__);
    for (auto const& term : lc.terms) {
      if (term.index == 0) {  // constants
        com_pub += sigma_g * term.coeff;
      } else {
        com_pub += com_w[term.index - 1] * term.coeff;
        com_sec += com_w_r[term.index - 1] * term.coeff;
      }
    }
  }

  static void BuildHpCom(
      int64_t m, int64_t n, std::vector<G1> const& com_w,
      std::vector<Fr> const& com_w_r,
      std::vector<libsnark::r1cs_constraint<Fr>> const& constraints,
      int64_t g_offset, typename Sec43::CommitmentPub& com_pub,
      typename Sec43::CommitmentSec& com_sec) {
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);
    com_sec.r.resize(m);
    FrZero(com_sec.r);
    com_sec.s.resize(m);
    FrZero(com_sec.s);
    com_sec.t.resize(m);
    FrZero(com_sec.t);

    auto pds_sigma_g = PcComputeSigmaG(g_offset, n);
    auto parallel_f = [&com_pub, &com_sec, &com_w, &com_w_r, &pds_sigma_g,
                       &constraints](int64_t i) {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto& com_sec_r = com_sec.r[i];
      auto& com_sec_s = com_sec.s[i];
      auto& com_sec_t = com_sec.t[i];
      auto const& constraint = constraints[i];
      BuildHpCom(com_w, com_w_r, constraint.a, com_pub_a, com_sec_r,
                 pds_sigma_g);
      BuildHpCom(com_w, com_w_r, constraint.b, com_pub_b, com_sec_s,
                 pds_sigma_g);
      BuildHpCom(com_w, com_w_r, constraint.c, com_pub_c, com_sec_t,
                 pds_sigma_g);
    };
    parallel::For(m, parallel_f);
  }

  static void DebugCheckHpCom(int64_t m, typename Sec43::ProveInput const& input,
                              typename Sec43::CommitmentPub const& com_pub,
                              typename Sec43::CommitmentSec const& com_sec) {
#ifdef _DEBUG
    Tick tick(__FUNCTION__);
    for (int64_t i = 0; i < m; ++i) {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto& com_sec_r = com_sec.r[i];
      auto& com_sec_s = com_sec.s[i];
      auto& com_sec_t = com_sec.t[i];

      auto const& xi = input.x[i];
      G1 check_com_pub_a =
          PcComputeCommitmentG(input.x_g_offset, xi, com_sec_r);
      assert(check_com_pub_a == com_pub_a);

      auto const& yi = input.y[i];
      G1 check_com_pub_b =
          PcComputeCommitmentG(input.y_g_offset, yi, com_sec_s);
      assert(check_com_pub_b == com_pub_b);

      auto const& zi = input.z[i];
      G1 check_com_pub_c =
          PcComputeCommitmentG(input.z_g_offset, zi, com_sec_t);
      assert(check_com_pub_c == com_pub_c);
    }
#else
    (void)m;
    (void)input;
    (void)com_pub;
    (void)com_sec;
#endif
  }

  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  static void BuildHpCom(std::vector<G1> const& com_w,
                         libsnark::linear_combination<Fr> const& lc,
                         G1& com_pub, G1 const& sigma_g) {
    // Tick tick(__FUNCTION__);
    for (auto const& term : lc.terms) {
      if (term.index == 0) {  // constants
        com_pub += sigma_g * term.coeff;
      } else {
        com_pub += com_w[term.index - 1] * term.coeff;
      }
    }
  }

  static void BuildHpCom(
      int64_t m, int64_t n, std::vector<G1> const& com_w,
      std::vector<libsnark::r1cs_constraint<Fr>> const& constraints,
      int64_t g_offset, typename Sec43::CommitmentPub& com_pub) {
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);

    auto pds_sigma_g = PcComputeSigmaG(g_offset, n);
    auto parallel_f = [&com_pub, &com_w, &pds_sigma_g,
                       &constraints](int64_t i) {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto const& constraint = constraints[i];
      BuildHpCom(com_w, constraint.a, com_pub_a, pds_sigma_g);
      BuildHpCom(com_w, constraint.b, com_pub_b, pds_sigma_g);
      BuildHpCom(com_w, constraint.c, com_pub_c, pds_sigma_g);
    };
    parallel::For(m, parallel_f);
  }

  static void UpdateSeed(h256_t& seed, std::vector<G1> const& com_w) {
    // update seed
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_w);
    hash.Final(seed.data());
  }
};

}  // namespace clink