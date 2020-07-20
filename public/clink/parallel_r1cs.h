#pragma once

#include <libsnark/gadgetlib1/protoboard.hpp>

#include "./details.h"
#include "groth09/groth09.h"

// r1cs_info: r1cs of one circuit. m=r1cs_info.num_constraints()
// s=r1cs_info.num_variables()
// w: matrix<Fr, s, n>, all var(witness) of the circuits.
// com_w_r: array<Fr, s>, random fr.
// com_w: array<Fr, s>, com_w[i] = com(w[i], com_w_r[i]).
// for open n instance of circuit, open com_w, prove consistency of the com_w

namespace clink {

struct R1csInfo {
  R1csInfo(libsnark::protoboard<Fr> const &pb)
      : num_constraints(pb.num_constraints()),
        num_variables(pb.num_variables()),
        constraint_system(pb.get_constraint_system()) {}
  int64_t num_constraints;
  int64_t num_variables;
  libsnark::r1cs_constraint_system<Fr> constraint_system;
  std::string to_string() const {
    return "R1csInfo: num_constraints: " + std::to_string(num_constraints) +
           ", num_variables: " + std::to_string(num_variables);
  }
};

template <typename Policy>
struct ParallelR1cs {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using Sec43 = typename Policy::Sec43;
  using Proof = typename Sec43::Proof;

  struct ProveInput {
    ProveInput(R1csInfo const &r1cs_info, std::string const &unique_tag,
               std::vector<std::vector<Fr>> &&w, std::vector<G1> const &com_w,
               std::vector<Fr> const &com_w_r, pc::GetRefG const &get_g)
        : r1cs_info(r1cs_info),
          unique_tag(unique_tag),
          com_w(com_w),
          com_w_r(com_w_r),
          get_g(get_g),
          m(r1cs_info.num_constraints),
          s(r1cs_info.num_variables),
          n((int64_t)w[0].size()) {
#ifdef _DEBUG
      if ((int64_t)w.size() != s) throw std::runtime_error("opps");
      if ((int64_t)com_w.size() != s) throw std::runtime_error("opps");
      if ((int64_t)com_w_r.size() != s) throw std::runtime_error("opps");
      for (size_t i = 0; i < constraint_system().primary_input_size; ++i) {
        if (com_w_r[i] != 0) throw std::runtime_error("opps");
      }

      for (int64_t i = 0; i < s; ++i) {
        if (com_w[i] != pc::PcComputeCommitmentG(get_g, w[i], com_w_r[i]))
          throw std::runtime_error("opps");
      }
#endif

      x.resize(m);
      for (auto &i : x) i.resize(n);

      y.resize(m);
      for (auto &i : y) i.resize(n);

      z.resize(m);
      for (auto &i : z) i.resize(n);

      auto parallel_f = [this, &w](int64_t j) {
        std::vector<Fr> witness(s);
        for (int64_t i = 0; i < s; ++i) {
          witness[i] = w[i][j];
        }
        for (int64_t i = 0; i < m; ++i) {
          auto const &constraint = constraints()[i];
          x[i][j] = constraint.a.evaluate(witness);
          y[i][j] = constraint.b.evaluate(witness);
          z[i][j] = constraint.c.evaluate(witness);
          assert(z[i][j] == x[i][j] * y[i][j]);
        }
      };
      parallel::For(n, parallel_f);
    }

    libsnark::r1cs_constraint_system<Fr> const &constraint_system() const {
      return r1cs_info.constraint_system;
    }

    std::vector<libsnark::r1cs_constraint<Fr>> const &constraints() const {
      return constraint_system().constraints;
    }

    R1csInfo const &r1cs_info;
    std::string unique_tag;
    std::vector<G1> const &com_w;
    std::vector<Fr> const &com_w_r;
    pc::GetRefG const &get_g;
    int64_t const m;
    int64_t const s;
    int64_t const n;
    std::vector<std::vector<Fr>> x;
    std::vector<std::vector<Fr>> y;
    std::vector<std::vector<Fr>> z;
  };

  // w: s*n
  static void Prove(Proof &proof, h256_t seed, ProveInput &&input) {
    Tick tick(__FN__);
    int64_t m = input.m;
    int64_t n = input.n;

    UpdateSeed(seed, input.com_w);

    // prove hadamard product
    typename Sec43::ProveInput input_43(m, std::move(input.x),
                                        std::move(input.y), std::move(input.z),
                                        input.get_g, input.get_g, input.get_g);

    typename Sec43::CommitmentPub com_pub;
    typename Sec43::CommitmentSec com_sec;
    BuildHpCom(m, n, input.com_w, input.com_w_r, input.constraints(),
               input.get_g, com_pub, com_sec);
    DebugCheckHpCom(m, input_43, com_pub, com_sec);

    Sec43::AlignData(input_43, com_pub, com_sec);
    Sec43::Prove(proof, seed, std::move(input_43), std::move(com_pub),
                 std::move(com_sec));
  }

  struct VerifyInput {
    VerifyInput(int64_t n, R1csInfo const &r1cs_info,
                std::string const &unique_tag,
                std::vector<G1> const &com_w,
                std::vector<std::vector<Fr>> const &public_w,
                pc::GetRefG const &get_g)
        : n(n),
          r1cs_info(r1cs_info),
          unique_tag(unique_tag),
          com_w(com_w),
          public_w(public_w),
          get_g(get_g),
          m(r1cs_info.num_constraints),
          s(r1cs_info.num_variables) {}

    libsnark::r1cs_constraint_system<Fr> const &constraint_system() const {
      return r1cs_info.constraint_system;
    }

    std::vector<libsnark::r1cs_constraint<Fr>> const &constraints() const {
      return constraint_system().constraints;
    }

    bool Check() const {
      if ((int64_t)com_w.size() != s) {
        assert(false);
        std::cerr << __FN__ << ":" << __LINE__ << " oops\n";
        return false;
      }

      for (auto const &i : public_w) {
        if ((int64_t)i.size() != n) {
          assert(false);
          std::cout << "public_w.size(): " << public_w.size() << "\n";
          std::cout << "constraint_system().primary_input_size: " <<
              constraint_system().primary_input_size << "\n";
          std::cerr << __FN__ << ":" << __LINE__ << " oops\n";
          return false;
        }
      }
      // check public_w
      if (public_w.size() != constraint_system().primary_input_size) {
        assert(false);
        std::cerr << __FN__ << ":" << __LINE__ << " oops\n";
        return false;
      }

      bool all_success = false;
      auto parallel_f = [this](int64_t i) {
        return com_w[i] ==
               pc::PcComputeCommitmentG(get_g, public_w[i], FrZero());
      };
      parallel::For(&all_success, (int64_t)public_w.size(), parallel_f);
      if (!all_success) {
        std::cerr << __FN__ << ":" << __LINE__ << " oops\n";
        return false;
      }
      return true;
    }

    int64_t const n;
    R1csInfo const &r1cs_info;
    std::string const unique_tag;
    std::vector<G1> const &com_w;
    std::vector<std::vector<Fr>> const &public_w;
    pc::GetRefG const &get_g;
    int64_t const m;
    int64_t const s;
  };

  static bool Verify(Proof const &proof, h256_t seed,
                     VerifyInput const &input) {
    Tick tick(__FN__);

    if (!input.Check()) {
      assert(false);
      return false;
    }

    UpdateSeed(seed, input.com_w);

    typename Sec43::CommitmentPub com_pub;
    BuildHpCom(input.m, input.n, input.com_w, input.constraints(), input.get_g,
               com_pub);
    com_pub.Align();
    typename Sec43::VerifyInput input_43(input.m, input.n, com_pub, input.get_g,
                                         input.get_g, input.get_g);
    return Sec43::Verify(proof, seed, input_43);
  }

  // see match.h
  static bool Test() { return true; }

 public:
  // for prove
  static void BuildHpCom(
      int64_t m, int64_t n, std::vector<G1> const &com_w,
      std::vector<Fr> const &com_w_r,
      std::vector<libsnark::r1cs_constraint<Fr>> const &constraints,
      pc::GetRefG const &get_g, typename Sec43::CommitmentPub &com_pub,
      typename Sec43::CommitmentSec &com_sec) {
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

    auto pds_sigma_g = pc::PcComputeSigmaG(get_g, n);
    auto parallel_f = [&com_pub, &com_sec, &com_w, &com_w_r, &pds_sigma_g,
                       &constraints](int64_t i) {
      auto &com_pub_a = com_pub.a[i];
      auto &com_pub_b = com_pub.b[i];
      auto &com_pub_c = com_pub.c[i];
      auto &com_sec_r = com_sec.r[i];
      auto &com_sec_s = com_sec.s[i];
      auto &com_sec_t = com_sec.t[i];
      auto const &constraint = constraints[i];
      BuildHpCom(com_w, com_w_r, constraint.a, com_pub_a, com_sec_r,
                 pds_sigma_g);
      BuildHpCom(com_w, com_w_r, constraint.b, com_pub_b, com_sec_s,
                 pds_sigma_g);
      BuildHpCom(com_w, com_w_r, constraint.c, com_pub_c, com_sec_t,
                 pds_sigma_g);
    };
    parallel::For(m, parallel_f);
  }

  static void DebugCheckHpCom(int64_t m,
                              typename Sec43::ProveInput const &input,
                              typename Sec43::CommitmentPub const &com_pub,
                              typename Sec43::CommitmentSec const &com_sec) {
#ifdef _DEBUG
    Tick tick(__FN__);
    for (int64_t i = 0; i < m; ++i) {
      auto &com_pub_a = com_pub.a[i];
      auto &com_pub_b = com_pub.b[i];
      auto &com_pub_c = com_pub.c[i];
      auto &com_sec_r = com_sec.r[i];
      auto &com_sec_s = com_sec.s[i];
      auto &com_sec_t = com_sec.t[i];

      auto const &xi = input.x[i];
      G1 check_com_pub_a =
          pc::PcComputeCommitmentG(input.get_gx, xi, com_sec_r);
      assert(check_com_pub_a == com_pub_a);

      auto const &yi = input.y[i];
      G1 check_com_pub_b =
          pc::PcComputeCommitmentG(input.get_gy, yi, com_sec_s);
      assert(check_com_pub_b == com_pub_b);

      auto const &zi = input.z[i];
      G1 check_com_pub_c =
          pc::PcComputeCommitmentG(input.get_gz, zi, com_sec_t);
      assert(check_com_pub_c == com_pub_c);
    }
#else
    (void)m;
    (void)input;
    (void)com_pub;
    (void)com_sec;
#endif
  }

  // for verify
  static void BuildHpCom(
      int64_t m, int64_t n, std::vector<G1> const &com_w,
      std::vector<libsnark::r1cs_constraint<Fr>> const &constraints,
      pc::GetRefG const &get_g, typename Sec43::CommitmentPub &com_pub) {
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);

    auto pds_sigma_g = pc::PcComputeSigmaG(get_g, n);
    auto parallel_f = [&com_pub, &com_w, &pds_sigma_g,
                       &constraints](int64_t i) {
      auto &com_pub_a = com_pub.a[i];
      auto &com_pub_b = com_pub.b[i];
      auto &com_pub_c = com_pub.c[i];
      auto const &constraint = constraints[i];
      BuildHpCom(com_w, constraint.a, com_pub_a, pds_sigma_g);
      BuildHpCom(com_w, constraint.b, com_pub_b, pds_sigma_g);
      BuildHpCom(com_w, constraint.c, com_pub_c, pds_sigma_g);
    };
    parallel::For(m, parallel_f);
  }

 private:
  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  static void BuildHpCom(std::vector<G1> const &com_w,
                         std::vector<Fr> const &com_w_r,
                         libsnark::linear_combination<Fr> const &lc,
                         G1 &com_pub, Fr &com_sec, G1 const &sigma_g) {
    // Tick tick(__FN__);
    for (auto const &term : lc.terms) {
      if (term.coeff == FrZero()) continue;
      if (term.coeff == FrOne()) {
        if (term.index == 0) {  // constants
          com_pub += sigma_g;
        } else {
          com_pub += com_w[term.index - 1];
          com_sec += com_w_r[term.index - 1];
        }
      } else {
        if (term.index == 0) {  // constants
          com_pub += sigma_g * term.coeff;
        } else {
          com_pub += com_w[term.index - 1] * term.coeff;
          com_sec += com_w_r[term.index - 1] * term.coeff;
        }
      }
    }
  }

  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  static void BuildHpCom(std::vector<G1> const &com_w,
                         libsnark::linear_combination<Fr> const &lc,
                         G1 &com_pub, G1 const &sigma_g) {
    // Tick tick(__FN__);
    for (auto const &term : lc.terms) {
      if (term.coeff == FrZero()) continue;
      if (term.coeff == FrOne()) {
        if (term.index == 0) {  // constants
          com_pub += sigma_g;
        } else {
          com_pub += com_w[term.index - 1];
        }
      } else {
        if (term.index == 0) {  // constants
          com_pub += sigma_g * term.coeff;
        } else {
          com_pub += com_w[term.index - 1] * term.coeff;
        }
      }
    }
  }



  static void UpdateSeed(h256_t &seed, std::vector<G1> const &com_w) {
    // update seed
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_w);
    hash.Final(seed.data());
  }
};

}  // namespace clink