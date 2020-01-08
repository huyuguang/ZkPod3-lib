#pragma once

#include <libsnark/gadgetlib1/protoboard.hpp>

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"
#include "groth09/groth09.h"

// pb: r1cs of one circuit. m=pb.num_constraints(),s=pb.num_variables()
// w: matrix<Fr, s, n>, all var(witness) of the circuits.
// com_w_r: array<Fr, s>, random fr.
// com_w: array<Fr, s>, com_w[i] = com(w[i], com_w_r[i]).
// for open n instance of circuit, open com_w, prove consistency of the com_w

namespace parallel_r1cs {

namespace details {

// com(<A,X>) or com(<B,X>) or com(<C,X>)
void BuildHpCom(std::vector<G1> const& com_w, std::vector<Fr> const& com_w_r,
                libsnark::linear_combination<Fr> const& lc, G1& com_pub,
                Fr& com_sec, G1 const& sigma_g) {
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

void BuildHpCom(int64_t m, int64_t n, std::vector<G1> const& com_w,
                std::vector<Fr> const& com_w_r,
                std::vector<libsnark::r1cs_constraint<Fr>> const& constraints,
                groth09::sec43::CommitmentPub& com_pub,
                groth09::sec43::CommitmentSec& com_sec) {
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

  auto pds_sigma_g = PcComputeSigmaG(n);
  auto parallel_f = [&com_pub, &com_sec, &com_w, &com_w_r, &pds_sigma_g,
                     &constraints](int64_t i) {
    auto& com_pub_a = com_pub.a[i];
    auto& com_pub_b = com_pub.b[i];
    auto& com_pub_c = com_pub.c[i];
    auto& com_sec_r = com_sec.r[i];
    auto& com_sec_s = com_sec.s[i];
    auto& com_sec_t = com_sec.t[i];
    auto const& constraint = constraints[i];
    BuildHpCom(com_w, com_w_r, constraint.a, com_pub_a, com_sec_r, pds_sigma_g);
    BuildHpCom(com_w, com_w_r, constraint.b, com_pub_b, com_sec_s, pds_sigma_g);
    BuildHpCom(com_w, com_w_r, constraint.c, com_pub_c, com_sec_t, pds_sigma_g);
  };
  parallel::For(m, parallel_f);
}

void DebugCheckHpCom(int64_t m, groth09::sec43::ProverInput const& input,
                     groth09::sec43::CommitmentPub const& com_pub,
                     groth09::sec43::CommitmentSec const& com_sec) {
#ifdef _DEBUG
  Tick tick(__FUNCTION__);
  for (int64_t i = 0; i < m; ++i) {
    auto& com_pub_a = com_pub.a[i];
    auto& com_pub_b = com_pub.b[i];
    auto& com_pub_c = com_pub.c[i];
    auto& com_sec_r = com_sec.r[i];
    auto& com_sec_s = com_sec.s[i];
    auto& com_sec_t = com_sec.t[i];

    auto const& xi = input.x(i);
    G1 check_com_pub_a = PcComputeCommitment(xi, com_sec_r);
    assert(check_com_pub_a == com_pub_a);

    auto const& yi = input.y(i);
    G1 check_com_pub_b = PcComputeCommitment(yi, com_sec_s);
    assert(check_com_pub_b == com_pub_b);

    auto const& zi = input.z(i);
    G1 check_com_pub_c = PcComputeCommitment(zi, com_sec_t);
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
void BuildHpCom(std::vector<G1> const& com_w,
                libsnark::linear_combination<Fr> const& lc, G1& com_pub,
                G1 const& sigma_g) {
  // Tick tick(__FUNCTION__);
  for (auto const& term : lc.terms) {
    if (term.index == 0) {  // constants
      com_pub += sigma_g * term.coeff;
    } else {
      com_pub += com_w[term.index - 1] * term.coeff;
    }
  }
}

void BuildHpCom(int64_t m, int64_t n, std::vector<G1> const& com_w,
                std::vector<libsnark::r1cs_constraint<Fr>> const& constraints,
                groth09::sec43::CommitmentPub& com_pub) {
  com_pub.a.resize(m);
  G1Zero(com_pub.a);
  com_pub.b.resize(m);
  G1Zero(com_pub.b);
  com_pub.c.resize(m);
  G1Zero(com_pub.c);

  auto pds_sigma_g = PcComputeSigmaG(n);
  auto parallel_f = [&com_pub, &com_w, &pds_sigma_g, &constraints](int64_t i) {
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

inline void UpdateSeed(h256_t& seed, std::vector<G1> const& com_w) {
  // update seed
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_w);
  hash.Final(seed.data());
}

}  // namespace details

// w: s*n
inline void Prove(groth09::sec43::RomProof& proof, h256_t seed,
                  libsnark::protoboard<Fr> const& pb,
                  std::vector<std::vector<Fr>> w, std::vector<G1> const& com_w,
                  std::vector<Fr> const& com_w_r) {
  Tick tick(__FUNCTION__);
  int64_t m = (int64_t)pb.num_constraints();
  int64_t s = (int64_t)pb.num_variables();
  assert((int64_t)w.size() == s);
  assert((int64_t)com_w.size() == s);
  assert((int64_t)com_w_r.size() == s);
  int64_t n = (int64_t)w[0].size();
  auto constraint_system = pb.get_constraint_system();
  auto const& constraints = constraint_system.constraints;

#ifdef _DEBUG
  for (int64_t i = 0; i < s; ++i) {
    assert(com_w[i] == PcComputeCommitment(w[i], com_w_r[i]));
  }
  for (auto i = 0; i < constraint_system.primary_input_size; ++i) {
    assert(com_w_r[i] == 0);
  }
#endif

  details::UpdateSeed(seed, com_w);

  // build hadamard product x, y, z, has x o y = z
  std::vector<std::vector<Fr>> x(m);
  std::vector<std::vector<Fr>> y(m);
  std::vector<std::vector<Fr>> z(m);
  for (auto& i : x) i.resize(n);
  for (auto& i : y) i.resize(n);
  for (auto& i : z) i.resize(n);

  auto parallel_f = [&constraints, &w, &x, &y, &z, m, s](int64_t j) {
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

  // now we do not need w
  w.clear();
  w.shrink_to_fit();

  // prove hadamard product
  groth09::sec43::ProverInput input(std::move(x), std::move(y), std::move(z));
  groth09::sec43::CommitmentPub com_pub;
  groth09::sec43::CommitmentSec com_sec;
  details::BuildHpCom(m, n, com_w, com_w_r, constraints, com_pub, com_sec);
  details::DebugCheckHpCom(m, input, com_pub, com_sec);

  groth09::sec43::AlignData(input, com_pub, com_sec);
  groth09::sec43::RomProve(proof, seed, std::move(input), std::move(com_pub),
                           std::move(com_sec));
}

inline bool Verify(groth09::sec43::RomProof const& proof, h256_t seed,
                   int64_t n, libsnark::protoboard<Fr> const& pb,
                   std::vector<G1> const& com_w,
                   std::vector<std::vector<Fr>> const& public_w) {
  Tick tick(__FUNCTION__);
  auto constraint_system = pb.get_constraint_system();
  auto const& constraints = constraint_system.constraints;
  int64_t m = (int64_t)pb.num_constraints();
  int64_t s = (int64_t)pb.num_variables();
  if ((int64_t)com_w.size() != s) {
    assert(false);
    return false;
  }
  if ((int64_t)misc::Pow2UB(m) != proof.m()) {
    assert(false);
    return false;
  }

  for (auto const& i : public_w) {
    if ((int64_t)i.size() != n) {
      assert(false);
      return false;
    }
  }

  details::UpdateSeed(seed, com_w);

  // check public_w
  if (public_w.size() != constraint_system.primary_input_size) {
    assert(false);
    return false;
  }

  bool all_success = false;
  auto parallel_f = [&com_w, &public_w](int64_t i) {
    return com_w[i] == PcComputeCommitment(public_w[i], FrZero());
  };
  parallel::For(&all_success, (int64_t)public_w.size(), parallel_f);
  if (!all_success) {
    std::cerr << "ASSERT: " << __FUNCTION__ << ": " << __LINE__ << "\n";
    return false;
  }

  groth09::sec43::CommitmentPub com_pub;
  details::BuildHpCom(m, n, com_w, constraints, com_pub);
  com_pub.Align();
  groth09::sec43::VerifierInput input(com_pub);
  return groth09::sec43::RomVerify(proof, seed, input);
}

// see match.h
inline bool Test() { return true; }
}  // namespace parallel_r1cs