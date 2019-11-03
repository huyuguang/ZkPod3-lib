#pragma once

#include "./details.h"

// a is a public vector<Fr>
// x is a secret vector<Fr>
// y = <x,a>
// xi = com(a), delta = com(y)
// open xi, delta, prove y=<x,a>

namespace hyrax::a2 {

struct ProverInput {
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& a) : x(x), a(a) {
    assert(x.size() == a.size() && !a.empty());
    y = InnerProduct(x, a);
  }
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& a, Fr const& y)
      : x(x), a(a), y(y) {
    assert(y == InnerProduct(x, a));
  }
  int64_t n() const { return (int64_t)x.size(); }
  std::vector<Fr> const& x;  // x.size = n
  std::vector<Fr> const& a;  // a.size = n
  Fr y;                      // y = <x, a>
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& xi, G1 const& tau) : xi(xi), tau(tau) {}
  G1 xi;   // com(x,r_xi)
  G1 tau;  // com(y,r_tau)
};

struct CommitmentSec {
  Fr r_xi;
  Fr r_tau;
};

struct CommitmentExtPub {
  CommitmentExtPub() {}
  CommitmentExtPub(G1 const& delta, G1 const& beta)
      : delta(delta), beta(beta) {}
  G1 delta;  // com(d, r_delta)
  G1 beta;   // com(<a,d>, r_beta)
};

struct CommitmentExtSec {
  std::vector<Fr> d;  // size = n
  Fr r_beta;
  Fr r_delta;
};

struct Proof {
  Proof() {}
  Proof(std::shared_ptr<std::vector<Fr>> const& z, Fr const& z_delta,
        Fr const& z_beta)
      : z(z), z_delta(z_delta), z_beta(z_beta) {}
  std::shared_ptr<std::vector<Fr>> z;  // z.size = n
  Fr z_delta;
  Fr z_beta;
  int64_t n() const { return (int64_t)z->size(); }
};

struct RomProof {
  CommitmentExtPub com_ext_pub;
  Proof proof;
  int64_t n() const { return proof.n(); }
};

struct VerifierInput {
  VerifierInput(std::vector<Fr> const& a, CommitmentPub const& com_pub)
      : a(a), com_pub(com_pub) {}
  std::vector<Fr> a;  // a.size = n
  CommitmentPub const& com_pub;
};

// com(n) + com(1) + ip(n)
inline bool VerifyInternal(PdsBase const& pds_base,
                           VerifierInput const& input, Fr const& challenge,
                           CommitmentExtPub const& com_ext_pub,
                           Proof const& proof) {
  Tick tick(__FUNCTION__);
  using details::ComputeCommitment;

  std::cout << Tick::GetIndentString() << "multiexp(" << proof.n() << ")\n";

  auto const& com_pub = input.com_pub;

  G1 left, right;
  auto const& xi = com_pub.xi;
  auto const& delta = com_ext_pub.delta;
  left = xi * challenge + delta;
  right = ComputeCommitment(pds_base, *proof.z, proof.z_delta);
  if (left != right) {
    assert(false);
    return false;
  }

  auto const& tau = com_pub.tau;
  auto const& beta = com_ext_pub.beta;
  left = tau * challenge + beta;
  auto ip_za = InnerProduct(*proof.z, input.a);
  right = ComputeCommitment(pds_base, ip_za, proof.z_beta);
  if (left != right) {
    assert(false);
    return false;
  }
  return true;
}

// com(n) + com(1)
inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input, PdsBase const& pds_base) {
  Tick tick(__FUNCTION__);
  using details::ComputeCommitment;
  com_sec.r_xi = FrRand();
  com_sec.r_tau = FrRand();
  com_pub.xi = ComputeCommitment(pds_base, input.x, com_sec.r_xi);
  com_pub.tau = ComputeCommitment(pds_base, input.y, com_sec.r_tau);
  std::cout << Tick::GetIndentString() << "multiexp(" << input.n() << ")\n";
}

// com(n) + com(1) + ip(n)
inline void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                 CommitmentExtSec& com_ext_sec,
                                 ProverInput const& input,
                                 PdsBase const& pds_base) {
  Tick tick(__FUNCTION__);
  using details::ComputeCommitment;
  auto n = input.n();
  com_ext_sec.d.resize(n);
  FrRand(com_ext_sec.d.data(), n);
  com_ext_sec.r_beta = FrRand();
  com_ext_sec.r_delta = FrRand();
  com_ext_pub.delta =
      ComputeCommitment(pds_base, com_ext_sec.d, com_ext_sec.r_delta);
  com_ext_pub.beta = ComputeCommitment(
      pds_base, InnerProduct(input.a, com_ext_sec.d), com_ext_sec.r_beta);
  std::cout << Tick::GetIndentString() << "multiexp(" << input.n() << ")\n";
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub) {
  using details::HashUpdate;
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.xi);
  HashUpdate(hash, com_pub.tau);
  HashUpdate(hash, com_ext_pub.beta);
  HashUpdate(hash, com_ext_pub.delta);
  hash.Final(seed.data());
}

inline void ComputeProof(Proof& proof, ProverInput const& input,
                         CommitmentSec const& com_sec,
                         CommitmentExtSec const& com_ext_sec,
                         Fr const& challenge) {
  // z = c * x + d
  proof.z = std::make_shared<std::vector<Fr>>();
  details::VectorMul(*proof.z, input.x, challenge);
  details::VectorInc(*proof.z, com_ext_sec.d);
  proof.z_delta = challenge * com_sec.r_xi + com_ext_sec.r_delta;
  proof.z_beta = challenge * com_sec.r_tau + com_ext_sec.r_beta;
}

inline void RomProve(RomProof& rom_proof, PdsBase const& pds_base,
                     h256_t const& common_seed, ProverInput const& input,
                     CommitmentPub const& com_pub,
                     CommitmentSec const& com_sec) {
  Tick tick(__FUNCTION__);

  assert(pds_base.size() >= input.n());

  CommitmentExtSec com_ext_sec;
  ComputeCommitmentExt(rom_proof.com_ext_pub, com_ext_sec, input, pds_base);

  auto seed = common_seed;
  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  ComputeProof(rom_proof.proof, input, com_sec, com_ext_sec, challenge);
}

inline bool RomVerify(RomProof const& rom_proof, PdsBase const& pds_base,
                      h256_t const& common_seed, VerifierInput const& input) {
  Tick tick(__FUNCTION__);
  assert(pds_base.size() >= rom_proof.n());
  if (input.a.size() != rom_proof.proof.z->size() || input.a.empty())
    return false;

  auto seed = common_seed;
  UpdateSeed(seed, input.com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  return VerifyInternal(pds_base, input, challenge, rom_proof.com_ext_pub,
                        rom_proof.proof);
}

inline bool TestRom(int64_t n) {
  auto pds_base = GetEccPub().GetPdsBase(n);
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  h256_t UpdateSeed = misc::RandH256();

  ProverInput prover_input(x, a);

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input, pds_base);

  RomProof rom_proof;
  RomProve(rom_proof, pds_base, UpdateSeed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(a, com_pub);
  return RomVerify(rom_proof, pds_base, UpdateSeed, verifier_input);
}
}  // namespace hyrax::a2