#pragma once

#include "./details.h"

// x, y, z: secret Fr
// open: com(x), com(y), com(z)
// prove: z = x*y
// proof size: 3 G1 and 5 Fr
// prove cost: 6 eccmul
// verify cost: 6 eccmul
namespace hyrax::a1 {
struct ProverInput {
  ProverInput(Fr const& x, Fr const& y, Fr const& z) : x(x), y(y), z(z) {
    assert(z == x * y);
  }
  ProverInput(Fr const& x, Fr const& y) : x(x), y(y) { z = x * y; }
  Fr x;
  Fr y;
  Fr z;
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& x, G1 const& y, G1 const& z) : x(x), y(y), z(z) {}
  G1 x;  // com(x,r_x)
  G1 y;  // com(y,r_y)
  G1 z;  // com(z,r_z)
};

struct CommitmentSec {
  CommitmentSec() {}
  CommitmentSec(Fr const& r_x, Fr const& r_y, Fr const& r_z)
      : r_x(r_x), r_y(r_y), r_z(r_z) {}
  Fr r_x;
  Fr r_y;
  Fr r_z;
};

struct CommitmentExtPub {
  CommitmentExtPub() {}
  CommitmentExtPub(G1 const& alpha, G1 const& beta, G1 const& delta)
      : alpha(alpha), beta(beta), delta(delta) {}
  G1 alpha;  // com(b1, b2)
  G1 beta;   // com(b3, b4)
  G1 delta;  // com_pub.x * b3 + h * b5
};

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.alpha == right.alpha && left.beta == right.beta &&
         left.delta == right.delta;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

struct CommitmentExtSec {
  Fr b1;
  Fr b2;
  Fr b3;
  Fr b4;
  Fr b5;
};

struct Proof {
  Fr z1;
  Fr z2;
  Fr z3;
  Fr z4;
  Fr z5;
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.z1 == right.z1 && left.z2 == right.z2 && left.z3 == right.z3 &&
         left.z4 == right.z4 && left.z5 == right.z5;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

struct RomProof {
  CommitmentExtPub com_ext_pub;  // 3 G1
  Proof proof;                   // 5 Fr
};

inline bool operator==(RomProof const& left, RomProof const& right) {
  return left.com_ext_pub == right.com_ext_pub && left.proof == right.proof;
}

inline bool operator!=(RomProof const& left, RomProof const& right) {
  return !(left == right);
}

struct VerifierInput {
  VerifierInput(CommitmentPub const& com_pub) : com_pub(com_pub) {}
  CommitmentPub const& com_pub;
};

inline bool VerifyInternal(VerifierInput const& input, Fr const& c,
                           CommitmentExtPub const& com_ext_pub,
                           Proof const& proof) {
  // Tick tick(__FUNCTION__);
  using details::ComputeCommitment;

  auto const& g = details::GetPdsG();
  auto const& h = details::GetPdsH();
  auto const& com_pub = input.com_pub;

  std::vector<parallel::Task> tasks(3);
  bool ret0 = false;
  tasks[0] = [&ret0, &com_pub, &com_ext_pub, &c, &proof, &g, &h]() {
    G1 left = com_ext_pub.alpha + com_pub.x * c;
    G1 right = g * proof.z1 + h * proof.z2;
    ret0 = left == right;
  };

  bool ret1 = false;
  tasks[1] = [&ret1, &com_pub, &com_ext_pub, &c, &proof, &g, &h]() {
    G1 left = com_ext_pub.beta + com_pub.y * c;
    G1 right = g * proof.z3 + h * proof.z4;
    ret1 = left == right;
  };

  bool ret2 = false;
  tasks[2] = [&ret2, &com_pub, &com_ext_pub, &c, &proof, &g, &h]() {
    G1 left = com_ext_pub.delta + com_pub.z * c;
    G1 right = com_pub.x * proof.z3 + h * proof.z5;
    ret2 = left == right;
  };

  parallel::Invoke(tasks, true);

  assert(ret0 && ret1 && ret2);
  return ret0 && ret1 && ret2;
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  auto const& g = details::GetPdsG();
  auto const& h = details::GetPdsH();

  com_sec.r_x = FrRand();
  com_sec.r_y = FrRand();
  com_sec.r_z = FrRand();

  std::vector<parallel::Task> tasks(3);
  tasks[0] = [&com_pub, &input, &com_sec, &g, &h]() {
    com_pub.x = g * input.x + h * com_sec.r_x;
  };
  tasks[1] = [&com_pub, &input, &com_sec, &g, &h]() {
    com_pub.y = g * input.y + h * com_sec.r_y;
  };
  tasks[2] = [&com_pub, &input, &com_sec, &g, &h]() {
    com_pub.z = g * input.z + h * com_sec.r_z;
  };

  parallel::Invoke(tasks, true);
}

inline void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                 CommitmentExtSec& com_ext_sec,
                                 CommitmentPub const& com_pub) {
  // Tick tick(__FUNCTION__);
  auto const& g = details::GetPdsG();
  auto const& h = details::GetPdsH();

  com_ext_sec.b1 = FrRand();
  com_ext_sec.b2 = FrRand();
  com_ext_sec.b3 = FrRand();
  com_ext_sec.b4 = FrRand();
  com_ext_sec.b5 = FrRand();

  std::vector<parallel::Task> tasks(3);
  tasks[0] = [&com_ext_pub, &com_ext_sec, &g, &h]() {
    com_ext_pub.alpha = g * com_ext_sec.b1 + h * com_ext_sec.b2;
  };
  tasks[1] = [&com_ext_pub, &com_ext_sec, &g, &h]() {
    com_ext_pub.beta = g * com_ext_sec.b3 + h * com_ext_sec.b4;
  };
  tasks[2] = [&com_ext_pub, &com_pub, &com_ext_sec, &h]() {
    com_ext_pub.delta = com_pub.x * com_ext_sec.b3 + h * com_ext_sec.b5;
  };
  parallel::Invoke(tasks, true);
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub) {
  using details::HashUpdate;
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.x);
  HashUpdate(hash, com_pub.y);
  HashUpdate(hash, com_pub.z);
  HashUpdate(hash, com_ext_pub.alpha);
  HashUpdate(hash, com_ext_pub.beta);
  HashUpdate(hash, com_ext_pub.delta);
  hash.Final(seed.data());
}

inline void ComputeProof(Proof& proof, ProverInput const& input,
                         CommitmentSec const& com_sec,
                         CommitmentExtSec const& com_ext_sec,
                         Fr const& c) {
  proof.z1 = com_ext_sec.b1 + c * input.x;
  proof.z2 = com_ext_sec.b2 + c * com_sec.r_x;
  proof.z3 = com_ext_sec.b3 + c * input.y;
  proof.z4 = com_ext_sec.b4 + c * com_sec.r_y;
  proof.z5 = com_ext_sec.b5 + c * (com_sec.r_z - com_sec.r_x * input.y);
}

inline void RomProve(RomProof& rom_proof, h256_t const& common_seed,
                     ProverInput input, CommitmentPub com_pub,
                     CommitmentSec com_sec) {
  //Tick tick(__FUNCTION__);

  CommitmentExtSec com_ext_sec;
  ComputeCommitmentExt(rom_proof.com_ext_pub, com_ext_sec, com_pub);

  auto seed = common_seed;
  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr c = H256ToFr(seed);

  ComputeProof(rom_proof.proof, input, com_sec, com_ext_sec, c);
}

inline bool RomVerify(RomProof const& rom_proof, h256_t const& common_seed,
                      VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  auto seed = common_seed;
  UpdateSeed(seed, input.com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);
  return VerifyInternal(input, challenge, rom_proof.com_ext_pub,
                        rom_proof.proof);
}

inline bool TestRom() {
  h256_t UpdateSeed = misc::RandH256();

  ProverInput prover_input(FrRand(), FrRand());

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, UpdateSeed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(com_pub);
  return RomVerify(rom_proof, UpdateSeed, verifier_input);
}
}  // namespace hyrax::a1