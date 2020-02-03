#pragma once

#include "./details.h"

// x, y, z: secret Fr
// open: com(gx,x), com(gy,y), com(gx,z), gx can equal gy
// prove: z = x*y
// proof size: 3 G1 and 5 Fr
// prove cost: 6 eccmul
// verify cost: 6 eccmul
namespace hyrax::a1 {
struct ProverInput {
  ProverInput(Fr const& x, Fr const& y, Fr const& z, int64_t x_g_offset,
              int64_t y_g_offset)
      : x(x), y(y), z(z), x_g_offset(x_g_offset), y_g_offset(y_g_offset) {
    assert(z == x * y);
  }
  Fr const x;
  Fr const y;
  Fr const z;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
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

struct SubProof {
  Fr z1;
  Fr z2;
  Fr z3;
  Fr z4;
  Fr z5;
};

inline bool operator==(SubProof const& left, SubProof const& right) {
  return left.z1 == right.z1 && left.z2 == right.z2 && left.z3 == right.z3 &&
         left.z4 == right.z4 && left.z5 == right.z5;
}

inline bool operator!=(SubProof const& left, SubProof const& right) {
  return !(left == right);
}

struct Proof {
  CommitmentExtPub com_ext_pub;  // 3 G1
  SubProof sub_proof;            // 5 Fr
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.com_ext_pub == right.com_ext_pub &&
         left.sub_proof == right.sub_proof;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

struct VerifierInput {
  VerifierInput(CommitmentPub const& com_pub, int64_t x_g_offset,
                int64_t y_g_offset)
      : com_pub(com_pub), x_g_offset(x_g_offset), y_g_offset(y_g_offset) {}
  CommitmentPub const& com_pub;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
};

inline bool VerifyInternal(VerifierInput const& input, Fr const& c,
                           CommitmentExtPub const& com_ext_pub,
                           SubProof const& sub_proof) {
  // Tick tick(__FUNCTION__);
  auto const& gx = PcG(input.x_g_offset);
  auto const& gy = PcG(input.y_g_offset);
  auto const& h = PcH();
  auto const& com_pub = input.com_pub;

  std::array<parallel::Task, 3> tasks;
  bool ret0 = false;
  tasks[0] = [&ret0, &com_pub, &com_ext_pub, &c, &sub_proof, &gx, &h]() {
    G1 left = com_ext_pub.alpha + com_pub.x * c;
    G1 right = gx * sub_proof.z1 + h * sub_proof.z2;
    ret0 = left == right;
  };

  bool ret1 = false;
  tasks[1] = [&ret1, &com_pub, &com_ext_pub, &c, &sub_proof, &gy, &h]() {
    G1 left = com_ext_pub.beta + com_pub.y * c;
    G1 right = gy * sub_proof.z3 + h * sub_proof.z4;
    ret1 = left == right;
  };

  bool ret2 = false;
  tasks[2] = [&ret2, &com_pub, &com_ext_pub, &c, &sub_proof, &h]() {
    G1 left = com_ext_pub.delta + com_pub.z * c;
    G1 right = com_pub.x * sub_proof.z3 + h * sub_proof.z5;
    ret2 = left == right;
  };

  parallel::Invoke(tasks, true);

  assert(ret0 && ret1 && ret2);
  return ret0 && ret1 && ret2;
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  auto const& gx = PcG(input.x_g_offset);
  auto const& gy = PcG(input.y_g_offset);
  auto const& gz = gx;
  auto const& h = PcH();

  com_sec.r_x = FrRand();
  com_sec.r_y = FrRand();
  com_sec.r_z = FrRand();

  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_pub, &input, &com_sec, &gx, &h]() {
    com_pub.x = gx * input.x + h * com_sec.r_x;
  };
  tasks[1] = [&com_pub, &input, &com_sec, &gy, &h]() {
    com_pub.y = gy * input.y + h * com_sec.r_y;
  };
  tasks[2] = [&com_pub, &input, &com_sec, &gz, &h]() {
    com_pub.z = gz * input.z + h * com_sec.r_z;
  };

  parallel::Invoke(tasks, true);
}

inline void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                 CommitmentExtSec& com_ext_sec,
                                 ProverInput const& input,
                                 CommitmentPub const& com_pub) {
  // Tick tick(__FUNCTION__);
  auto const& gx = PcG(input.x_g_offset);
  auto const& gy = PcG(input.y_g_offset);
  auto const& h = PcH();

  com_ext_sec.b1 = FrRand();
  com_ext_sec.b2 = FrRand();
  com_ext_sec.b3 = FrRand();
  com_ext_sec.b4 = FrRand();
  com_ext_sec.b5 = FrRand();

  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_ext_pub, &com_ext_sec, &gx, &h]() {
    com_ext_pub.alpha = gx * com_ext_sec.b1 + h * com_ext_sec.b2;
  };
  tasks[1] = [&com_ext_pub, &com_ext_sec, &gy, &h]() {
    com_ext_pub.beta = gy * com_ext_sec.b3 + h * com_ext_sec.b4;
  };
  tasks[2] = [&com_ext_pub, &com_pub, &com_ext_sec, &h]() {
    com_ext_pub.delta = com_pub.x * com_ext_sec.b3 + h * com_ext_sec.b5;
  };
  parallel::Invoke(tasks, true);
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub) {
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

inline void ComputeSubProof(SubProof& sub_proof, ProverInput const& input,
                            CommitmentSec const& com_sec,
                            CommitmentExtSec const& com_ext_sec, Fr const& c) {
  sub_proof.z1 = com_ext_sec.b1 + c * input.x;
  sub_proof.z2 = com_ext_sec.b2 + c * com_sec.r_x;
  sub_proof.z3 = com_ext_sec.b3 + c * input.y;
  sub_proof.z4 = com_ext_sec.b4 + c * com_sec.r_y;
  sub_proof.z5 = com_ext_sec.b5 + c * (com_sec.r_z - com_sec.r_x * input.y);
}

inline void Prove(Proof& proof, h256_t seed, ProverInput input,
                  CommitmentPub com_pub, CommitmentSec com_sec) {
  // Tick tick(__FUNCTION__);

  CommitmentExtSec com_ext_sec;
  ComputeCommitmentExt(proof.com_ext_pub, com_ext_sec, input, com_pub);

  UpdateSeed(seed, com_pub, proof.com_ext_pub);
  Fr c = H256ToFr(seed);

  ComputeSubProof(proof.sub_proof, input, com_sec, com_ext_sec, c);
}

inline bool Verify(Proof const& proof, h256_t seed,
                   VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  UpdateSeed(seed, input.com_pub, proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);
  return VerifyInternal(input, challenge, proof.com_ext_pub, proof.sub_proof);
}

inline bool Test() {
  Tick tick(__FUNCTION__);
  h256_t UpdateSeed = misc::RandH256();

  int64_t x_g_offset = 1;
  int64_t y_g_offset = 10;
  Fr x = FrRand();
  Fr y = FrRand();
  Fr z = x * y;
  ProverInput prover_input(x, y, z, x_g_offset, y_g_offset);

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  Proof proof;
  Prove(proof, UpdateSeed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(com_pub, x_g_offset, y_g_offset);
  bool success = Verify(proof, UpdateSeed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace hyrax::a1