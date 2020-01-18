#pragma once

#include "./details.h"

// a: public vector<Fr>, size = n
// x: secret vector<Fr>, size = n
// y: secret Fr, = <x,a>
// open: com(gx,x), com(gy,y)
// prove: y=<x,a>
// proof size: 2 G1 and n+2 Fr
// prove cost: mulexp(n)
// verify cost: mulexp(n)
namespace hyrax::a2 {

struct ProverInput {
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& a, Fr const& y,
              int64_t x_g_offset, int64_t y_g_offset)
      : x(x), a(a), y(y), x_g_offset(x_g_offset), y_g_offset(y_g_offset) {
    assert(y == InnerProduct(x, a));
  }
  int64_t n() const { return (int64_t)x.size(); }
  std::vector<Fr> const& x;  // x.size = n
  std::vector<Fr> const& a;  // a.size = n
  Fr const y;                // y = <x, a>
  int64_t const x_g_offset;
  int64_t const y_g_offset;
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& xi, G1 const& tau) : xi(xi), tau(tau) {}
  G1 xi;   // com(gx, x, r_xi)
  G1 tau;  // com(gy, y, r_tau)
};
inline bool operator==(CommitmentPub const& left, CommitmentPub const& right) {
  return left.xi == right.xi && left.tau == right.tau;
}

inline bool operator!=(CommitmentPub const& left, CommitmentPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentPub const& t) {
  ar& YAS_OBJECT_NVP("a2.cp", ("xi", t.xi), ("tau", t.tau));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentPub& t) {
  ar& YAS_OBJECT_NVP("a2.cp", ("xi", t.xi), ("tau", t.tau));
}

struct CommitmentSec {
  CommitmentSec() {}
  CommitmentSec(Fr const& x, Fr const& t) : r_xi(x), r_tau(t) {}
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

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.delta == right.delta && left.beta == right.beta;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

struct CommitmentExtSec {
  std::vector<Fr> d;  // size = n
  Fr r_beta;
  Fr r_delta;
};

struct Proof {
  std::vector<Fr> z;  // z.size = n
  Fr z_delta;
  Fr z_beta;
  int64_t n() const { return (int64_t)z.size(); }
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.z == right.z && left.z_delta == right.z_delta &&
         left.z_beta == right.z_beta;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("a2.pf", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("a2.pf", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

struct RomProof {
  CommitmentExtPub com_ext_pub;  // 2 G1
  Proof proof;                   // n+2 Fr
  int64_t n() const { return proof.n(); }
};

inline bool operator==(RomProof const& left, RomProof const& right) {
  return left.com_ext_pub == right.com_ext_pub && left.proof == right.proof;
}

inline bool operator!=(RomProof const& left, RomProof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, RomProof const& t) {
  ar& YAS_OBJECT_NVP("a2.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, RomProof& t) {
  ar& YAS_OBJECT_NVP("a2.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

struct VerifierInput {
  VerifierInput(std::vector<Fr> const& a, CommitmentPub const& com_pub,
                int64_t x_g_offset, int64_t y_g_offset)
      : a(a),
        com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset) {}
  std::vector<Fr> const& a;  // a.size = n
  CommitmentPub const& com_pub;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
};

// com(n) + com(1) + ip(n)
inline bool VerifyInternal(VerifierInput const& input, Fr const& challenge,
                           CommitmentExtPub const& com_ext_pub,
                           Proof const& proof) {
  // Tick tick(__FUNCTION__);
  auto const& com_pub = input.com_pub;

  std::array<parallel::Task, 2> tasks;
  bool ret1 = false;
  tasks[0] = [&ret1, &com_pub, &com_ext_pub, &challenge, &proof,&input]() {
    auto const& xi = com_pub.xi;
    auto const& delta = com_ext_pub.delta;
    G1 left = xi * challenge + delta;
    G1 right = PcComputeCommitmentG(input.x_g_offset, proof.z, proof.z_delta);
    ret1 = left == right;
  };

  bool ret2 = false;
  tasks[1] = [&ret2, &com_pub, &com_ext_pub, &challenge, &proof, &input]() {
    auto const& tau = com_pub.tau;
    auto const& beta = com_ext_pub.beta;
    G1 left = tau * challenge + beta;
    auto ip_za = InnerProduct(proof.z, input.a);
    G1 right = PcComputeCommitmentG(input.y_g_offset, ip_za, proof.z_beta);
    ret2 = left == right;
  };
  parallel::Invoke(tasks, true);

  assert(ret1 && ret2);
  return ret1 && ret2;
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec const& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&com_pub, &input, &com_sec]() {
    com_pub.xi = PcComputeCommitmentG(input.x_g_offset, input.x, com_sec.r_xi);
  };
  tasks[1] = [&com_pub, &input, &com_sec]() {
    com_pub.tau = PcComputeCommitmentG(input.y_g_offset, input.y, com_sec.r_tau);
  };
  parallel::Invoke(tasks, true);
}

// com(n) + com(1) + ip(n)
inline void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                 CommitmentExtSec& com_ext_sec,
                                 ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  auto n = input.n();
  com_ext_sec.d.resize(n);
  FrRand(com_ext_sec.d.data(), n);
  com_ext_sec.r_beta = FrRand();
  com_ext_sec.r_delta = FrRand();

  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&com_ext_pub, &com_ext_sec, &input]() {
    com_ext_pub.delta = PcComputeCommitmentG(input.x_g_offset, com_ext_sec.d,
                                             com_ext_sec.r_delta);
  };
  tasks[1] = [&com_ext_pub, &input, &com_ext_sec]() {
    com_ext_pub.beta = PcComputeCommitmentG(input.y_g_offset,
        InnerProduct(input.a, com_ext_sec.d), com_ext_sec.r_beta);
  };
  parallel::Invoke(tasks, true);

  // std::cout << Tick::GetIndentString() << "multiexp(" << input.n() << ")\n";
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub) {
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
  proof.z = input.x * challenge + com_ext_sec.d;
  proof.z_delta = challenge * com_sec.r_xi + com_ext_sec.r_delta;
  proof.z_beta = challenge * com_sec.r_tau + com_ext_sec.r_beta;
}

// TODO: should ProverInput const& input?
inline void RomProve(RomProof& rom_proof, h256_t seed, ProverInput input,
                     CommitmentPub com_pub, CommitmentSec com_sec) {
  // Tick tick(__FUNCTION__);

  assert(PcBase::kGSize >= input.n());

  CommitmentExtSec com_ext_sec;
  ComputeCommitmentExt(rom_proof.com_ext_pub, com_ext_sec, input);

  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  ComputeProof(rom_proof.proof, input, com_sec, com_ext_sec, challenge);
}

inline bool RomVerify(RomProof const& rom_proof, h256_t seed,
                      VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  assert(PcBase::kGSize >= rom_proof.n());
  if (input.a.size() != rom_proof.proof.z.size() || input.a.empty())
    return false;

  UpdateSeed(seed, input.com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  return VerifyInternal(input, challenge, rom_proof.com_ext_pub,
                        rom_proof.proof);
}

inline bool TestRom(int64_t n) {
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  h256_t UpdateSeed = misc::RandH256();

  int64_t x_g_offset = 30;
  int64_t y_g_offset = -1;
  auto z = InnerProduct(x, a);
  ProverInput prover_input(x, a, z, x_g_offset, y_g_offset);

  CommitmentPub com_pub;
  CommitmentSec com_sec(FrRand(), FrRand());
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, UpdateSeed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(a, com_pub, x_g_offset, y_g_offset);
  bool success = RomVerify(rom_proof, UpdateSeed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace hyrax::a2