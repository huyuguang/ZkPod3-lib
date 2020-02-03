#pragma once

#include "groth09/details.h"
#include "hyrax/a3.h"
#include "pc_utils/equal_ip.h"
#include "utils/fst.h"
#include "bp/bp.h"

// t: public vector<Fr>, size = n
// X, Y: secret vector<Fr>, size = n
// z: secret Fr
// open: com(gx, X), com(gy, Y), com(gz, z). gx, gy maybe overlapped, gz can not
// prove: z = <X, Y o t>
// proof size: 2log(n) Fr and 4 G1
// prove cost: 2*mulexp(n)
// verify cost: mulexp(n)

namespace groth09::sec51c {

inline int64_t SelectYtOffset(int64_t n, int64_t x_g_offset) {
  auto kGSize = PcBase::kGSize;
  if (x_g_offset >= n) return x_g_offset - n;
  if (x_g_offset + 2 * n < kGSize) return x_g_offset + n;
  // since kMaxUnitPerZkp <= (int64_t)PcBase::kGSize/3,
  // the exception should never been throw
  std::cout << __FUNCTION__ << " oops, n: " << n
            << ", x_g_offset: " << x_g_offset << "\n";
  throw std::runtime_error("invalid x_g_offset");
}

struct ProverInput {
  std::vector<Fr> const& x;   // size = n
  std::vector<Fr> const& y;   // size = n
  std::vector<Fr> const& t;   // size = n
  std::vector<Fr> const& yt;  // yt = y o t
  Fr const z;                 // z = <x, y o t>
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;

  int64_t n() const { return (int64_t)x.size(); }
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& y,
              std::vector<Fr> const& t, std::vector<Fr> const& yt, Fr const& z,
              int64_t x_g_offset, int64_t y_g_offset, int64_t z_g_offset)
      : x(x),
        y(y),
        t(t),
        yt(yt),
        z(z),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {
    assert(!x.empty());
    assert(x.size() == y.size());
    assert(x.size() == t.size());
    assert(x.size() == yt.size());
    assert(yt == HadamardProduct(y, t));
    assert(z == InnerProduct(x, yt));
  }
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& a, G1 const& b, G1 const& c) : a(a), b(b), c(c) {}
  G1 a;
  G1 b;
  G1 c;
};

struct CommitmentSec {
  CommitmentSec() {}
  CommitmentSec(Fr const& r, Fr const& s, Fr const& t) : r(r), s(s), t(t) {}
  Fr r;
  Fr s;
  Fr t;
};

struct Proof {
  G1 com_yt;
  pc_utils::equal_ip::Proof proof_eip;
  bp::p3::Proof proof_bp;
  bool CheckFormat() const { return true; }
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.com_yt == right.com_yt && left.proof_eip == right.proof_eip &&
         left.proof_bp == right.proof_bp;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("51c.pf", ("c", t.com_yt), ("eip", t.proof_eip),
                     ("bp", t.proof_bp));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("51c.pf", ("c", t.com_yt), ("eip", t.proof_eip),
                     ("bp", t.proof_bp));
}

void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                std::vector<Fr> const& t) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, com_pub.a);
  HashUpdate(hash, com_pub.b);
  HashUpdate(hash, com_pub.c);
  HashUpdate(hash, t);
  hash.Final(seed.data());
}

inline void ProveYt(Proof& proof, h256_t seed, ProverInput const& input,
                    CommitmentPub const& com_pub, CommitmentSec const& com_sec,
                    Fr const& com_yt_r, int64_t yt_g_offset) {
  std::vector<Fr> e(input.n());
  ComputeFst1(seed, "sec51c", e);

  Fr yte = InnerProduct(input.yt, e);
  auto et = HadamardProduct(e, input.t);
  assert(yte == InnerProduct(input.y, et));

  // prove <yt, e> == <y, et>
  namespace eip = pc_utils::equal_ip;
  eip::ProverInput input_eip(input.yt, e, proof.com_yt, com_yt_r, yt_g_offset,
                             input.y, et, com_pub.b, com_sec.s,
                             input.y_g_offset, yte);
  eip::Prove(proof.proof_eip, seed, input_eip);
}

inline void Prove(Proof& proof, h256_t seed, ProverInput const& input,
                  CommitmentPub const& com_pub, CommitmentSec const& com_sec) {
  UpdateSeed(seed, com_pub, input.t);

  // prove yt = <y o t>
  Fr com_yt_r = FrRand();
  int64_t yt_g_offset = SelectYtOffset(input.n(), input.x_g_offset);
  proof.com_yt = PcComputeCommitmentG(yt_g_offset, input.yt, com_yt_r);
  ProveYt(proof, seed, input, com_pub, com_sec, com_yt_r, yt_g_offset);

  // now we have com(yt_g_offset, yt, com_yt_r), want to prove <x, yt> == z
  std::vector<G1> g1 = GetPcBase().CopyG(input.x_g_offset, input.n());
  std::vector<G1> g2 = GetPcBase().CopyG(yt_g_offset, input.n());
  std::vector<Fr> x_copy = input.x;
  std::vector<Fr> yt_copy = input.yt;
  bp::p3::ProverInput input_bp(std::move(g1), std::move(g2), PcH(),
                               PcG(input.z_g_offset), std::move(x_copy),
                               std::move(yt_copy), input.z);
  bp::p3::CommitmentPub com_pub_bp;
  bp::p3::CommitmentSec com_sec_bp;
  com_sec_bp.alpha = com_sec.r + com_yt_r;
  com_sec_bp.beta = com_sec.t;
  com_pub_bp.p = com_pub.a + proof.com_yt;
  com_pub_bp.q = com_pub.c;
  bp::p3::Prove(proof.proof_bp, seed, std::move(input_bp), com_pub_bp,
                com_sec_bp);
}

struct VerifierInput {
  VerifierInput(std::vector<Fr> const& t, CommitmentPub const& com_pub,
                int64_t x_g_offset, int64_t y_g_offset, int64_t z_g_offset)
      : t(t),
        com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {}
  int64_t n() const { return (int64_t)t.size(); }
  std::vector<Fr> const& t;
  CommitmentPub const& com_pub;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
};

inline bool VerifyYt(Proof const& proof, h256_t seed,
                     VerifierInput const& input, int64_t yt_g_offset) {
  std::vector<Fr> e(input.n());
  ComputeFst1(seed, "sec51c", e);
  auto et = HadamardProduct(e, input.t);

  namespace eip = pc_utils::equal_ip;
  eip::VerifierInput input_eip(e, proof.com_yt, yt_g_offset, et,
                               input.com_pub.b, input.y_g_offset);
  return eip::Verify(seed, proof.proof_eip, input_eip);
}

inline bool Verify(Proof const& proof, h256_t seed,
                   VerifierInput const& input) {
  UpdateSeed(seed, input.com_pub, input.t);

  int64_t yt_g_offset = SelectYtOffset(input.n(), input.x_g_offset);
  if (!VerifyYt(proof, seed, input, yt_g_offset)) {
    assert(false);
    return false;
  }

  std::vector<G1> g1 = GetPcBase().CopyG(input.x_g_offset, input.n());
  std::vector<G1> g2 = GetPcBase().CopyG(yt_g_offset, input.n());
  bp::p3::CommitmentPub com_pub_bp;
  com_pub_bp.p = input.com_pub.a + proof.com_yt;
  com_pub_bp.q = input.com_pub.c;
  bp::p3::VerifierInput input_p3(std::move(g1), std::move(g2), PcH(),
                                 PcG(input.z_g_offset), com_pub_bp);

  return bp::p3::Verify(proof.proof_bp, seed, std::move(input_p3));
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  com_sec.r = FrRand();
  com_sec.s = FrRand();
  com_sec.t = FrRand();

  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_pub, &input, &com_sec]() {
    com_pub.a = PcComputeCommitmentG(input.x_g_offset, input.x, com_sec.r);
  };
  tasks[1] = [&com_pub, &input, &com_sec]() {
    com_pub.b = PcComputeCommitmentG(input.y_g_offset, input.y, com_sec.s);
  };
  tasks[2] = [&com_pub, &input, &com_sec]() {
    com_pub.c = PcComputeCommitmentG(input.z_g_offset, input.z, com_sec.t);
  };
  parallel::Invoke(tasks);
}

inline bool Test(int64_t n) {
  Tick tick(__FUNCTION__);
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> y(n);
  FrRand(y.data(), n);
  std::vector<Fr> t(n);
  FrRand(t.data(), t.size());

  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 10;
  int64_t y_g_offset = 30;
  int64_t z_g_offset = -1;
  std::vector<Fr> yt(n);
  yt = HadamardProduct(y, t);
  Fr z = InnerProduct(x, yt);
  ProverInput prover_input(x, y, t, yt, z, x_g_offset, y_g_offset, z_g_offset);

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  Proof proof;
  Prove(proof, seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(t, com_pub, x_g_offset, y_g_offset, z_g_offset);
  bool success = Verify(proof, seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace groth09::sec51c
