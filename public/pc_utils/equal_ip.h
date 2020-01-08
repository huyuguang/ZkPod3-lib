#pragma once

#include "./details.h"
#include "./types.h"
#include "hyrax/hyrax.h"

namespace pc_utils::equal_ip {

// x: secret vector<Fr>, size = n
// y: secret vector<Fr>, size = m
// a: public vector<Fr>, size = n
// b: public vector<Fr>, size = m
// open: com(x), com(y)
// prove: <x,a> == <y, b>
// two hyrax::a2

struct Proof {
  hyrax::a2::RomProof p1;
  hyrax::a2::RomProof p2;
  G1 com_z;
};
inline bool operator==(Proof const& a, Proof const& b) {
  return a.p1 == b.p1 && a.p2 == b.p2 && a.com_z == b.com_z;
}
inline bool operator!=(Proof const& a, Proof const& b) { return !(a == b); }

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

inline void Prove(Proof& proof, h256_t const& seed, std::vector<Fr> const& x,
                  std::vector<Fr> const& a, G1 const& com_x, Fr const& com_x_r,
                  std::vector<Fr> const& y, std::vector<Fr> const& b,
                  G1 const& com_y, Fr const& com_y_r, Fr const& z) {
  assert(!a.empty() && a.size() == x.size());
  assert(!b.empty() && b.size() == y.size());

#ifdef _DEBUG
  assert(com_x == PcComputeCommitment(x, com_x_r));
  assert(com_y == PcComputeCommitment(y, com_y_r));
  assert(z == InnerProduct(x, a));
  assert(z == InnerProduct(y, b));
#endif

  Fr r_tau = FrRand();
  proof.com_z = PcComputeCommitment(z, r_tau);

  hyrax::a2::ProverInput input1(x, a, z);
  hyrax::a2::CommitmentSec com_sec1;
  com_sec1.r_xi = com_x_r;
  com_sec1.r_tau = r_tau;
  hyrax::a2::CommitmentPub com_pub1;
  com_pub1.xi = com_x;
  com_pub1.tau = proof.com_z;

  hyrax::a2::ProverInput input2(y, b, z);
  hyrax::a2::CommitmentSec com_sec2;
  com_sec2.r_xi = com_y_r;
  com_sec2.r_tau = r_tau;
  hyrax::a2::CommitmentPub com_pub2;
  com_pub2.xi = com_y;
  com_pub2.tau = proof.com_z;

  std::array<parallel::Task, 2> tasks;

  tasks[0] = [&proof, &seed, &input1, &com_pub1, &com_sec1]() {
    hyrax::a2::RomProve(proof.p1, seed, input1, com_pub1, com_sec1);
  };

  tasks[1] = [&proof, &seed, &input2, &com_pub2, &com_sec2]() {
    hyrax::a2::RomProve(proof.p2, seed, input2, com_pub2, com_sec2);
  };

  parallel::Invoke(tasks);
}

inline bool Verify(h256_t const& seed, std::vector<Fr> const& a,
                   G1 const& com_x, std::vector<Fr> const& b, G1 const& com_y,
                   Proof const& proof) {
  std::array<std::atomic<bool>, 2> rets;
  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&proof, &com_x, &rets, &seed, &a]() {
    hyrax::a2::CommitmentPub com_pub;
    com_pub.xi = com_x;
    com_pub.tau = proof.com_z;
    hyrax::a2::VerifierInput input(a, com_pub);
    rets[0] = hyrax::a2::RomVerify(proof.p1, seed, input);
  };
  tasks[1] = [&proof, &com_y, &rets, &seed, &b]() {
    hyrax::a2::CommitmentPub com_pub;
    com_pub.xi = com_y;
    com_pub.tau = proof.com_z;
    hyrax::a2::VerifierInput input(b, com_pub);
    rets[1] = hyrax::a2::RomVerify(proof.p2, seed, input);
  };

  parallel::Invoke(tasks);

  if (!rets[0] || !rets[1]) {
    assert(false);
    return false;
  }
  return true;
}

inline bool Test() {
  auto seed = misc::RandH256();
  std::vector<Fr> x(10);
  FrRand(x);
  std::vector<Fr> a(10);
  FrRand(a);
  std::vector<Fr> y(3);
  FrRand(y);
  std::vector<Fr> b(3);
  FrRand(b);
  auto z = InnerProduct(x, a);
  auto temp = InnerProduct(y, b);
  y.push_back(FrRand());
  b.push_back((z - temp) * FrInv(y.back()));
  assert(z == InnerProduct(y, b));

  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitment(x, com_x_r);

  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitment(y, com_y_r);

  Proof proof;
  Prove(proof, seed, x, a, com_x, com_x_r, y, b, com_y, com_y_r, z);

  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);

  Proof proof2;
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream> ia(is);
  ia.serialize(proof2);

  assert(proof == proof2);

  return Verify(seed, a, com_x, b, com_y, proof2);
}
}  // namespace pc_utils::equal_ip