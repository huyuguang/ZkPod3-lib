#pragma once

#include "./details.h"
#include "./types.h"
#include "hyrax/hyrax.h"

namespace pc_utils {

// x: secret vector<Fr>, size = n
// y: secret vector<Fr>, size = m
// a: public vector<Fr>, size = n
// b: public vector<Fr>, size = m
// open: com(gx,x), com(gy,y)
// prove: <x,a> == <y, b>
// two HyraxA2

// HyraxA: hyrax::A2 or hyrax::A3
template <typename HyraxA>
struct EqualIp {
  struct Proof {
    typename HyraxA::Proof p1;
    typename HyraxA::Proof p2;
    G1 com_z;
    bool operator==(Proof const& b) const {
      return p1 == b.p1 && p2 == b.p2 && com_z == b.com_z;
    }
    bool operator!=(Proof const& b) const { return !(*this == b); }
  };

  struct ProverInput {
    std::vector<Fr> const& x;
    std::vector<Fr> const& a;
    G1 const& com_x;
    Fr const& com_x_r;
    int64_t const x_g_offset;
    std::vector<Fr> const& y;
    std::vector<Fr> const& b;
    G1 const& com_y;
    Fr const& com_y_r;
    int64_t const y_g_offset;
    Fr const& z;
    int64_t const z_g_offset = -1;

    ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& a,
                G1 const& com_x, Fr const& com_x_r, int64_t x_g_offset,
                std::vector<Fr> const& y, std::vector<Fr> const& b,
                G1 const& com_y, Fr const& com_y_r, int64_t y_g_offset,
                Fr const& z)
        : x(x),
          a(a),
          com_x(com_x),
          com_x_r(com_x_r),
          x_g_offset(x_g_offset),
          y(y),
          b(b),
          com_y(com_y),
          com_y_r(com_y_r),
          y_g_offset(y_g_offset),
          z(z) {
#ifdef _DEBUG
      assert(!a.empty() && a.size() == x.size());
      assert(!b.empty() && b.size() == y.size());
      assert(com_x == PcComputeCommitmentG(x_g_offset, x, com_x_r));
      assert(com_y == PcComputeCommitmentG(y_g_offset, y, com_y_r));
      assert(z == InnerProduct(x, a));
      assert(z == InnerProduct(y, b));
#endif
    }
  };

  static void Prove(Proof& proof, h256_t const& seed,
                    ProverInput const& input) {
    Fr r_tau = FrRand();
    proof.com_z = PcComputeCommitmentG(input.z_g_offset, input.z, r_tau);

    HyraxA::ProverInput input1(input.x, input.a, input.z, input.x_g_offset,
                               input.z_g_offset);
    HyraxA::CommitmentSec com_sec1;
    com_sec1.r_xi = input.com_x_r;
    com_sec1.r_tau = r_tau;
    HyraxA::CommitmentPub com_pub1;
    com_pub1.xi = input.com_x;
    com_pub1.tau = proof.com_z;

    HyraxA::ProverInput input2(input.y, input.b, input.z, input.y_g_offset,
                               input.z_g_offset);
    HyraxA::CommitmentSec com_sec2;
    com_sec2.r_xi = input.com_y_r;
    com_sec2.r_tau = r_tau;
    HyraxA::CommitmentPub com_pub2;
    com_pub2.xi = input.com_y;
    com_pub2.tau = proof.com_z;

    std::array<parallel::Task, 2> tasks;

    tasks[0] = [&proof, &seed, &input1, &com_pub1, &com_sec1]() {
      HyraxA::Prove(proof.p1, seed, input1, com_pub1, com_sec1);
    };

    tasks[1] = [&proof, &seed, &input2, &com_pub2, &com_sec2]() {
      HyraxA::Prove(proof.p2, seed, input2, com_pub2, com_sec2);
    };

    parallel::Invoke(tasks);
  }

  struct VerifierInput {
    VerifierInput(std::vector<Fr> const& a, G1 const& com_x, int64_t x_g_offset,
                  std::vector<Fr> const& b, G1 const& com_y, int64_t y_g_offset)
        : a(a),
          com_x(com_x),
          x_g_offset(x_g_offset),
          b(b),
          com_y(com_y),
          y_g_offset(y_g_offset) {}
    std::vector<Fr> const& a;
    G1 const& com_x;
    int64_t const x_g_offset;
    std::vector<Fr> const& b;
    G1 const& com_y;
    int64_t const y_g_offset;
    int64_t const z_g_offset = -1;
  };

  static bool Verify(h256_t const& seed, Proof const& proof,
                     VerifierInput const& input) {
    std::array<std::atomic<bool>, 2> rets;
    std::array<parallel::Task, 2> tasks;
    tasks[0] = [&proof, &input, &rets, &seed]() {
      HyraxA::CommitmentPub com_pub;
      com_pub.xi = input.com_x;
      com_pub.tau = proof.com_z;
      HyraxA::VerifierInput a_input(input.a, com_pub, input.x_g_offset,
                                    input.z_g_offset);
      rets[0] = HyraxA::Verify(proof.p1, seed, a_input);
    };
    tasks[1] = [&proof, &input, &rets, &seed]() {
      HyraxA::CommitmentPub com_pub;
      com_pub.xi = input.com_y;
      com_pub.tau = proof.com_z;
      HyraxA::VerifierInput a_input(input.b, com_pub, input.y_g_offset,
                                    input.z_g_offset);
      rets[1] = HyraxA::Verify(proof.p2, seed, a_input);
    };

    parallel::Invoke(tasks);

    if (!rets[0] || !rets[1]) {
      assert(false);
      return false;
    }
    return true;
  }

  static bool Test(int64_t xn, int64_t yn);
};

//// save to bin
//template <typename Ar, typename HyraxA>
//void serialize(Ar& ar, typename EqualIp<HyraxA>::Proof const& t) {
//  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
//}
//
//// load from bin
//template <typename Ar, typename HyraxA>
//void serialize(Ar& ar, typename EqualIp<HyraxA>::Proof& t) {
//  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
//}

// save to bin
template <typename Ar>
void serialize(Ar& ar, typename EqualIp<hyrax::A3>::Proof const& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, typename EqualIp<hyrax::A3>::Proof& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, typename EqualIp<hyrax::A2>::Proof const& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, typename EqualIp<hyrax::A2>::Proof& t) {
  ar& YAS_OBJECT_NVP("ei.p", ("p1", t.p1), ("p2", t.p2), ("z", t.com_z));
}

template <typename HyraxA>
bool EqualIp<HyraxA>::Test(int64_t xn, int64_t yn) {
  auto seed = misc::RandH256();
  std::vector<Fr> x(xn);
  FrRand(x);
  std::vector<Fr> a(xn);
  FrRand(a);
  std::vector<Fr> y(yn);
  FrRand(y);
  std::vector<Fr> b(yn);
  FrRand(b);
  auto z = InnerProduct(x, a);
  auto temp = InnerProduct(y, b);
  y.push_back(FrRand());
  b.push_back((z - temp) * FrInv(y.back()));
  assert(z == InnerProduct(y, b));

  int64_t x_g_offset = 10;
  int64_t y_g_offset = 40;

  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(x_g_offset, x, com_x_r);

  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitmentG(y_g_offset, y, com_y_r);

  ProverInput prover_input(x, a, com_x, com_x_r, x_g_offset, y, b, com_y,
                           com_y_r, y_g_offset, z);
  Proof proof;
  Prove(proof, seed, prover_input);

  // test serialize
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);

  Proof proof2;
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  ia.serialize(proof2);

  assert(proof == proof2);

  VerifierInput verifier_input(a, com_x, x_g_offset, b, com_y, y_g_offset);
  bool success = Verify(seed, proof, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace pc_utils
