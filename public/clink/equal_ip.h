#pragma once

#include "./details.h"
#include "hyrax/hyrax.h"

namespace clink {

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

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("ei.p", ("p1", p1), ("p2", p2), ("z", com_z));
    }

    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("ei.p", ("p1", p1), ("p2", p2), ("z", com_z));
    }
  };

  struct ProveInput {
    std::vector<Fr> const& x;
    std::vector<Fr> const& a;
    G1 const& com_x;
    Fr const& com_x_r;
    GetRefG1 const& get_gx;
    std::vector<Fr> const& y;
    std::vector<Fr> const& b;
    G1 const& com_y;
    Fr const& com_y_r;
    GetRefG1 const& get_gy;
    Fr const& z;
    G1 const& gz = pc::PcU();

    ProveInput(std::vector<Fr> const& x, std::vector<Fr> const& a,
               G1 const& com_x, Fr const& com_x_r, GetRefG1 const& get_gx,
               std::vector<Fr> const& y, std::vector<Fr> const& b,
               G1 const& com_y, Fr const& com_y_r, GetRefG1 const& get_gy,
               Fr const& z)
        : x(x),
          a(a),
          com_x(com_x),
          com_x_r(com_x_r),
          get_gx(get_gx),
          y(y),
          b(b),
          com_y(com_y),
          com_y_r(com_y_r),
          get_gy(get_gy),
          z(z) {
#ifdef _DEBUG
      assert(!a.empty() && a.size() == x.size());
      assert(!b.empty() && b.size() == y.size());
      assert(com_x == pc::ComputeCom(get_gx, x, com_x_r));
      assert(com_y == pc::ComputeCom(get_gy, y, com_y_r));
      assert(z == InnerProduct(x, a));
      assert(z == InnerProduct(y, b));
#endif
    }
  };

  static void Prove(Proof& proof, h256_t const& seed, ProveInput const& input) {
    Fr r_tau = FrRand();
    proof.com_z = pc::ComputeCom(input.gz, input.z, r_tau);

    typename HyraxA::ProveInput input1("eip", input.x, input.a, input.z,
                                       input.get_gx, input.gz);
    typename HyraxA::CommitmentSec com_sec1;
    com_sec1.r_xi = input.com_x_r;
    com_sec1.r_tau = r_tau;
    typename HyraxA::CommitmentPub com_pub1;
    com_pub1.xi = input.com_x;
    com_pub1.tau = proof.com_z;

    typename HyraxA::ProveInput input2("eip", input.y, input.b, input.z,
                                       input.get_gy, input.gz);
    typename HyraxA::CommitmentSec com_sec2;
    com_sec2.r_xi = input.com_y_r;
    com_sec2.r_tau = r_tau;
    typename HyraxA::CommitmentPub com_pub2;
    com_pub2.xi = input.com_y;
    com_pub2.tau = proof.com_z;

    std::array<parallel::VoidTask, 2> tasks;

    tasks[0] = [&proof, &seed, &input1, &com_pub1, &com_sec1]() {
      HyraxA::Prove(proof.p1, seed, input1, com_pub1, com_sec1);
    };

    tasks[1] = [&proof, &seed, &input2, &com_pub2, &com_sec2]() {
      HyraxA::Prove(proof.p2, seed, input2, com_pub2, com_sec2);
    };

    parallel::Invoke(tasks);
  }

  struct VerifyInput {
    VerifyInput(std::vector<Fr> const& a, G1 const& com_x,
                GetRefG1 const& get_gx, std::vector<Fr> const& b,
                G1 const& com_y, GetRefG1 const& get_gy)
        : a(a),
          com_x(com_x),
          get_gx(get_gx),
          b(b),
          com_y(com_y),
          get_gy(get_gy) {}
    std::vector<Fr> const& a;
    G1 const& com_x;
    GetRefG1 const& get_gx;
    std::vector<Fr> const& b;
    G1 const& com_y;
    GetRefG1 const& get_gy;
    G1 const& gz = pc::PcU();
  };

  static bool Verify(h256_t const& seed, Proof const& proof,
                     VerifyInput const& input) {
    std::array<std::atomic<bool>, 2> rets;
    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&proof, &input, &rets, &seed]() {
      typename HyraxA::CommitmentPub com_pub;
      com_pub.xi = input.com_x;
      com_pub.tau = proof.com_z;
      typename HyraxA::VerifyInput a_input("eip", input.a, com_pub,
                                           input.get_gx, input.gz);
      rets[0] = HyraxA::Verify(proof.p1, seed, a_input);
    };
    tasks[1] = [&proof, &input, &rets, &seed]() {
      typename HyraxA::CommitmentPub com_pub;
      com_pub.xi = input.com_y;
      com_pub.tau = proof.com_z;
      typename HyraxA::VerifyInput a_input("eip", input.b, com_pub,
                                           input.get_gy, input.gz);
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

template <typename HyraxA>
bool EqualIp<HyraxA>::Test(int64_t xn, int64_t yn) {
  Tick tick(__FN__);
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
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  GetRefG1 get_gy = [y_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[y_g_offset + i];
  };

  Fr com_x_r = FrRand();
  G1 com_x = pc::ComputeCom(get_gx, x, com_x_r);

  Fr com_y_r = FrRand();
  G1 com_y = pc::ComputeCom(get_gy, y, com_y_r);

  ProveInput prove_input(x, a, com_x, com_x_r, get_gx, y, b, com_y, com_y_r,
                         get_gy, z);
  Proof proof;
  Prove(proof, seed, prove_input);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << Tick::GetIndentString()
            << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  CHECK(proof == proof2, "");
#endif

  VerifyInput verify_input(a, com_x, get_gx, b, com_y, get_gy);
  bool success = Verify(seed, proof, verify_input);
  std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink
