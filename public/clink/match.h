#pragma once

#include "./details.h"
#include "./parallel_r1cs.h"
#include "circuit/match_gadget.h"

// x: vector<Fr>, size = n
// k: Fr
// y: vector<Fr>, size = n, {0,1}
// open com(gx, x), k, com(gx, y)
// prove y[i] = (k == x[i])? 1:0

namespace clink {

template <typename Policy>
struct Match {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using R1cs = typename clink::ParallelR1cs<Policy>;

  struct Proof {
    typename R1cs::Proof r1cs_proof;
    std::vector<G1> com_w;
    G1 const& com_x() const { return com_w.front(); }
    G1 const& com_y() const { return com_w.back(); }
    bool operator==(Proof const& b) const {
      return r1cs_proof == b.r1cs_proof && com_w == b.com_w;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("m.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("m.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
  };

  struct ProveInput {
    ProveInput(Fr const& k, std::vector<Fr> const& x, G1 const& com_x,
               Fr const& com_x_r, std::vector<Fr> const& y, G1 const& com_y,
               Fr const& com_y_r, int64_t g_offset)
        : k(k),
          x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          y(y),
          com_y(com_y),
          com_y_r(com_y_r),
          g_offset(g_offset),
          n((int64_t)x.size()),
          gadget(circuit::MatchGadget(pb, k)),
          s((int64_t)pb.num_variables()) {
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);

      w.resize(s);
      for (auto& i : w) i.resize(n);

      for (int64_t j = 0; j < n; ++j) {
        gadget.Assign(x[j]);
        assert(pb.is_satisfied());
        auto v = pb.full_variable_assignment();
        for (int64_t i = 0; i < s; ++i) {
          w[i][j] = v[i];
        }
        assert(w.front()[j] == x[j]);
        assert(w.back()[j] == y[j]);
      }
    }
    Fr const& k;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    std::vector<Fr> const& y;
    G1 const& com_y;
    Fr const& com_y_r;
    int64_t const g_offset;

    int64_t const n;
    libsnark::protoboard<Fr> pb;
    circuit::MatchGadget gadget;
    int64_t const s;
    std::vector<std::vector<Fr>> mutable w;
  };

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input) {
    Tick tick(__FUNCTION__);

    std::vector<G1> com_w(input.s);
    std::vector<Fr> com_w_r(input.s);

    com_w.front() = input.com_x;
    com_w_r.front() = input.com_x_r;
    com_w.back() = input.com_y;
    com_w_r.back() = input.com_y_r;

    std::cout << "compute com(witness)\n";
    auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
      com_w_r[i] = FrRand();
      com_w[i] = PcComputeCommitmentG(input.g_offset, input.w[i], com_w_r[i]);
    };
    parallel::For<int64_t>(1, input.s - 1, parallel_f);

    typename R1cs::ProveInput r1cs_input(input.pb, std::move(input.w), com_w,
                                         com_w_r, input.g_offset);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    proof.com_w = std::move(com_w);
  }

  struct VerifyInput {
    VerifyInput(int64_t n, Fr const& k, int64_t g_offset)
        : n(n),
          k(k),
          g_offset(g_offset),
          gadget(circuit::MatchGadget(pb, k)),
          m((int64_t)pb.num_constraints()),
          s((int64_t)pb.num_variables()) {
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);
    }
    int64_t n;
    Fr const& k;
    int64_t const g_offset;
    libsnark::protoboard<Fr> pb;
    circuit::MatchGadget gadget;
    int64_t const m;
    int64_t const s;
    // since primary_input_size = 0, public_w is empty
    std::vector<std::vector<Fr>> public_w;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FUNCTION__);

    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }
    // if (proof.r1cs_proof.m() != input.m) {
    //  assert(false);
    //  return false;
    //}

    typename R1cs::VerifyInput pr_input(input.n, input.pb, proof.com_w,
                                        input.public_w, input.g_offset);
    return R1cs::Verify(proof.r1cs_proof, seed, pr_input);
  }

  static bool Test();
};

template <typename Policy>
bool Match<Policy>::Test() {
  auto seed = misc::RandH256();
  int64_t n = 100;
  std::vector<Fr> x(n);
  FrRand(x);
  Fr k = x[(int64_t)rand() % n];
  std::vector<Fr> y(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    y[i] = x[i] == k ? 1 : 0;
  }
  std::cout << "compute com(x), com(y)\n";
  int64_t g_offset = 30;
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(g_offset, x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitmentG(g_offset, y, com_y_r);

  Tick tick(__FUNCTION__);

  ProveInput prove_input(k, x, com_x, com_x_r, y, com_y, com_y_r, g_offset);
  Proof proof;
  Prove(proof, seed, prove_input);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  if (proof != proof2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  VerifyInput verify_input(n, k, g_offset);
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink