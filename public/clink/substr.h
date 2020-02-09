#pragma once

#include "./details.h"
#include "./parallel_r1cs.h"
#include "circuit/has0_gadget.h"
#include "circuit/substr_gadget.h"

// x: vector<Fr>, size = n
// k: Fr
// y: vector<Fr>, size = n, {0,1}
// open com(gx,x), k, com(gx,y)
// prove y[i] = IsSubStr(x[i],k)? 1:0

namespace clink {

template <typename Policy>
struct Substr {
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
      ar& YAS_OBJECT_NVP("s.p", ("hp", r1cs_proof), ("w", com_w));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("s.p", ("hp", r1cs_proof), ("w", com_w));
    }
  };

  struct ProveInput {
    ProveInput(std::string const& k, std::vector<Fr> const& x, G1 const& com_x,
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
          gadget(circuit::SubstrGadget(pb, k)),
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
    std::string const& k;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    std::vector<Fr> const& y;
    G1 const& com_y;
    Fr const& com_y_r;
    int64_t const g_offset;

    int64_t const n;
    libsnark::protoboard<Fr> pb;
    circuit::SubstrGadget gadget;
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

    {
      // Tick tick2("compute com(witness)");
      auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
        com_w_r[i] = FrRand();
        com_w[i] =
            PcComputeCommitmentG(input.g_offset, input.w[i], com_w_r[i], true);
      };
      parallel::For<int64_t>(1LL, input.s - 1, parallel_f);
    }

    typename R1cs::ProveInput r1cs_input(input.pb, std::move(input.w), com_w,
                                         com_w_r, input.g_offset);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    proof.com_w = std::move(com_w);
  }

  struct VerifyInput {
    VerifyInput(int64_t n, std::string const& k, int64_t g_offset)
        : n(n),
          k(k),
          g_offset(g_offset),
          gadget(circuit::SubstrGadget(pb, k)),
          m((int64_t)pb.num_constraints()),
          s((int64_t)pb.num_variables()) {
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);
    }
    int64_t n;
    std::string const& k;
    int64_t const g_offset;
    libsnark::protoboard<Fr> pb;
    circuit::SubstrGadget gadget;
    int64_t const m;
    int64_t const s;
    // since primary_input_size = 0, public_w is empty
    std::vector<std::vector<Fr>> public_w;
  };

  // NOTE: com_x and com_y can get by proof.com_x() and proof.com_y()
  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FUNCTION__);
    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }

    // if (proof.r1cs_proof.m() != (int64_t)misc::Pow2UB(input.m)) {
    //  assert(false);
    //  return false;
    //}

    typename ParallelR1cs<Policy>::VerifyInput pr_input(
        input.n, input.pb, proof.com_w, input.public_w, input.g_offset);
    return ParallelR1cs<Policy>::Verify(proof.r1cs_proof, seed, pr_input);
  }

  static bool Test();
};

template <typename Policy>
bool Substr<Policy>::Test() {
  auto seed = misc::RandH256();
  int64_t g_offset = 30;

  std::vector<std::string> x_str{{"1234567890123456789012345678901", "",
                                  "abcdefg", "34356356", "234qasdfaq44uUU",
                                  "AQW34352", "ASGDFUR7erweq", "asdfq45343",
                                  "q34asdfasfa9", "4534asq34afa"}};

#ifndef _DEBUG
  while (x_str.size() < 100000) {
    uint8_t buf[33];
    int len = rand() % 32;
    for (int i = 0; i < len; ++i) {
      buf[i] = (uint8_t)rand();
    }
    buf[len] = 0;
    x_str.push_back((char*)buf);
  }
#endif

  int64_t n = (int64_t)x_str.size();
  std::vector<Fr> x(n);
  for (int64_t i = 0; i < n; ++i) {
    x[i] = PackStrToFr(x_str[i].c_str());
  }
  std::string k = "343";
  std::vector<Fr> y(x.size());
  int64_t find_count = 0;
  for (size_t i = 0; i < x_str.size(); ++i) {
    bool find = x_str[i].find(k) != std::string::npos;
    if (find) ++find_count;
    y[i] = find ? 1 : 0;
  }
  std::cout << "find_count: " << find_count << "\n";

  Tick tick(__FUNCTION__);
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(g_offset, x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitmentG(g_offset, y, com_y_r, true);  // y is {0,1}

  ProveInput prove_input(k, x, com_x, com_x_r, y, com_y, com_y_r, g_offset);
  Proof proof;
  Prove(proof, seed, prove_input);

  if (proof.com_x() != com_x || proof.com_y() != com_y) {
    assert(false);
    return false;
  }

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