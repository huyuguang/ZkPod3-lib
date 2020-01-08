#pragma once

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./details.h"
#include "./equal_ip.h"
#include "./parallel_r1cs.h"
#include "./types.h"

// x: vector<Fr>, size = n
// k: Fr
// y: vector<Fr>, size = n, {0,1}
// open com(x), k, com(y)
// prove y[i] = k==x[i]? 1:0

namespace pc_utils::match {

// x2 = (k==x)? 0:1/(x-k);
// assert(y*(x-k))==0)
// assert(x2 * (x-k) == y - 1
// 3 vars, 2 constraints
class MatchGadget : public libsnark::gadget<Fr> {
 public:
  MatchGadget(libsnark::protoboard<Fr>& pb, Fr const& k)
      : libsnark::gadget<Fr>(pb, "match"), k_(k) {
    x_.allocate(pb, "x");
    x2_.allocate(pb, "x2");
    y_.allocate(pb, "y");
    generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> const& y() const { return y_; }

  void Assign(Fr const& x) {
    this->pb.val(x_) = x;
    generate_r1cs_witness();
    assert(this->pb.is_satisfied());
  }

 private:
  void generate_r1cs_constraints() {
    // y * (x-k)) == 0
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(y_, x_ - k_, FrZero()),
        "y * (x-k)) == 0");

    // x2 * (x-k) == 1-ret
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(x2_, x_ - k_, -y_ + FrOne()),
        "x2 * (x-k) == 1-ret");
  }

  void generate_r1cs_witness() {
    auto const& x = this->pb.val(x_);
    bool equal = x == k_;
    this->pb.val(y_) = equal ? 1 : 0;
    this->pb.val(x2_) = equal ? 0 : FrInv(x - k_);
  }

 private:
  Fr const& k_;
  libsnark::pb_variable<Fr> x_;
  libsnark::pb_variable<Fr> x2_; // k==x? 0:1/(x-k)
  libsnark::pb_variable<Fr> y_;
};

struct Proof {
  groth09::sec43::RomProof hp_proof;
  std::vector<G1> com_w;
  G1 const& com_x() const { return com_w.front(); }
  G1 const& com_y() const { return com_w.back(); }
};

inline bool operator==(Proof const& a, Proof const& b) {
  return a.hp_proof == b.hp_proof && a.com_w == b.com_w;
}

inline bool operator!=(Proof const& a, Proof const& b) { return !(a == b); }

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("m.p", ("hp", t.hp_proof), ("w", t.com_w));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("m.p", ("hp", t.hp_proof), ("w", t.com_w));
}

inline void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                  Fr const& k, G1 const& com_x, Fr const& com_x_r,
                  G1 const& com_y, Fr const& com_y_r) {
  Tick tick(__FUNCTION__);
  int64_t n = (int64_t)x.size();
  int64_t const primary_input_size = 0;
  libsnark::protoboard<Fr> pb;
  MatchGadget gadget(pb, k);
  pb.set_input_sizes(primary_input_size);
  int64_t s = (int64_t)pb.num_variables();

  std::vector<std::vector<Fr>> w(s);
  for (auto& i : w) i.resize(n);

  for (int64_t j = 0; j < n; ++j) {
    gadget.Assign(x[j]);
    assert(pb.is_satisfied());
    auto v = pb.full_variable_assignment();
    for (int64_t i = 0; i < s; ++i) {
      w[i][j] = v[i];
    }
    assert(w[0][j] == x[j]);
  }

  std::vector<G1> com_w(s);
  std::vector<Fr> com_w_r(s);

  com_w[0] = com_x;
  com_w_r[0] = com_x_r;
  com_w[s - 1] = com_y;
  com_w_r[s - 1] = com_y_r;

  std::cout << "compute com(witness)\n";
  auto parallel_f = [&com_w_r, &com_w, &w](int64_t i) {
    com_w_r[i] = FrRand();
    com_w[i] = PcComputeCommitment(w[i], com_w_r[i]);
  };
  parallel::For(1LL, s - 1, parallel_f);

  parallel_r1cs::Prove(proof.hp_proof, seed, pb, w, com_w, com_w_r);
  proof.com_w = std::move(com_w);
}

inline bool Verify(Proof const& proof, h256_t seed, int64_t n, Fr const& k) {
  Tick tick(__FUNCTION__);
  int64_t const primary_input_size = 0;
  libsnark::protoboard<Fr> pb;
  MatchGadget gadget(pb, k);
  pb.set_input_sizes(primary_input_size);
  int64_t m = (int64_t)pb.num_constraints();
  int64_t s = (int64_t)pb.num_variables();
  if ((int64_t)proof.com_w.size() != s) {
    assert(false);
    return false;
  }
  if (proof.hp_proof.m() != m) {
    assert(false);
    return false;
  }
  std::vector<std::vector<Fr>> public_w;  // empty public_w
  return parallel_r1cs::Verify(proof.hp_proof, seed, n, pb, proof.com_w,
                               public_w);
}

inline bool Test() {
  auto seed = misc::RandH256();
  int64_t n = 100000;
  std::vector<Fr> x(n);
  FrRand(x);
  Fr k = x[(int64_t)rand()%n];
  std::vector<Fr> y(x.size());
  for (auto i = 0; i < x.size(); ++i) {
    y[i] = x[i] == k ? 1 : 0;
  }
  std::cout << "compute com(x), com(y)\n";
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitment(x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitment(y, com_y_r);

  Tick tick(__FUNCTION__);
  Proof proof;
  Prove(proof, seed, x, k, com_x, com_x_r, com_y, com_y_r);

  return Verify(proof, seed, n, k);
}
}  // namespace pc_utils::match