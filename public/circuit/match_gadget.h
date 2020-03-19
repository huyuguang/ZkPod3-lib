#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>
#include <vector>

#include "ecc/ecc.h"

namespace circuit {
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
    //assert(this->pb.is_satisfied());
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
  libsnark::pb_variable<Fr> x2_;  // k==x? 0:1/(x-k)
  libsnark::pb_variable<Fr> y_;
};

}  // namespace circuit