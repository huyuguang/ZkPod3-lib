#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./mimc5.h"

namespace circuit {
class TestGadget : public libsnark::gadget<Fr> {
 public:
  TestGadget(libsnark::protoboard<Fr>& pb, const std::string& annotation_prefix,
             size_t m)
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    x_.allocate(pb, m, FMT(annotation_prefix, " rounds_x2"));

    generate_r1cs_constraints();
  }

  void Assign(Fr const& x) { generate_r1cs_witness(x); }

 private:
  void generate_r1cs_constraints() {
    for (size_t i = 1; i < x_.size(); ++i) {
      std::string tag =
          "x" + std::to_string(i - 1) + "^2=x" + std::to_string(i);
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(x_[i - 1], x_[i - 1], x_[i]), "xx");
    }
  }

  void generate_r1cs_witness(Fr const& x) {
    this->pb.val(x_[0]) = x;
    for (size_t i = 1; i < x_.size(); ++i) {
      this->pb.val(x_[i]) = this->pb.val(x_[i - 1]) * this->pb.val(x_[i - 1]);
    }
  }

 private:
  libsnark::pb_variable_array<Fr> x_;
};

}  // namespace circuit
