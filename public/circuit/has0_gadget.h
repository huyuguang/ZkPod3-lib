#pragma once

#include <stdlib.h>
#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>
#include <vector>

#include "ecc/ecc.h"

namespace circuit {
  // if (is_any(x,0)) ret = 1; else ret = 0;
  class HasZeroGadget : public libsnark::gadget<Fr> {
   public:
    HasZeroGadget(libsnark::protoboard<Fr>& pb,
                  std::vector<libsnark::pb_linear_combination<Fr>> const& x)
        : libsnark::gadget<Fr>(pb, "HasZeroGadget"), x_(x) {
      assert(!x_.empty());
      y_.resize(x_.size() - 1);
      for (auto& i : y_) i.allocate(pb, " y");
      inv_.allocate(pb, " inv");
      ret_.allocate(pb, " ret");
    }

    libsnark::pb_variable<Fr> const& ret() const { return ret_; }

    // ret = (x0*x1...)? 0 : 1
    void generate_r1cs_constraints() {
      libsnark::pb_linear_combination<Fr> pie_of_x;
      if (x_.size() == 1) {
        pie_of_x = x_.back();
      } else {
        // x[0] * x[1] == y_[0]
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(x_[0], x_[1], y_[0]),
            "y[0] = x[0] * x[1]");
        for (size_t i = 1; i < x_.size() - 1; ++i) {
          // y[i-1] * x[i+1] == y[i]
          this->pb.add_r1cs_constraint(
              libsnark::r1cs_constraint<Fr>(y_[i - 1], x_[i + 1], y_[i]),
              "y[i] = y[i-1] * x[i+1]");
        }

        pie_of_x = y_.back();
      }

      // ret * pie_of_x == 0
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(ret_, pie_of_x, FrZero()),
          "ret_ * pie_of_x == 0");
      // inv * pie_of_x == 1 - ret
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(inv_, pie_of_x, -ret_ + FrOne()),
          "inv * pie_of_x == 1 - ret");
    }

    void generate_r1cs_witness() {
      Fr pie_of_x;
      if (x_.size() == 1) {
        // y_ is empty
        pie_of_x = this->pb.lc_val(x_[0]);
      } else {
        // y_[0] = x_[0] * x_[1];
        this->pb.val(y_[0]) = this->pb.lc_val(x_[0]) * this->pb.lc_val(x_[1]);
        for (size_t i = 1; i < x_.size() - 1; ++i) {
          // y[i] = y[i-1] * x[i+1]
          this->pb.val(y_[i]) =
              this->pb.val(y_[i - 1]) * this->pb.lc_val(x_[i + 1]);
        }

        pie_of_x = this->pb.val(y_.back());
      }
      this->pb.val(inv_) = pie_of_x == FrZero() ? FrZero() : FrInv(pie_of_x);
      this->pb.val(ret_) = pie_of_x == FrZero() ? FrOne() : FrZero();
    }

   private:
    std::vector<libsnark::pb_linear_combination<Fr>> const& x_;
    std::vector<libsnark::pb_variable<Fr>> y_;  // (x1*x2...)
    libsnark::pb_variable<Fr> inv_;             // y.back()? 1/y.back():0
    libsnark::pb_variable<Fr> ret_;             // y.back() == 0? 1 : 0
  };
}