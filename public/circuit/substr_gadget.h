#pragma once

#include "./has0_gadget.h"

namespace circuit {
class SubstrGadget : public libsnark::gadget<Fr> {
 public:
  SubstrGadget(libsnark::protoboard<Fr>& pb, std::string const& k)
      : libsnark::gadget<Fr>(pb, "SubstrGadget"), k_(PackStrToFr(k.c_str())) {
    // NOTE: 1 Fr pack 31 char
    x_.allocate(pb, " x");
    dual_.reset(
        new libsnark::dual_variable_gadget<Fr>(this->pb, x_, 31 * 8, "dual"));

    int64_t k_len = k.size();
    sub_duals_.resize(31 - k_len + 1);
    sub_minus_k_.resize(sub_duals_.size());
    for (size_t i = 0; i < sub_duals_.size(); ++i) {
      auto begin = dual_->bits.begin() + i * 8;
      auto end = begin + k_len * 8;
      libsnark::pb_variable_array<Fr> bits(begin, end);
      sub_duals_[i].reset(
          new libsnark::dual_variable_gadget<Fr>(this->pb, bits, "sub_dual"));
      sub_minus_k_[i].assign(this->pb, sub_duals_[i]->packed - k_);
    }

    ret_.reset(new HasZeroGadget(pb, sub_minus_k_));

    generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> const& ret() const { return ret_->ret(); }

  void Assign(Fr const& x) {
    this->pb.val(x_) = x;
    generate_r1cs_witness();
    //assert(this->pb.is_satisfied());
  }

 private:
  void generate_r1cs_constraints() {
    dual_->generate_r1cs_constraints(true);

    for (auto& i : sub_duals_) {
      i->generate_r1cs_constraints(false);
    }

    ret_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    dual_->generate_r1cs_witness_from_packed();

    for (auto& i : sub_duals_) {
      i->generate_r1cs_witness_from_bits();
    }

    for (auto& i : sub_minus_k_) {
      i.evaluate(this->pb);
    }

    ret_->generate_r1cs_witness();
  }

 private:
  Fr k_;

  libsnark::pb_variable<Fr> x_;

  // 248 var, 0 or 1
  std::unique_ptr<libsnark::dual_variable_gadget<Fr>> dual_;

  // 31 - k_len + 1 var
  std::vector<std::unique_ptr<libsnark::dual_variable_gadget<Fr>>> sub_duals_;

  std::vector<libsnark::pb_linear_combination<Fr>> sub_minus_k_;

  std::unique_ptr<HasZeroGadget> ret_;
};
}  // namespace circuit