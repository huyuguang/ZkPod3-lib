#pragma once

#include "./sign_gadget.h"

namespace circuit::fixed_point {

// insure x >=-kFrDN && x < kFrDN, that is [-2^(D+N), 2^(D+N)-1]
// ret = max(a)
// assert(sum(sign(a-x[i]))==x.size())
// assert((a-x1)(a-x2)...==0)
// num_constraints:
// num_variables:

template <size_t D, size_t N>
class MaxGadget : public libsnark::gadget<Fr> {
  static_assert(D + N < 253, "invalid D,N");

 public:
  MaxGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination_array<Fr> const& x,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        x_(x),
        diffs_(x.size()) {
    assert(x.size() > 1);

    max_.allocate(this->pb, " max");

    for (size_t i = 0; i < x_.size(); ++i) {
      diffs_[i].assign(this->pb, max_ - x_[i]);
    }

    sign_gadgets_.resize(x_.size());
    for (size_t i = 0; i < sign_gadgets_.size(); ++i) {
      sign_gadgets_[i].reset(new SignGadget<D, N>(
          this->pb, diffs_[i],
          FMT(this->annotation_prefix, " sign of max - x[%uz]", i)));
    }
    diff_prods_.allocate(this->pb, x.size() - 1, " diff_prods");
  }

  void generate_r1cs_constraints() {
    libsnark::linear_combination<Fr> sum_of_sign;
    for (size_t i = 0; i < sign_gadgets_.size(); ++i) {
      sign_gadgets_[i]->generate_r1cs_constraints();
      sum_of_sign = sum_of_sign + sign_gadgets_[i]->ret();
    }
    sum_of_sign = sum_of_sign - Fr(x_.size());

    libsnark::pb_linear_combination<Fr> lc_sum_of_sign;
    lc_sum_of_sign.assign(this->pb, sum_of_sign);
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(lc_sum_of_sign, Fr(1), Fr(0)),
        FMT(this->annotation_prefix, " sum_of_sign = x.size"));

    // (diffs_[0]) * (diffs_[1]) == diff_prods_[0]
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(diffs_[0], diffs_[1], diff_prods_[0]),
        " diffs_[0] * diffs_[1 == diff_prods_[0]");
    for (size_t i = 1; i < x_.size() - 1; ++i) {
      // (diff_prods_[i-1]) * (diffs_[i+1]) == diff_prods_[i]
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(diff_prods_[i - 1], diffs_[i + 1],
                                        diff_prods_[i]),
          " diff_prods_[i-1] * diffs_[i+1] == diff_prods_[i]");
    }
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(diff_prods_[diff_prods_.size() - 1],
                                      Fr(1), Fr(0)),
        FMT(this->annotation_prefix, " diff_prods.back = 0"));
  }

  void generate_r1cs_witness() {
    for (size_t i = 0; i < x_.size(); ++i) {
      x_[i].evaluate(this->pb);
    }
    Fr max = this->pb.lc_val(x_[0]);
    for (size_t i = 1; i < x_.size(); ++i) {
      Fr x = this->pb.lc_val(x_[i]);
      Fr diff = max - x;
      if (diff.isNegative()) max = x;
    }

    this->pb.val(max_) = max;

    diffs_.evaluate(this->pb);

    for (size_t i = 0; i < sign_gadgets_.size(); ++i) {
      sign_gadgets_[i]->generate_r1cs_witness();
    }

    // diff_prods_[0] = diffs_[0] * diffs_[1];
    this->pb.val(diff_prods_[0]) =
        this->pb.lc_val(diffs_[0]) * this->pb.lc_val(diffs_[1]);
    for (size_t i = 1; i < x_.size() - 1; ++i) {
      // diff_prods_[i] = diff_prods_[i-1] * diffs_[i+1]
      this->pb.val(diff_prods_[i]) =
          this->pb.val(diff_prods_[i - 1]) * this->pb.lc_val(diffs_[i + 1]);
    }
  }

  libsnark::pb_variable<Fr> ret() const { return max_; }

  static bool Test(std::vector<Fr> const& x, Fr max) {
    libsnark::protoboard<Fr> pb;
    libsnark::pb_variable_array<Fr> pb_x;
    pb_x.allocate(pb, x.size(), "Test x");

    MaxGadget<D, N> gadget(pb, pb_x, "MaxGadget");

    gadget.generate_r1cs_constraints();
    for (size_t i = 0; i < x.size(); ++i) {
      pb.val(pb_x[i]) = x[i];
    }

    gadget.generate_r1cs_witness();
    assert(pb.is_satisfied());
    assert(pb.val(gadget.ret()) == max);
    return pb.is_satisfied() && pb.val(gadget.ret()) == max;
  }

 private:
  libsnark::pb_linear_combination_array<Fr> x_;
  libsnark::pb_variable<Fr> max_;
  libsnark::pb_linear_combination_array<Fr> diffs_;
  libsnark::pb_variable_array<Fr> diff_prods_;
  std::vector<std::unique_ptr<SignGadget<D, N>>> sign_gadgets_;
};

inline bool TestMax() {
  Tick tick(__FN__);
  constexpr size_t D = 8;
  constexpr size_t N = 2;
  std::vector<Fr> x;
  std::vector<bool> rets;

  rets.push_back(MaxGadget<D, N>::Test({100, 0, -11}, 100));
  rets.push_back(MaxGadget<D, N>::Test({1, 3, -11}, 3));
  rets.push_back(MaxGadget<D, N>::Test({-32, 31, -1}, 31));
  rets.push_back(MaxGadget<D, N>::Test({-32, 31, -1}, 31));
  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point