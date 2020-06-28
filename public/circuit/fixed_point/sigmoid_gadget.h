#pragma once

#include "./exp_gadget.h"

namespace circuit::fixed_point {

// sigmoid(x) = {\frac {1}{1+e^{-x}}}

template <size_t D, size_t N>
class SigmoidGadget : public libsnark::gadget<Fr> {
 public:
  SigmoidGadget(libsnark::protoboard<Fr>& pb,
                libsnark::pb_linear_combination<Fr> const& x,
                size_t level, const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), x_(x) {

    neg_x_.assign(this->pb, -x);

    exp_gadget_.reset(new ExpGadget<D, N>(
        pb, neg_x_, level, FMT(this->annotation_prefix, " exp_gadget")));
    exp_plus_1_.assign(this->pb,
                       exp_gadget_->ret() + RationalConst<D, N>().kFrN);

    // since e^x always >=0, check_sign = false
    libsnark::pb_linear_combination<Fr> pb_lc_sign;
    pb_lc_sign.assign(pb, FrOne());  // exp_plus_1 always > 0
    inv_gadget_.reset(
        new InvGadget<D, N>(pb, exp_plus_1_, &pb_lc_sign,
                                FMT(this->annotation_prefix, " inv_gadget")));
  }

  void generate_r1cs_constraints() {
    exp_gadget_->generate_r1cs_constraints();
    inv_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    neg_x_.evaluate(this->pb);

    exp_gadget_->generate_r1cs_witness();

    exp_plus_1_.evaluate(this->pb);
    inv_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_linear_combination<Fr> ret() const { return inv_gadget_->ret(); }

  static bool Test(double const& double_x);

 private:
  libsnark::pb_linear_combination<Fr> x_;
  std::unique_ptr<ExpGadget<D, N>> exp_gadget_;
  std::unique_ptr<InvGadget<D, N>> inv_gadget_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> neg_x_;
  libsnark::pb_linear_combination<Fr> exp_plus_1_;
};

template <size_t D, size_t N>
bool SigmoidGadget<D, N>::Test(double const& double_x) {
  auto double_ret = 1.0 / (1.0 + std::exp(-double_x));
  std::cout << "double_ret: " << double_ret << "\n";

  Fr x = DoubleToRational<D, N>(double_x);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_x;
  pb_x.allocate(pb, "SigmoidGadget::Test x");

  SigmoidGadget<D, N> gadget(pb, pb_x, 10, "SigmoidGadget::Test");
  gadget.generate_r1cs_constraints();
  pb.val(pb_x) = x;
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t"
            << RationalToDouble<D, N>(fr_ret) << "\n";
  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline static bool TestSigmoid() {
  Tick tick(__FN__);
  bool ret;
  std::vector<bool> rets;
  double double_x = -1.3;
  ret = SigmoidGadget<32, 32>::Test(double_x);
  rets.push_back(ret);

  double_x = 0.5;
  ret = SigmoidGadget<32, 32>::Test(double_x);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit
