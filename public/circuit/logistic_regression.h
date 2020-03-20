#pragma once

#include "./fixed_point/fixed_point.h"

namespace circuit {

// h(theta,x) = 1/(1+e^-<theta,x>)

template <size_t D, size_t N>
class LrFuncHGadget : public libsnark::gadget<Fr> {
 public:
  LrFuncHGadget(libsnark::protoboard<Fr>& pb,
                libsnark::pb_linear_combination_array<Fr> const& theta,
                libsnark::pb_linear_combination_array<Fr> const& x,
                size_t level, const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), theta_(theta), x_(x) {
    ip_gadget_.reset(new fp::IpGadget<D, N>(
        pb, theta, x, FMT(this->annotation_prefix, " ip_gadget")));

    neg_ip_.assign(this->pb, -ip_gadget_->ret());    

    exp_gadget_.reset(
        new fp::Exp2Gadget<D, N>(pb, neg_ip_, level,
                                FMT(this->annotation_prefix, " exp_gadget")));    
    exp_plus_1_.assign(this->pb,
                      exp_gadget_->ret() + fp::RationalConst<D, N>().kFrN);

    // since e^x always >=0, check_sign = false
    inv_gadget_.reset(new fp::InvGadget<D, N>(
        pb, exp_plus_1_, FMT(this->annotation_prefix, " inv_gadget")));
  }

  void generate_r1cs_constraints() {
    ip_gadget_->generate_r1cs_constraints();
    exp_gadget_->generate_r1cs_constraints();
    inv_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    ip_gadget_->generate_r1cs_witness();
    neg_ip_.evaluate(this->pb);

    exp_gadget_->generate_r1cs_witness();

    exp_plus_1_.evaluate(this->pb);
    inv_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return inv_gadget_->ret(); }

  static bool Test(std::vector<double> const& double_theta,
                   std::vector<double> const& double_x);

 private:
  libsnark::pb_linear_combination_array<Fr> theta_;
  libsnark::pb_linear_combination_array<Fr> x_;
  std::unique_ptr<fp::IpGadget<D, N>> ip_gadget_;
  std::unique_ptr<fp::Exp2Gadget<D, N>> exp_gadget_;
  std::unique_ptr<fp::InvGadget<D, N>> inv_gadget_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> neg_ip_;
  libsnark::pb_linear_combination<Fr> exp_plus_1_;
};

template <size_t D, size_t N>
bool LrFuncHGadget<D, N>::Test(std::vector<double> const& double_theta,
                               std::vector<double> const& double_x) {
  assert(double_theta.size() == double_x.size());
  size_t count = double_theta.size();
  auto double_ip = std::inner_product(double_theta.begin(), double_theta.end(),
                                      double_x.begin(), 0.0);
  auto double_ret = 1.0 / (1.0 + std::exp(-double_ip));
  std::cout << "double_ret: " << double_ret << "\n";

  std::vector<Fr> theta(count);
  std::vector<Fr> x(count);
  for (size_t i = 0; i < count; ++i) {
    theta[i] = fp::DoubleToRational<D, N>(double_theta[i]);
    x[i] = fp::DoubleToRational<D, N>(double_x[i]);
  }

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable_array<Fr> pb_theta;
  libsnark::pb_variable_array<Fr> pb_x;
  pb_theta.allocate(pb, count, "LrFuncHGadget::Test theta");
  pb_x.allocate(pb, count, "LrFuncHGadget::Test x");

  LrFuncHGadget<D, N> gadget(pb, pb_theta, pb_x, 10, "LrFuncHGadget::Test");
  gadget.generate_r1cs_constraints();
  pb_theta.fill_with_field_elements(pb, theta);
  pb_x.fill_with_field_elements(pb, x);
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t" << fp::RationalToDouble<D, N>(fr_ret)
            << "\n";
  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline static bool TestFuncH() {
  Tick tick(__FN__);
  std::vector<double> double_theta{1, 2, 3};
  std::vector<double> double_x{3, -5, 2};
  return LrFuncHGadget<32, 32>::Test(double_theta, double_x);
}
}  // namespace circuit
