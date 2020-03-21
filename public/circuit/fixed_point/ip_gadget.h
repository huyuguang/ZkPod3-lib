#pragma once

#include "./misc.h"

namespace circuit::fixed_point {

// ret = <a, b>
template <size_t D, size_t N>
class IpGadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  IpGadget(libsnark::protoboard<Fr>& pb,
           libsnark::pb_linear_combination_array<Fr> const& a,
           libsnark::pb_linear_combination_array<Fr> const& b,
           const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    product_.allocate(pb, FMT(this->annotation_prefix, " product"));
    ip_gadget_.reset(new libsnark::inner_product_gadget<Fr>(
        this->pb, a, b, product_, FMT(this->annotation_prefix, " ip")));

    auto constances = RationalConst<D, N>();
    product_off_.assign(this->pb, libsnark::linear_combination<Fr>(product_) +
                                      constances.kFrD2N);

    bits_.allocate(this->pb, D + 2 * N + 1,
                   FMT(this->annotation_prefix, " bits"));
    p1_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb, bits_, product_off_, FMT(this->annotation_prefix, " p1")));

    packed_.allocate(pb, FMT(this->annotation_prefix, " packed"));
    p2_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb,
        libsnark::pb_variable_array<Fr>(bits_.begin() + N, bits_.end()),
        packed_, FMT(this->annotation_prefix, " p2")));

    ret_.assign(this->pb,
                libsnark::linear_combination<Fr>(packed_) - constances.kFrDN);
    sign_.assign(this->pb, p2_gadget_->bits[D + N]);
  }

  void generate_r1cs_constraints() {
    ip_gadget_->generate_r1cs_constraints();
    p1_gadget_->generate_r1cs_constraints(true);
    p2_gadget_->generate_r1cs_constraints(false);
  }

  void generate_r1cs_witness() {
    ip_gadget_->generate_r1cs_witness();
    product_off_.evaluate(this->pb);
    p1_gadget_->generate_r1cs_witness_from_packed();
    p2_gadget_->generate_r1cs_witness_from_bits();
    ret_.evaluate(this->pb);
    sign_.evaluate(this->pb);
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

  // 1: >=0; 0: <=0
  libsnark::pb_linear_combination<Fr> sign() const { return sign_; };

  static bool Test(std::vector<double> const& double_a,
                   std::vector<double> const& double_b);

 private:
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> sign_;
  libsnark::pb_variable<Fr> product_;
  libsnark::pb_linear_combination<Fr> product_off_;  // = product_ + 2^(D+2N)
  std::unique_ptr<libsnark::inner_product_gadget<Fr>> ip_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p1_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p2_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
  libsnark::pb_variable<Fr> packed_;
};

template <size_t D, size_t N>
bool IpGadget<D, N>::Test(std::vector<double> const& double_a,
                          std::vector<double> const& double_b) {
  assert(double_a.size() == double_b.size());
  size_t count = double_a.size();
  auto double_ret = std::inner_product(double_a.begin(), double_a.end(),
                                       double_b.begin(), 0.0);
  std::cout << "double_ret: " << double_ret << "\n";

  std::vector<Fr> a(count);
  std::vector<Fr> b(count);
  for (size_t i = 0; i < count; ++i) {
    a[i] = DoubleToRational<D, N>(double_a[i]);
    b[i] = DoubleToRational<D, N>(double_b[i]);
  }

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable_array<Fr> pb_a;
  libsnark::pb_variable_array<Fr> pb_b;
  pb_a.allocate(pb, count, "TestIp a");
  pb_b.allocate(pb, count, "TestIp b");
  IpGadget<D, N> gadget(pb, pb_a, pb_b, "TestIp");
  gadget.generate_r1cs_constraints();
  pb_a.fill_with_field_elements(pb, a);
  pb_b.fill_with_field_elements(pb, b);
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t"
            << fp::RationalToDouble<D, N>(fr_ret) << "\n";
  Fr fr_sign = pb.lc_val(gadget.sign());
  std::cout << "sign: " << fr_sign << "\n";
  assert(fr_sign == (double_ret >= 0 ? 1 : 0));

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestIp() {
  Tick tick(__FN__);
  bool ret;
  std::vector<bool> rets;
  std::vector<double> double_a;
  std::vector<double> double_b;

  double_a = std::vector<double> {0.1, -2.2, 3};
  double_b = std::vector<double>{3, 5, 2.7};
  ret =  IpGadget<32, 32>::Test(double_a, double_b);
  rets.push_back(ret);

  double_a = std::vector<double> {10.1};
  double_b = std::vector<double>{0};
  ret =  IpGadget<32, 32>::Test(double_a, double_b);
  rets.push_back(ret);

  double_a = std::vector<double> {10.1};
  double_b = std::vector<double>{0.002};
  ret =  IpGadget<32, 32>::Test(double_a, double_b);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point