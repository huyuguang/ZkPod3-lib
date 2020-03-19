#pragma once

#include "./div_gadget.h"

namespace circuit::fixed_point {

// ret = a / b, c = a / b ... d
// C * B + D * 2^N == A * 2^N && |D| < |B| && C is valid

template <size_t D, size_t N>
class Div2Gadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  Div2Gadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& a,
             libsnark::pb_linear_combination<Fr> const& sign_a,
             libsnark::pb_linear_combination<Fr> const& b,             
             libsnark::pb_linear_combination<Fr> const& sign_b,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        a_(a),
        sign_a_(sign_a),
        b_(b),
        sign_b_(sign_b) {
    abs_a_.allocate(pb, FMT(this->annotation_prefix, " abs_a"));
    abs_b_.allocate(pb, FMT(this->annotation_prefix, " abs_b"));
    sign_ret_.allocate(pb, FMT(this->annotation_prefix, " sign_ret"));
    ret_.allocate(pb, FMT(this->annotation_prefix, " ret"));

    div_gadget_.reset(
        new DivGadget<D, N>(this->pb, abs_a_, abs_b_, false,
                            FMT(this->annotation_prefix, " div_gadget")));
  }

  void generate_r1cs_constraints() {
    // abs_a = sign_a? a:-a
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_a_, a_ * 2, a_ + abs_a_),
        FMT(this->annotation_prefix, " abs_a = sign_a? a:-a"));
    // abs_b = sign_b? b:-b
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_b_, b_ * 2, b_ + abs_b_),
        FMT(this->annotation_prefix, " abs_b = sign_b? b:-b"));

    // sign_ret = 2 * sign_a * sign_b + 1 - sign_a - sign_b
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_a_, sign_b_ * 2,
                                      sign_a_ + sign_b_ - 1 + sign_ret_),
        FMT(this->annotation_prefix,
            " sign_ret = 2 * sign_a * sign_b + 1 - sign_a - sign_b"));

    // ret = sign_ret? positive_ret: -positive_ret
    auto positive_ret = div_gadget_->ret();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_ret_, positive_ret * 2,
                                      positive_ret + ret_),
        FMT(this->annotation_prefix,
            " ret = sign_ret? positive_ret: -positive_ret"));

    div_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);
    sign_a_.evaluate(this->pb);
    sign_b_.evaluate(this->pb);

    Fr a = this->pb.lc_val(a_);
    Fr sign_a = this->pb.lc_val(sign_a_);
    Fr check_sign_a = a.isNegative() ? 0 : 1;
    std::cout << "a: " << a << "\n";
    std::cout << "sign_a: " << sign_a << "\n";
    std::cout << "check_sign_a: " << check_sign_a << "\n";

    assert(sign_a == check_sign_a);
    this->pb.val(abs_a_) = a.isNegative() ? -a : a;

    Fr b = this->pb.lc_val(b_);
    Fr sign_b = this->pb.lc_val(sign_b_);
    Fr check_sign_b = b.isNegative() ? 0 : 1;
    std::cout << "b: " << b << "\n";
    std::cout << "sign_b: " << sign_b << "\n";
    std::cout << "check_sign_b: " << check_sign_b << "\n";

    assert(sign_b == check_sign_b);
    this->pb.val(abs_b_) = b.isNegative() ? -b : b;

    this->pb.val(sign_ret_) = sign_a * sign_b * 2 + 1 - sign_a - sign_b;

    div_gadget_->generate_r1cs_witness();

    Fr positive_ret = this->pb.val(div_gadget_->ret());
    this->pb.val(ret_) =
        this->pb.val(sign_ret_) == 1 ? positive_ret : -positive_ret;
    std::cout << "positive_ret: " << positive_ret << "\n";
    std::cout << "ret: " << this->pb.val(ret_) << "\n";
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

  libsnark::pb_variable<Fr> sign() const { return sign_ret_; }

  static bool Test(double double_a, double double_b);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> sign_a_;
  libsnark::pb_linear_combination<Fr> b_;
  libsnark::pb_linear_combination<Fr> sign_b_;
  libsnark::pb_variable<Fr> abs_a_;
  libsnark::pb_variable<Fr> abs_b_;
  libsnark::pb_variable<Fr> sign_ret_;
  libsnark::pb_variable<Fr> ret_;
  std::unique_ptr<DivGadget<D, N>> div_gadget_;
};

template <size_t D, size_t N>
bool Div2Gadget<D, N>::Test(double double_a, double double_b) {
  auto double_ret = double_a / double_b;
  std::cout << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<N>(double_a);
  Fr b = DoubleToRational<N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  libsnark::pb_variable<Fr> pb_sign_a;
  libsnark::pb_variable<Fr> pb_sign_b;
  pb_a.allocate(pb, "TestDiv2 a");
  pb_b.allocate(pb, "TestDiv2 b");
  pb_sign_a.allocate(pb, "TestDiv2 sign_a");
  pb_sign_b.allocate(pb, "TestDiv2 sign_b");
  Div2Gadget<D, N> gadget(pb, pb_a, pb_sign_a, pb_b, pb_sign_b, "TestDiv2");
  gadget.generate_r1cs_constraints();
  pb.val(pb_a) = a;
  pb.val(pb_b) = b;
  pb.val(pb_sign_a) = Fr(double_a >= 0 ? 1 : 0);
  pb.val(pb_sign_b) = Fr(double_b >= 0 ? 1 : 0);
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t" << fp::RationalToDouble<N>(fr_ret)
            << "\n";
  Fr fr_sign = pb.lc_val(gadget.sign());
  std::cout << "sign: " << fr_sign << "\n";
  assert(fr_sign == (double_ret >= 0 ? 1 : 0));

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestDiv2() {
  double double_a = -7.3;
  double double_b = -32.1;
  return Div2Gadget<32, 32>::Test(double_a, double_b);
}
}  // namespace circuit::fixed_point