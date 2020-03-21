#pragma once

#include "./div2_gadget.h"
#include "./sign_gadget.h"

namespace circuit::fixed_point {

// ret = a / b, c = a / b ... d
// C * B + D * 2^N == A * 2^N && |D| < |B| && C is valid

template <size_t D, size_t N>
class Div3Gadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  Div3Gadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& a,
             libsnark::pb_linear_combination<Fr> const& b,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),a_(a),b_(b) {
    sign_a_gadget_.reset(new SignGadget<D + N>(
        this->pb, a, FMT(this->annotation_prefix, " sign_a_gadget")));
    sign_b_gadget_.reset(new SignGadget<D + N>(
        this->pb, b, FMT(this->annotation_prefix, " sign_b_gadget")));

    div2_gadget_.reset(new Div2Gadget<D, N>(
        this->pb, a, sign_a_gadget_->ret(), b, sign_b_gadget_->ret(),
        FMT(this->annotation_prefix, " div2_gadget")));
  }

  void generate_r1cs_constraints() {
    sign_a_gadget_->generate_r1cs_constraints();
    sign_b_gadget_->generate_r1cs_constraints();
    div2_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);

    sign_a_gadget_->generate_r1cs_witness();
    sign_b_gadget_->generate_r1cs_witness();

    div2_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return div2_gadget_->ret(); }
  
  libsnark::pb_variable<Fr> sign() const { return div2_gadget_->sign(); }

  static bool Test(double double_a, double double_b);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> b_;
  std::unique_ptr<SignGadget<D + N>> sign_a_gadget_;
  std::unique_ptr<SignGadget<D + N>> sign_b_gadget_;
  std::unique_ptr<Div2Gadget<D, N>> div2_gadget_;
};

template <size_t D, size_t N>
bool Div3Gadget<D, N>::Test(double double_a, double double_b) {
  auto double_ret = double_a / double_b;
  std::cout << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<D, N>(double_a);
  Fr b = DoubleToRational<D, N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  pb_a.allocate(pb, "TestDiv a");
  pb_b.allocate(pb, "TestDiv b");
  Div3Gadget<D, N> gadget(pb, pb_a, pb_b, "Test3Div");
  gadget.generate_r1cs_constraints();
  pb.val(pb_a) = a;
  pb.val(pb_b) = b;
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t" << fp::RationalToDouble<D, N>(fr_ret)
            << "\n";
  Fr fr_sign = pb.lc_val(gadget.sign());
  std::cout << "sign: " << fr_sign << "\n";
  assert(fr_sign == (double_ret >= 0 ? 1 : 0));

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestDiv3() {
  Tick tick(__FN__);
  double double_a = -7.3;
  double double_b = -32.1;
  return Div3Gadget<32, 32>::Test(double_a, double_b);
}
}  // namespace circuit::fixed_point