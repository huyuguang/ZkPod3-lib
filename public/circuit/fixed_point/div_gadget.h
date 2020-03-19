#pragma once

#include "./abs_gadget.h"
#include "./type_gadget.h"

namespace circuit::fixed_point {

// ret = a / b, c = a / b ... d
// C * B + D * 2^N == A * 2^N && |D| < |B| && C is valid

template <size_t D, size_t N>
class DivGadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  DivGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& a,
            libsnark::pb_linear_combination<Fr> const& b, bool check_sign,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        a_(a),
        b_(b),
        check_sign_(check_sign) {
    auto constances = RationalConst<D, N>();
    c_.allocate(pb, FMT(this->annotation_prefix, " c"));
    dn_.allocate(pb, FMT(this->annotation_prefix, " dn"));

    libsnark::pb_linear_combination<Fr> bn;
    bn.assign(this->pb, b_ * constances.kFrN);

    if (check_sign_) {
      abs_bn_gadget_.reset(new AbsGadget<D + N + N>(
          this->pb, bn, FMT(this->annotation_prefix, " abs_bn")));
      abs_dn_gadget_.reset(new AbsGadget<D + N + N>(
          this->pb, dn_, FMT(this->annotation_prefix, " abs_dn")));
    }

    libsnark::pb_linear_combination<Fr> abs_bn;
    abs_bn.assign(this->pb, check_sign_ ? abs_bn_gadget_->ret() : bn);

    libsnark::pb_linear_combination<Fr> abs_dn;
    abs_dn.assign(this->pb, check_sign_ ? abs_dn_gadget_->ret() : dn_);

    less_.allocate(pb, FMT(this->annotation_prefix, " less"));
    less_or_eq_.allocate(pb, FMT(this->annotation_prefix, " less_or_eq"));
    comparison_gadget_.reset(new libsnark::comparison_gadget<Fr>(
        pb, D + N + N, abs_dn, abs_bn, less_, less_or_eq_,
        FMT(this->annotation_prefix, " comparison")));

    // check if c is valid
    sign_gadget_.reset(new SignGadget<D + N>(
        pb, ret(), FMT(this->annotation_prefix, " type_gadget")));
  }

  void generate_r1cs_constraints() {
    auto constances = RationalConst<D, N>();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(c_, b_, a_ * constances.kFrN - dn_),
        FMT(this->annotation_prefix, " C * B + D * 2^N == A * 2^N"));
    if (check_sign_) {
      abs_bn_gadget_->generate_r1cs_constraints();
      abs_dn_gadget_->generate_r1cs_constraints();
    }
    comparison_gadget_->generate_r1cs_constraints();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(less_, libsnark::pb_variable<Fr>(0),
                                      libsnark::pb_variable<Fr>(0)),
        FMT(this->annotation_prefix, " less == 1"));
    sign_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    auto constances = RationalConst<D, N>();
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);

    mpz_class mpz_a = FrToSignedMpz(this->pb.lc_val(a_));    
    mpz_class mpz_b = FrToSignedMpz(this->pb.lc_val(b_));
    mpz_class mpz_c = mpz_a * constances.kMpzN / mpz_b;
    mpz_class mpz_dn = mpz_a * constances.kMpzN - mpz_b * mpz_c;
    this->pb.val(c_) = SignedMpzToFr(mpz_c);
    this->pb.val(dn_) = SignedMpzToFr(mpz_dn);

    if (check_sign_) {
      abs_bn_gadget_->generate_r1cs_witness();
      abs_dn_gadget_->generate_r1cs_witness();
    }

    comparison_gadget_->generate_r1cs_witness();

    sign_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return c_; }
  
  libsnark::pb_variable<Fr> sign() const { return sign_gadget_->ret(); }

  static bool Test(double double_a, double double_b);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> b_;
  bool check_sign_;
  libsnark::pb_variable<Fr> c_;
  libsnark::pb_variable<Fr> dn_;
  std::unique_ptr<AbsGadget<D + N + N>> abs_bn_gadget_;
  std::unique_ptr<AbsGadget<D + N + N>> abs_dn_gadget_;
  libsnark::pb_variable<Fr> less_;
  libsnark::pb_variable<Fr> less_or_eq_;
  std::unique_ptr<libsnark::comparison_gadget<Fr>> comparison_gadget_;
  std::unique_ptr<SignGadget<D + N>> sign_gadget_;
};

template <size_t D, size_t N>
bool DivGadget<D, N>::Test(double double_a, double double_b) {
  auto double_ret = double_a / double_b;
  std::cout << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<N>(double_a);
  Fr b = DoubleToRational<N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  pb_a.allocate(pb, "TestDiv a");
  pb_b.allocate(pb, "TestDiv b");
  DivGadget<D, N> gadget(pb, pb_a, pb_b, true, "TestDiv");
  gadget.generate_r1cs_constraints();
  pb.val(pb_a) = a;
  pb.val(pb_b) = b;
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

inline bool TestDiv() {
  double double_a = -7.3;
  double double_b = -32.1;
  return DivGadget<32, 32>::Test(double_a, double_b);
}
}  // namespace circuit::fixed_point