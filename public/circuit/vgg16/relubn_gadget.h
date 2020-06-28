#pragma once

#include "circuit/fixed_point/fixed_point.h"

namespace circuit::vgg16 {

// insure a >=-kFrD2N && a < kFrD2N, that is [-2^(D+2N), 2^(D+2N)-1]
// a' = a/2^N (reduce precision)
// b = a >= 0? a : 0
// c = alpha*(b-mu) + beta*2^(2N)
// ret = c/2^N (reduce precision)
// num_constraints:
// num_variables:

// must same
template <size_t D, size_t N>
class ReluBnGadget : public libsnark::gadget<Fr> {
  static_assert(D + N * 2 < 253, "invalid D,N");

 public:
  ReluBnGadget(libsnark::protoboard<Fr>& pb,
               const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    namespace fp = circuit::fp;
    // dnot change the order of the following variables allocation
    a_.allocate(this->pb, FMT(this->annotation_prefix, " a"));
    ret_.allocate(this->pb, FMT(this->annotation_prefix, " ret"));
    alpha_.allocate(this->pb, FMT(this->annotation_prefix, " alpha"));
    beta_.allocate(this->pb, FMT(this->annotation_prefix, " beta"));
    mu_.allocate(this->pb, FMT(this->annotation_prefix, " mu"));    

    b_.allocate(this->pb, FMT(this->annotation_prefix, " b"));
    c_.allocate(this->pb, FMT(this->annotation_prefix, " c"));

    precision_gadget1_.reset(new fp::PrecisionGadget<D, 2 * N, N>(
        this->pb, a_, FMT(this->annotation_prefix, " precision_gadget a")));

    precision_gadget2_.reset(new fp::PrecisionGadget<D, 2 * N, N>(
        this->pb, c_, FMT(this->annotation_prefix, " precision_gadget c")));    
    generate_r1cs_constraints();
  }

  void Assign(Fr const& a, Fr const& alpha, Fr const& beta, Fr const& mu) {
    this->pb.val(a_) = a;
    this->pb.val(alpha_) = alpha;
    this->pb.val(beta_) = beta;
    this->pb.val(mu_) = mu;
    generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

 private:
  void generate_r1cs_constraints() {
    // a' = a/(2^N)
    precision_gadget1_->generate_r1cs_constraints();
    auto a2 = precision_gadget1_->ret();
    auto sign = precision_gadget1_->sign();

    // b = sign_ * a'
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(a2, sign, b_),
        FMT(this->annotation_prefix, " b = sign? a':0"));

    // c = alpha*(b-mu) + beta*2^N
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(
            alpha_, b_ - mu_, c_ - beta_ * fp::RationalConst<D, N>().kFrN),
        FMT(this->annotation_prefix, " c = alpha*(b-mu) + beta*2^N"));

    // ret = c/(2^N)
    precision_gadget2_->generate_r1cs_constraints();

    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(precision_gadget2_->ret(), 1, ret_),
        FMT(this->annotation_prefix, " ret"));

    // precision_gadget2_->ret()
  }

  void generate_r1cs_witness() {
    namespace fp = circuit::fp;
    auto alpha = this->pb.val(alpha_);
    auto beta = this->pb.val(beta_);
    auto mu = this->pb.val(mu_);
    //std::cout << "gadget alpha:" << alpha << "\n";
    //std::cout << "gadget beta:" << beta << "\n";
    //std::cout << "gadget mu:" << mu << "\n";

    precision_gadget1_->generate_r1cs_witness();
    auto a = this->pb.val(a_);
    //std::cout << "gadget a: " << a << "\n";

    precision_gadget1_->ret().evaluate(this->pb);
    auto a2 = this->pb.lc_val(precision_gadget1_->ret());
    //std::cout << "gadget a2: " << a2 << "\n";

    precision_gadget1_->sign().evaluate(this->pb);
    auto sign = this->pb.lc_val(precision_gadget1_->sign());
    auto b = sign == Fr(0) ? Fr(0) : a2;
    this->pb.val(b_) = b;
    //std::cout << "gadget b: " << b << "\n";

    auto c = alpha * (b - mu) + beta * fp::RationalConst<D, N>().kFrN;
    this->pb.val(c_) = c;
    //std::cout << "gadget c: " << c << "\n";

    precision_gadget2_->generate_r1cs_witness();    
    precision_gadget2_->ret().evaluate(this->pb);
    auto ret = this->pb.lc_val(precision_gadget2_->ret());
    this->pb.val(ret_) = ret;
    //std::cout << "gadget ret: " << ret << "\n";
  }

 private:
  libsnark::pb_variable<Fr> a_;
  libsnark::pb_variable<Fr> b_;
  libsnark::pb_variable<Fr> c_;
  libsnark::pb_variable<Fr> alpha_;
  libsnark::pb_variable<Fr> beta_;
  libsnark::pb_variable<Fr> mu_;
  libsnark::pb_variable<Fr> ret_;
  std::unique_ptr<fixed_point::PrecisionGadget<D, N * 2, N>> precision_gadget1_;
  std::unique_ptr<fixed_point::PrecisionGadget<D, N * 2, N>> precision_gadget2_;
};

}  // namespace circuit::vgg16