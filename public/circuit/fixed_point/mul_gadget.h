#pragma once

#include "./misc.h"

namespace circuit::fixed_point {

// ret = a * b / 2^N
template <size_t D, size_t N>
class MulGadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  MulGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& a,
            libsnark::pb_linear_combination<Fr> const& b,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a), b_(b) {
    auto constances = RationalConst<D, N>();

    product_.allocate(pb, FMT(this->annotation_prefix, " product"));

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
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(a_, b_, product_),
        FMT(this->annotation_prefix, " a * b = product"));
    p1_gadget_->generate_r1cs_constraints(true);
    p2_gadget_->generate_r1cs_constraints(false);
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);
    Fr a = this->pb.lc_val(a_);
    Fr b = this->pb.lc_val(b_);
    Fr c = a * b;
    this->pb.val(product_) = c;

    //std::cout << "mul a: " << a << "\t" << RationalToDouble<D,N>(a) <<"\n";
    //std::cout << "mul b: " << b << "\t" << RationalToDouble<D,N>(b) <<"\n";
    //std::cout << "mul c: " << c << "\t" << RationalToDouble<D,N>(c) <<"\n";

    product_off_.evaluate(this->pb);
    p1_gadget_->generate_r1cs_witness_from_packed();
    p2_gadget_->generate_r1cs_witness_from_bits();
    ret_.evaluate(this->pb);
    sign_.evaluate(this->pb);
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

  // 1: >=0; 0: <=0
  libsnark::pb_linear_combination<Fr> sign() const { return sign_; };

  static bool Test(double double_a, double double_b);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> b_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> sign_;
  libsnark::pb_variable<Fr> product_;
  libsnark::pb_linear_combination<Fr> product_off_;  // = product_ + 2^(D+2N)
  std::unique_ptr<libsnark::packing_gadget<Fr>> p1_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p2_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
  libsnark::pb_variable<Fr> packed_;
};

template <size_t D, size_t N>
bool MulGadget<D, N>::Test(double double_a, double double_b) {
  auto double_ret = double_a * double_b;
  std::cout << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<D, N>(double_a);
  Fr b = DoubleToRational<D, N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  pb_a.allocate(pb, "TestMul a");
  pb_b.allocate(pb, "TestMul b");
  MulGadget<D, N> gadget(pb, pb_a, pb_b, "TestMul");
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

inline bool TestMul() {
  Tick tick(__FN__);
  double double_a = 7.3;
  double double_b = 32.1;
  return MulGadget<32, 32>::Test(double_a, double_b);
}
}  // namespace circuit::fixed_point