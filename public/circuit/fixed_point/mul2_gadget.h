#pragma once

#include "./abs_gadget.h"

namespace circuit::fixed_point {

// ret = a[0]*a[1]...*c, c is a const
template <size_t D, size_t N>
class Mul2Gadget : public libsnark::gadget<Fr> {
 public:
  Mul2Gadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination_array<Fr> const& a,
             Fr const& c, const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a), c_(c) {
    CheckMaxNumOfBits();

    auto constances = RationalConst<D, N>();

    products_.allocate(pb, a_.size() - 1, FMT(this->annotation_prefix, " products"));

    Fr kFrDxN = constances.GetFrDxN(a_.size() + 1);
    product_off_.assign(this->pb, libsnark::linear_combination<Fr>(
                                      products_[products_.size() - 1]) +
                                      kFrDxN);

    bits_.allocate(this->pb, D + (a_.size() + 1) * N + 1,
                   FMT(this->annotation_prefix, " bits"));
    p1_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb, bits_, product_off_, FMT(this->annotation_prefix, " p1")));

    packed_.allocate(pb, FMT(this->annotation_prefix, " packed"));
    p2_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb,
        libsnark::pb_variable_array<Fr>(bits_.begin() + a_.size() * N,
                                        bits_.end()),
        packed_, FMT(this->annotation_prefix, " p2")));

    ret_.assign(this->pb,
                libsnark::linear_combination<Fr>(packed_) - constances.kFrDN);
    sign_.assign(this->pb, p2_gadget_->bits[D + N]);
  }

  void generate_r1cs_constraints() {
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(a_[0], a_[1] * c_, products_[0]),
        FMT(this->annotation_prefix, " a[0]*a[1]*c = products[0]"));
    for (size_t i = 2; i < a_.size(); ++i) {
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(a_[i], products_[i - 2],
                                        products_[i - 1]),
          FMT(this->annotation_prefix, " a[%uz]*a[%zu]*c = products[%zu]", i,
              i - 2, i - 1));
    }

    p1_gadget_->generate_r1cs_constraints(true);
    p2_gadget_->generate_r1cs_constraints(false);
  }

  void generate_r1cs_witness() {
    for (auto& i : a_) {
      i.evaluate(this->pb);
    }

    this->pb.val(products_[0]) =
        this->pb.lc_val(a_[0]) * this->pb.lc_val(a_[1]) * c_;

    for (size_t i = 2; i < a_.size(); ++i) {
      this->pb.val(products_[i - 1]) =
          this->pb.lc_val(a_[i]) * this->pb.val(products_[i - 2]);
    }

    product_off_.evaluate(this->pb);

    for (auto& i : products_) {
      Fr f = this->pb.val(i);
    }

    Fr fr_product_off = this->pb.lc_val(product_off_);
    p1_gadget_->generate_r1cs_witness_from_packed();
    p2_gadget_->generate_r1cs_witness_from_bits();
    ret_.evaluate(this->pb);
    sign_.evaluate(this->pb);
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

  // 1: >=0; 0: <=0
  libsnark::pb_linear_combination<Fr> sign() const { return sign_; };

  static bool Test(std::vector<double> const& double_a, double double_c);

 private:
  void CheckMaxNumOfBits() {
    if (a_.size() < 2) throw std::runtime_error(__FN__);

    Fr abs_c = c_.isNegative() ? -c_ : c_;
    auto num =
        a_.size() * (D + N) +
        libff::bigint<Fr::num_limbs>(abs_c.getMpz().get_mpz_t()).num_bits();
        
    assert(num < 253);
    if (num >= 253) throw std::runtime_error(__FN__);
  }
 private:
  libsnark::pb_linear_combination_array<Fr> a_;
  Fr c_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> sign_;
  libsnark::pb_variable_array<Fr> products_;
  libsnark::pb_linear_combination<Fr> product_off_;  // = product_ + 2^(D+2N)
  std::unique_ptr<libsnark::packing_gadget<Fr>> p1_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p2_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
  libsnark::pb_variable<Fr> packed_;
};

template <size_t D, size_t N>
bool Mul2Gadget<D, N>::Test(std::vector<double> const& double_a,
                            double double_c) {
  auto double_ret = double_c;
  for (auto& i : double_a) {
    double_ret = double_ret * i;
  }

  std::cout << "double_ret: " << double_ret << "\n";

  std::vector<Fr> fr_a(double_a.size());
  for (size_t i = 0; i < double_a.size();++i) {
    fr_a[i] = DoubleToRational<D, N>(double_a[i]);
  }

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable_array<Fr> pb_a;
  pb_a.allocate(pb, double_a.size(), "TestMul2 a");

  Fr c = DoubleToRational<D, N>(double_c);
  Mul2Gadget<D, N> gadget(pb, pb_a, c, "TestMul2 gadget");
  
  gadget.generate_r1cs_constraints();

  for (size_t i = 0; i < double_a.size(); ++i) {
    pb.val(pb_a[i]) = fr_a[i];
  }

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

inline bool TestMul2() {
  Tick tick(__FN__);
  std::vector<double> double_a{-7.3,21.1,-1.4};
  double double_c = -2.4;
  return Mul2Gadget<32, 32>::Test(double_a, double_c);
}

}  // namespace circuit::fixed_point