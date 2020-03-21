#pragma once

#include "./abs_gadget.h"
#include "./div_gadget.h"
#include "./mul_gadget.h"

namespace circuit::fixed_point {

// e^x = limit(n->infinte) 1+x+x^2/2!+x^3/3!+...x^n/n!
// X = x * N (N=2^n)
// X1 = X
// X2 = x^2 * N/2! = X1 * X * (N/2) / N^2
// X3 = x^3 * N/3! = X2 * X * (N/3) / N^2
// X4 = x^4 * N/4! = X3 * X * (N/4) / N^2
// f(n) = f(n-1) * X * (N/n) / N^2

// exp(x) * N = N + X1 + X2 + ...
// x can be negative
// We do not need to check sign_x, but if not, the e^negative wants larger
// series_num

template <size_t D, size_t N>
class ExpGadget : public libsnark::gadget<Fr> {
 public:
  ExpGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& x,
            libsnark::pb_linear_combination<Fr> const& sign_x,
            size_t series_num, const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), x_(x), sign_x_(sign_x) {
    auto constances = RationalConst<D, N>();
    abs_x_.allocate(pb, FMT(this->annotation_prefix, " abs_x"));
    ret_.allocate(pb, FMT(this->annotation_prefix, " ret"));

    libsnark::linear_combination<Fr> e_abs_x(abs_x_ + constances.kFrN);
    mul_gadgets_.resize(series_num);
    for (size_t i = 0; i < series_num; ++i) {
      libsnark::pb_linear_combination<Fr> last =
          i == 0 ? abs_x_ : mul_gadgets_[i - 1]->ret();
      mul_gadgets_[i].reset(new MulVVCGadget<D, N>(
          pb, abs_x_, last, 1.0/(i+2),
          FMT(this->annotation_prefix, " MulGadget_%zu", i)));
      e_abs_x = e_abs_x + mul_gadgets_[i]->ret();
    }
    e_abs_x_.assign(this->pb, e_abs_x);

    libsnark::pb_linear_combination<Fr> pb_lc_sign;
    pb_lc_sign.assign(pb, FrOne()); // e_abs_x_ always > 0
    inv_gadget_.reset(new InvGadget<D, N>(
        this->pb, e_abs_x_, &pb_lc_sign, FMT(this->annotation_prefix, " inv_gadget")));
  }

  void generate_r1cs_constraints() {
    // abs_x = sign_x? x:-x
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_x_, x_ * 2, x_ + abs_x_),
        FMT(this->annotation_prefix, " abs_x = sign_x? x:-x"));

    for (auto& i : mul_gadgets_) {
      i->generate_r1cs_constraints();
    }

    inv_gadget_->generate_r1cs_constraints();

    // ret = sign_x? e_abs_x : 1/e_abs_x
    auto e_abs_x_inv = inv_gadget_->ret();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_x_, e_abs_x_ - e_abs_x_inv,
                                      ret_ - e_abs_x_inv),
        FMT(this->annotation_prefix, " ret = sign_x? e_abs_x : 1/e_abs_x"));
  }

  void generate_r1cs_witness() {
    auto constances = RationalConst<D, N>();
    x_.evaluate(this->pb);
    sign_x_.evaluate(this->pb);

    Fr x = this->pb.lc_val(x_);
    Fr sign_x = this->pb.lc_val(sign_x_);
    Fr check_sign_x = x.isNegative() ? 0 : 1;
    assert(sign_x == check_sign_x);
    this->pb.val(abs_x_) = x.isNegative() ? -x : x;

    for (auto& i : mul_gadgets_) {
      i->generate_r1cs_witness();
    }

    e_abs_x_.evaluate(this->pb);
    Fr e_abs_x = this->pb.lc_val(e_abs_x_);

    inv_gadget_->generate_r1cs_witness();
    Fr e_abs_x_inv = this->pb.lc_val(inv_gadget_->ret());

    this->pb.val(ret_) = x.isNegative() ? e_abs_x_inv : e_abs_x;
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

  static bool Test(double double_x);

 private:
  libsnark::pb_linear_combination<Fr> x_;
  libsnark::pb_linear_combination<Fr> sign_x_;
  libsnark::pb_variable<Fr> abs_x_;
  std::vector<std::unique_ptr<MulGadget<D, N>>> mul_gadgets_;
  std::unique_ptr<InvGadget<D, N>> inv_gadget_;
  libsnark::pb_linear_combination<Fr> e_abs_x_;
  libsnark::pb_variable<Fr> ret_;
};

template <size_t D, size_t N>
bool ExpGadget<D, N>::Test(double double_x) {
  std::cout << "exp(" << double_x << ")\n";
  auto double_ret = std::exp(double_x);
  std::cout << "should be: " << double_ret << "\n";

  Fr x = DoubleToRational<D, N>(double_x);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_x;
  libsnark::pb_variable<Fr> pb_sign_x;
  pb_x.allocate(pb, "ExpGadget::Test a");
  pb_sign_x.allocate(pb, "ExpGadget::Test sign_x");
  size_t kSeriesNum = 9;
  ExpGadget<D, N> gadget(pb, pb_x, pb_sign_x, kSeriesNum, "ExpGadget::Test");
  gadget.generate_r1cs_constraints();
  pb.val(pb_x) = x;
  pb.val(pb_sign_x) = Fr(double_x >= 0 ? 1 : 0);
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t"
            << fp::RationalToDouble<D, N>(fr_ret) << "\n";

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestExp() {
  Tick tick(__FN__);
  double double_x = -7.3;
  return ExpGadget<32, 32>::Test(double_x);
}
}  // namespace circuit::fixed_point