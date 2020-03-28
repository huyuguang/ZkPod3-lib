#pragma once

#include "./abs_gadget.h"
#include "./div_gadget.h"
#include "./mul_gadget.h"

namespace circuit::fixed_point {

// e^x = limit(n->infinte) (1+x/n)^n
// let n = 1024, we have
// double my_exp(double x) {
//  x = 1.0 + x / 1024;
//  x *= x; x *= x; x *= x; x *= x;
//  x *= x; x *= x; x *= x; x *= x;
//  x *= x; x *= x;
//  return x;
//}

// x can be negative

template <size_t D, size_t N>
class ExpGadget : public libsnark::gadget<Fr> {
 public:
  ExpGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& x, size_t accuracy_level,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), x_(x) {
    assert(accuracy_level >= 6 && accuracy_level <= 20);
    auto constances = RationalConst<D, N>();

    mulc_gadget_.reset(
        new MulVCGadget<D, N>(this->pb, x_, 1.0 / (1ULL << accuracy_level),
                              FMT(this->annotation_prefix, " mulc_gadget")));

    base_x_.assign(this->pb, mulc_gadget_->ret() + constances.kFrN);

    mul_gadgets_.resize(accuracy_level);
    for (size_t i = 0; i < accuracy_level; ++i) {
      libsnark::pb_linear_combination<Fr> last =
          i == 0 ? base_x_ : mul_gadgets_[i - 1]->ret();
      mul_gadgets_[i].reset(new MulVVGadget<D, N>(
          pb, last, last, FMT(this->annotation_prefix, " mul_gadgets_%zu", i)));
    }
  }

  void generate_r1cs_constraints() {
    mulc_gadget_->generate_r1cs_constraints();
    for (auto& i : mul_gadgets_) {
      i->generate_r1cs_constraints();
    }
  }

  void generate_r1cs_witness() {
    mulc_gadget_->generate_r1cs_witness();

    base_x_.evaluate(this->pb);

    for (auto& i : mul_gadgets_) {
      i->generate_r1cs_witness();
    }
  }

  libsnark::pb_linear_combination<Fr> ret() const {
    return mul_gadgets_.back()->ret();
  }

  static bool Test(double double_x);

 private:
  libsnark::pb_linear_combination<Fr> x_;
  std::unique_ptr<MulVCGadget<D, N>> mulc_gadget_;
  libsnark::pb_linear_combination<Fr> base_x_;
  std::vector<std::unique_ptr<MulVVGadget<D, N>>> mul_gadgets_;
};

template <size_t D, size_t N>
bool ExpGadget<D, N>::Test(double double_x) {
  std::cout << "e^" << double_x << "\n";
  auto double_ret = std::exp(double_x);
  std::cout << "double_ret: " << double_ret << "\n";

  Fr x = DoubleToRational<D, N>(double_x);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_x;
  pb_x.allocate(pb, "ExpGadget::Test a");
  size_t kAccuracyLevel = 10;
  ExpGadget<D, N> gadget(pb, pb_x, kAccuracyLevel, "ExpGadget::Test");
  gadget.generate_r1cs_constraints();
  pb.val(pb_x) = x;
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t"
            << fp::RationalToDouble<D, N>(fr_ret) << "\n";

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestExp() {
  Tick tick(__FN__);
  bool ret;
  std::vector<bool> rets;
  double double_x;

  double_x = -9.8;
  ret = ExpGadget<32, 32>::Test(double_x);
  rets.push_back(ret);

  double_x = -1.0;
  ret = ExpGadget<32, 32>::Test(double_x);
  rets.push_back(ret);

  double_x = 9.5;
  ret = ExpGadget<32, 32>::Test(double_x);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point