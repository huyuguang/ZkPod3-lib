#pragma once

#include "./abs_gadget.h"
#include "./inv_gadget.h"
#include "./mulc_gadget.h"

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
class Exp2Gadget : public libsnark::gadget<Fr> {
 public:
  Exp2Gadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& x, size_t level,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), x_(x) {
    assert(level >= 6 && level <= 20);
    auto constances = RationalConst<D, N>();

    mulc_gadget_.reset(
        new MulCGadget<D, N>(this->pb, x_, 1.0 / (1ULL << level),
                             FMT(this->annotation_prefix, " mulc_gadget")));

    base_x_.assign(this->pb, mulc_gadget_->ret() + constances.kFrN);

    mul_gadgets_.resize(level);
    for (size_t i = 0; i < level; ++i) {
      libsnark::pb_linear_combination<Fr> last =
          i == 0 ? base_x_ : mul_gadgets_[i - 1]->ret();
      mul_gadgets_[i].reset(new MulGadget<D, N>(
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
  std::unique_ptr<MulCGadget<D, N>> mulc_gadget_;
  libsnark::pb_linear_combination<Fr> base_x_;
  std::vector<std::unique_ptr<MulGadget<D, N>>> mul_gadgets_;
};

template <size_t D, size_t N>
bool Exp2Gadget<D, N>::Test(double double_x) {
  auto double_ret = std::exp(double_x);
  std::cout << "double_ret: " << double_ret << "\n";

  Fr x = DoubleToRational<D, N>(double_x);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_x;
  pb_x.allocate(pb, "Exp2Gadget::Test a");
  size_t kLevel = 10;
  Exp2Gadget<D, N> gadget(pb, pb_x, kLevel, "Exp2Gadget::Test");
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

inline bool TestExp2() {
  Tick tick(__FN__);
  double double_x = -10;
  std::cout << "e^" << double_x << "\n";
  Exp2Gadget<32, 32>::Test(double_x);

  double_x = -1.0;
  std::cout << "e^" << double_x << "\n";
  Exp2Gadget<32, 32>::Test(double_x);

  double_x = 10;
  std::cout << "e^" << double_x << "\n";
  Exp2Gadget<32, 32>::Test(double_x);

  return true;
}
}  // namespace circuit::fixed_point