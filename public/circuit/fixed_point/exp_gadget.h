#pragma once

#include "./abs_gadget.h"
#include "./mul2_gadget.h"

namespace circuit::fixed_point {

// X = x * N (N=2^n)
// X1 = X
// X2 = x^2 * N/2! = X1 * X * (N/2) / N^2
// X3 = x^3 * N/3! = X2 * X * (N/3) / N^2
// X4 = x^4 * N/4! = X3 * X * (N/4) / N^2
// f(n) = f(n-1) * X * (N/n) / N^2

// exp(x) * N = N + X1 + X2 + ...
// x can be negative

template <size_t D, size_t N>
class ExpGadget : public libsnark::gadget<Fr> {

 public:
  ExpGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& x,
            size_t series_num,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        x_(x) {
    auto constances = RationalConst<D, N>();
    libsnark::linear_combination<Fr> ret(x_ + constances.kFrN);

    mul2_gadgets_.resize(series_num);
    for (size_t i = 0; i < series_num; ++i) {
      auto mpz_coef = constances.kMpzN / (i + 2);
      Fr fr_coef;
      fr_coef.setMpz(mpz_coef);
      libsnark::pb_linear_combination<Fr> last =
          i == 0 ? x_ : mul2_gadgets_[i - 1]->ret();
      mul2_gadgets_[i].reset(new Mul2Gadget<D, N>(
          pb, x, last, fr_coef,
          FMT(this->annotation_prefix, " Mul2Gadget_%zu", i)));
      ret = ret + mul2_gadgets_[i]->ret();
    }
    ret_.assign(this->pb, ret);
  }

  void generate_r1cs_constraints() {
    for (auto& i : mul2_gadgets_) {
      i->generate_r1cs_constraints();
    }
  }

  void generate_r1cs_witness() {
    auto constances = RationalConst<D, N>();
    x_.evaluate(this->pb);
    
    for (auto& i : mul2_gadgets_) {
      i->generate_r1cs_witness();
    }

    ret_.evaluate(this->pb);

#ifdef _DEBUG
    Fr x = this->pb.lc_val(x_);
    Fr ret = this->pb.lc_val(ret_);
    std::cout << "exp(" << fp::RationalToDouble<N>(x)
              << ") = " << fp::RationalToDouble<N>(ret) << "\n";
#endif
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

 private:
  libsnark::pb_linear_combination<Fr> x_;
  std::vector<std::unique_ptr<Mul2Gadget<D,N>>> mul2_gadgets_;
  libsnark::pb_linear_combination<Fr> ret_;
};

inline bool TestExp() {
  constexpr size_t D = 32;
  constexpr size_t N = 24;
  double double_a = -0.5;  
  std::cout << "e^" << double_a << "\n";

  Fr a = DoubleToRational<N>(double_a);
  // (uint64_t)(double_a * (1ULL << N));

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_x;
  pb_x.allocate(pb, "TestExp x");
  ExpGadget<D, N> gadget(pb, pb_x, 9, "ExpGadget");
  gadget.generate_r1cs_constraints();
  pb.val(pb_x) = a;
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());

  Fr ret = pb.lc_val(gadget.ret());
  std::cout << ret << "\n";
  std::cout << RationalToDouble<N>(ret) << "\n";

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}
}  // namespace circuit::fixed_point