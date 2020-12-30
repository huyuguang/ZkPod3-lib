#pragma once

#include "./misc.h"

namespace circuit::fixed_point {

// precision reduction
// insure a >=-kFrDN && a < kFrDN, that is [-2^(D+N), 2^(D+N)-1]
// convert to <D,M> which M < N
// sign = a >= 0? 1 : 0
// num_constraints:
// num_variables:

template <size_t D, size_t N, size_t M>
class PrecisionGadget : public libsnark::gadget<Fr> {
  static_assert(D + N < 253, "invalid D or N");
  static_assert(N > M, "invalid N or M");

 public:
  PrecisionGadget(libsnark::protoboard<Fr>& pb,
                  libsnark::pb_linear_combination<Fr> const& a,
                  const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a) {
    a_off_.assign(this->pb, a_ + RationalConst<D, N>().kFrDN);
    // 2kFrDN > a_ + kFrDN >=0, so use D+N+1 bits
    bits_.allocate(this->pb, D + N + 1, FMT(this->annotation_prefix, " bits"));
    p1_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb, bits_, a_off_, FMT(this->annotation_prefix, " p1")));

    packed_.allocate(pb, FMT(this->annotation_prefix, " packed"));
    p2_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb,
        libsnark::pb_variable_array<Fr>(bits_.begin() + N - M, bits_.end()),
        packed_, FMT(this->annotation_prefix, " p2")));

    ret_.assign(this->pb, libsnark::linear_combination<Fr>(packed_) -
                              RationalConst<D, M>().kFrDN);
    sign_.assign(this->pb, p2_gadget_->bits[D + M]);
  }

  void generate_r1cs_constraints() {
    p1_gadget_->generate_r1cs_constraints(true);
    p2_gadget_->generate_r1cs_constraints(false);
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    a_off_.evaluate(this->pb);

    // std::cout << "a: " << this->pb.lc_val(a_) << "\n";
    // std::cout << "a_off: " << this->pb.lc_val(a_off_) << "\n";

    p1_gadget_->generate_r1cs_witness_from_packed();
    p2_gadget_->generate_r1cs_witness_from_bits();
    ret_.evaluate(this->pb);
    sign_.evaluate(this->pb);

#ifdef _DEBUG
    Fr a = this->pb.lc_val(a_);
    Fr b = this->pb.lc_val(ret_);
    Fr sign = this->pb.lc_val(sign_);
    assert(a.isNegative() == b.isNegative());
    assert(Fr(a.isNegative() ? 0 : 1) == sign);
    double da = RationalToDouble<D, N>(a);
    double db = RationalToDouble<D, M>(b);
    CHECK(std::abs(da - db) < 0.001, "");
    auto val = ReducePrecision<D, N, M>(this->pb.lc_val(a_));
    CHECK(val == this->pb.lc_val(ret_), "");
#endif
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

  libsnark::pb_linear_combination<Fr> sign() const { return sign_; }

  static bool Test(double const& dx) {
    Fr x = DoubleToRational<D, N>(dx);
    libsnark::protoboard<Fr> pb;
    libsnark::pb_variable<Fr> pb_x;
    pb_x.allocate(pb, "Test x");
    PrecisionGadget<D, N, M> gadget(pb, pb_x, "PrecisionGadget");
    gadget.generate_r1cs_constraints();
    pb.val(pb_x) = x;
    gadget.generate_r1cs_witness();
    assert(pb.is_satisfied());
    return pb.is_satisfied();
  }

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> a_off_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p1_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p2_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
  libsnark::pb_variable<Fr> packed_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> sign_;
};

inline bool TestPrecision() {
  Tick tick(__FN__);
  constexpr size_t D = 12;
  constexpr size_t N = 32;
  constexpr size_t M = 24;
  std::vector<bool> rets;
  rets.push_back(PrecisionGadget<D, N, M>::Test(3.124));
  rets.push_back(PrecisionGadget<D, N, M>::Test(-22.212));
  rets.push_back(PrecisionGadget<D, N, M>::Test(0.00123));
  rets.push_back(PrecisionGadget<D, N, M>::Test(-0.00123));
  rets.push_back(PrecisionGadget<D, N, M>::Test(234.123));
  rets.push_back(PrecisionGadget<D, N, M>::Test(-23.11224));

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point