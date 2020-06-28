#pragma once

#include "./precision_gadget.h"

namespace circuit::fixed_point {

// insure a >=-kFrDN && a < kFrDN, that is [-2^(D+N), 2^(D+N)-1]
// ret = a >= 0? a : 0
// num_constraints:
// num_variables:

template <size_t D, size_t N, size_t M>
class Relu2Gadget : public libsnark::gadget<Fr> {
  static_assert(D + N < 253, "invalid D,N");
  static_assert(N > M, "invalid N or M");

 public:
  Relu2Gadget(libsnark::protoboard<Fr>& pb,
              const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    a_.allocate(this->pb, FMT(this->annotation_prefix, " a"));
    precision_gadget_.reset(new PrecisionGadget<D, N, M>(
        this->pb, a_, FMT(this->annotation_prefix, " precision_gadget")));
    ret_.allocate(pb, FMT(this->annotation_prefix, " ret"));
    generate_r1cs_constraints();
  }

  void Assign(Fr const& a) {
    this->pb.val(a_) = a;
    generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

 private:
  void generate_r1cs_constraints() {
    precision_gadget_->generate_r1cs_constraints();

    // ret = sign_*b
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(precision_gadget_->ret(),
                                      precision_gadget_->sign(), ret_),
        FMT(this->annotation_prefix, " ret = sign? b:0"));
  }

  void generate_r1cs_witness() {
    precision_gadget_->generate_r1cs_witness();

    auto a = this->pb.lc_val(a_);
    auto b = this->pb.lc_val(precision_gadget_->ret());
    auto sign = this->pb.lc_val(precision_gadget_->sign());
    this->pb.val(ret_) = sign == Fr(0) ? Fr(0) : b;
  }  

 public:
  static bool Test(double const& dx) {
    Fr x = DoubleToRational<D, N>(dx);
    std::unique_ptr<Relu2Gadget<D, N, M>> gadget;
    libsnark::protoboard<Fr> pb;
    gadget = std::make_unique<Relu2Gadget<D, N, M>>(pb, "Relu2Gadget");
    gadget->Assign(x);
    assert(pb.is_satisfied());

#ifdef _DEBUG
    double dr = RationalToDouble<D, M>(pb.val(gadget->ret()));
    if (dx < 0) {
      assert(dr == 0);
    } else {
      assert(std::abs(dx - dr) < 0.001);
    }
#endif
    return pb.is_satisfied();
  }

 private:
  libsnark::pb_variable<Fr> a_;
  std::unique_ptr<PrecisionGadget<D, N, M>> precision_gadget_;
  libsnark::pb_variable<Fr> ret_;
};

inline bool TestRelu2() {
  Tick tick(__FN__);
  constexpr size_t D = 18;
  constexpr size_t N = 32;
  constexpr size_t M = 24;
  std::vector<bool> rets;
  rets.push_back(Relu2Gadget<D, N, M>::Test(3.124));
  rets.push_back(Relu2Gadget<D, N, M>::Test(-22.212));
  rets.push_back(Relu2Gadget<D, N, M>::Test(0.00123));
  rets.push_back(Relu2Gadget<D, N, M>::Test(-0.00123));
  rets.push_back(Relu2Gadget<D, N, M>::Test(234.123));
  rets.push_back(Relu2Gadget<D, N, M>::Test(-23.11224));
  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point