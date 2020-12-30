#pragma once

#include "./sign_gadget.h"

namespace circuit::fixed_point {

// ret = abs(a), a >= 0? a : -a
template <size_t D, size_t N>
class AbsGadget : public libsnark::gadget<Fr> {
  static_assert(D + N < 253, "invalid D,N");

 public:
  AbsGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& a,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a) {
    ret_.allocate(pb, FMT(this->annotation_prefix, " ret"));
    sign_gadget_.reset(new SignGadget<D, N>(
        this->pb, a_, FMT(this->annotation_prefix, " sign")));
  }

  void generate_r1cs_constraints() {
    sign_gadget_->generate_r1cs_constraints();
    // sign * 2a = a + ret
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(sign_gadget_->ret(), a_ * 2, ret_ + a_),
        FMT(this->annotation_prefix, " sign? a : -a"));
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    Fr fr_a = this->pb.lc_val(a_);
    this->pb.val(ret_) = fr_a.isNegative() ? -fr_a : fr_a;
    sign_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

 private:
  libsnark::pb_linear_combination<Fr> a_;
  std::unique_ptr<SignGadget<D, N>> sign_gadget_;
  libsnark::pb_variable<Fr> ret_;
};

inline bool TestAbs() {
  Tick tick(__FN__);
  constexpr size_t D = 8;
  constexpr size_t N = 2;
  Fr a = -100;
  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  pb_a.allocate(pb, " a");
  AbsGadget<D, N> gadget(pb, pb_a, "AbsGadget");
  gadget.generate_r1cs_constraints();
  pb.val(pb_a) = a;
  gadget.generate_r1cs_witness();
  CHECK(pb.is_satisfied(), "");

  std::cout << Tick::GetIndentString() << pb.val(gadget.ret()) << "\n";
  std::cout << Tick::GetIndentString()
            << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << Tick::GetIndentString()
            << "num_variables: " << pb.num_variables() << "\n";
  return true;
}
}  // namespace circuit::fixed_point