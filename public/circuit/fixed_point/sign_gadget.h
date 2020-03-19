#pragma once

#include "./misc.h"

namespace circuit::fixed_point {

// kFrW = 2^W
// insure a >=-kFrW && a < kFrW, that is [-2^W, 2^W-1]
// ret = a >= 0? 1 : 0
// num_constraints:
// num_variables:

template <size_t W>
class SignGadget : public libsnark::gadget<Fr> {
  static_assert(W < 253, "invalid W");

 public:
  SignGadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& a,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a) {
    a_off_.assign(this->pb, a_ + BigIntConst<W>().kFrW);
    // 2kFrW > a_ + kFrW >=0, so use W+1 bits
    bits_.allocate(this->pb, W + 1, FMT(this->annotation_prefix, " bits"));
    packing_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb, bits_, a_off_, FMT(this->annotation_prefix, " packing")));
  }

  void generate_r1cs_constraints() {
    packing_gadget_->generate_r1cs_constraints(true);
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    a_off_.evaluate(this->pb);

    // std::cout << "a: " << this->pb.lc_val(a_) << "\n";
    // std::cout << "a_off: " << this->pb.lc_val(a_off_) << "\n";

    packing_gadget_->generate_r1cs_witness_from_packed();

    if (this->pb.lc_val(a_).isNegative()) {
      assert(this->pb.val(ret()) == 0);
    } else {
      assert(this->pb.val(ret()) == 1);
    }
  }

  libsnark::pb_variable<Fr> ret() const { return bits_[W]; }

  static bool Test(Fr const& x, bool except) {
    libsnark::protoboard<Fr> pb;
    libsnark::pb_variable<Fr> pb_x;
    pb_x.allocate(pb, "Test x");
    SignGadget<W> gadget(pb, pb_x, "SignGadget");
    gadget.generate_r1cs_constraints();
    pb.val(pb_x) = x;
    gadget.generate_r1cs_witness();
    assert(pb.is_satisfied() == except);
    return pb.is_satisfied() == except;
  }

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> a_off_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> packing_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
};

inline bool TestSign() {
  constexpr size_t W = 10;
  std::vector<bool> rets;
  rets.push_back(SignGadget<W>::Test(Fr(100), true));
  rets.push_back(SignGadget<W>::Test(Fr(0), true));
  rets.push_back(SignGadget<W>::Test(Fr(-1), true));
  rets.push_back(SignGadget<W>::Test(Fr(-100), true));
  rets.push_back(SignGadget<W>::Test(Fr(-1024), true));

#ifndef _DEBUG  // would trigger assert in basic_gadgets
  rets.push_back(SignGadget<W>::Test(Fr(1024), false));
  rets.push_back(SignGadget<W>::Test(Fr(-1025), false));
  rets.push_back(SignGadget<W>::Test(Fr(1025), false));
  rets.push_back(SignGadget<W>::Test(Fr(122222244), false));
#endif

  // std::cout << pb.val(gadget.ret()) << "\n";
  // std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  // std::cout << "num_variables: " << pb.num_variables() << "\n";
  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::fixed_point