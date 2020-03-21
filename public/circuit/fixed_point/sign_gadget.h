#pragma once

#include "./misc.h"

namespace circuit::fixed_point {

// insure a >=-kFrDN && a < kFrDN, that is [-2^(D+N), 2^(D+N)-1]
// ret = a >= 0? 1 : 0
// num_constraints:
// num_variables:

template <size_t D, size_t N>
class SignGadget : public libsnark::gadget<Fr> {
  static_assert(D + N < 253, "invalid D,N");

 public:
  SignGadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& a,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a) {
    a_off_.assign(this->pb, a_ + RationalConst<D, N>().kFrDN);
    // 2kFrDN > a_ + kFrDN >=0, so use D+N+1 bits
    bits_.allocate(this->pb, D + N + 1, FMT(this->annotation_prefix, " bits"));
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

  libsnark::pb_variable<Fr> ret() const { return bits_[D + N]; }

  static bool Test(Fr const& x, bool except) {
    libsnark::protoboard<Fr> pb;
    libsnark::pb_variable<Fr> pb_x;
    pb_x.allocate(pb, "Test x");
    SignGadget<D, N> gadget(pb, pb_x, "SignGadget");
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
  Tick tick(__FN__);
  constexpr size_t D = 8;
  constexpr size_t N = 2;
  std::vector<bool> rets;
  rets.push_back(SignGadget<D, N>::Test(Fr(100), true));
  rets.push_back(SignGadget<D, N>::Test(Fr(0), true));
  rets.push_back(SignGadget<D, N>::Test(Fr(-1), true));
  rets.push_back(SignGadget<D, N>::Test(Fr(-100), true));
  rets.push_back(SignGadget<D, N>::Test(Fr(-1024), true));

#ifndef _DEBUG  // would trigger assert in basic_gadgets
  rets.push_back(SignGadget<D, N>::Test(Fr(1024), false));
  rets.push_back(SignGadget<D, N>::Test(Fr(-1025), false));
  rets.push_back(SignGadget<D, N>::Test(Fr(1025), false));
  rets.push_back(SignGadget<D, N>::Test(Fr(122222244), false));
#endif

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}

template <size_t D, size_t N>
using TypeGadget = SignGadget<D, N>;
}  // namespace circuit::fixed_point