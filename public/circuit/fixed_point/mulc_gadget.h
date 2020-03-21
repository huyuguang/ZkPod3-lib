#pragma once

#include "./mul_gadget.h"

namespace circuit::fixed_point {

template <size_t D, size_t N>
class MulCGadget : public MulGadget<D, N> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  MulCGadget(libsnark::protoboard<Fr>& pb,
             libsnark::pb_linear_combination<Fr> const& a, double const& b,
             const std::string& annotation_prefix = "")
      : MulGadget<D, N>(pb, a, lc_b(pb, b), annotation_prefix) {}

 private:
  libsnark::pb_linear_combination<Fr> lc_b(libsnark::protoboard<Fr>& pb,
                                           double const& b) {
    libsnark::pb_linear_combination<Fr> lc_b;
    Fr fr_b = DoubleToRational<D, N>(b);
    lc_b.assign(pb, fr_b);
    return lc_b;
  }

 private:
  std::unique_ptr<MulGadget<D, N>> mul_gadget_;
};
}  // namespace circuit::fixed_point
