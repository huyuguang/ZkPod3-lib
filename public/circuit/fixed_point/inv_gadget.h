#pragma once

#include "./div_gadget.h"

namespace circuit::fixed_point {

// NOTE: b must >=0
template <size_t D, size_t N>
class InvGadget : public DivGadget<D, N> {
 public:
  InvGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& b, const std::string&
                annotation_prefix = "")
      : DivGadget<D, N>(pb, GetA(pb), b, annotation_prefix) {}

 private:
  libsnark::pb_linear_combination<Fr> GetA(libsnark::protoboard<Fr>& pb) {
    libsnark::pb_linear_combination<Fr> a;
    a.assign(pb, libsnark::pb_variable<Fr>(0) * RationalConst<D, N>().kFrN);
    return a;
  }
};

inline bool TestInv() {
  Tick tick(__FN__);
  constexpr size_t D = 10;
  constexpr size_t N = 2;
  Fr b = -7;

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_b;
  pb_b.allocate(pb, "TestInv b");
  InvGadget<D, N> gadget(pb, pb_b, "TestInv");
  gadget.generate_r1cs_constraints();
  pb.val(pb_b) = b;
  gadget.generate_r1cs_witness();
  assert(pb.is_satisfied());

  std::cout << pb.val(gadget.ret()) << "\n";
  std::cout << -pb.val(gadget.ret()) << "\n";
  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";
  return true;
}
}  // namespace circuit::fixed_point