#pragma once

#include "./permutation_gadget.h"

namespace circuit {
// prove input is permutaion of output
// input always 1~D
class SudokuGadget {
 private:
  std::vector<libsnark::pb_variable<Fr>> input_vars_;
  std::vector<libsnark::pb_variable<Fr>> output_vars_;
  std::unique_ptr<PermutationGadget> sub_gadget_;

 public:
  SudokuGadget(libsnark::protoboard<Fr>& pb, size_t D) {
    input_vars_.resize(D);
    output_vars_.resize(D);
    for (size_t packet_idx = 0; packet_idx < D; ++packet_idx) {
      input_vars_[packet_idx].allocate(pb, FMT("", "input_%zu", packet_idx));
    }
    for (size_t packet_idx = 0; packet_idx < D; ++packet_idx) {
      output_vars_[packet_idx].allocate(pb, FMT("", "output_%zu", packet_idx));
    }

    sub_gadget_.reset(
        new PermutationGadget(pb, input_vars_, output_vars_, "sudoku"));
    sub_gadget_->generate_r1cs_constraints();
  }

  void Assign(std::vector<Fr> const& output) {
    std::vector<Fr> input(output.size());
    for (size_t i = 0; i < input.size(); ++i) {
      input[i] = i + 1;
    }
    sub_gadget_->generate_r1cs_witness(input, output);
  }
};
}  // namespace circuit