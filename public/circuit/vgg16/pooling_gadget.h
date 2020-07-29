#pragma once

#include "../fixed_point/fixed_point.h"

// a simple wrapper for fixed_point::MaxGadget

namespace circuit::vgg16 {

template <size_t D, size_t N, size_t S = 4>
class PoolingGadget : public libsnark::gadget<Fr> {
 public:
  PoolingGadget(libsnark::protoboard<Fr>& pb,
                const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    data_.allocate(this->pb, S, FMT(this->annotation_prefix, " data"));

    max_.reset(new fixed_point::MaxGadget<D, N>(this->pb, data_, "max"));

    max_->generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> ret() const { return max_->ret(); }

  void Assign(std::array<Fr const*, S> const& data) {
    std::vector<Fr> vec_data(S);
    for (size_t i = 0; i < S; ++i) {
      vec_data[i] = *data[i];
    }
    data_.fill_with_field_elements(this->pb, vec_data);
    max_->generate_r1cs_witness();
  }

 private:
  libsnark::pb_variable_array<Fr> data_;
  std::unique_ptr<fixed_point::MaxGadget<D, N>> max_;
};

inline bool TestPooling() { return true; }
};  // namespace circuit::vgg16