#pragma once

#include "../fixed_point/fixed_point.h"

// a simple wrapper for libsnark::inner_product_gadget

namespace circuit::vgg16 {

class IpGadget : public libsnark::gadget<Fr> {
 public:
  IpGadget(libsnark::protoboard<Fr>& pb,
           const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    data_.allocate(this->pb, 9, FMT(this->annotation_prefix, " data"));
    para_.allocate(this->pb, 9, FMT(this->annotation_prefix, " para"));
    ret_.allocate(this->pb, " ret");

    ip_.reset(new libsnark::inner_product_gadget<Fr>(this->pb, data_, para_,
                                                     ret_, "ip"));

    ip_->generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> ret() const { return ret_; }

  void Assign(std::array<Fr const*, 9> const& data,
              std::array<Fr const*, 9> const& para) {
    std::vector<Fr> vec_data(9);
    std::vector<Fr> vec_para(9);
    for (size_t i = 0; i < 9; ++i) {
      vec_data[i] = *data[i];
      vec_para[i] = *para[i];
    }
    data_.fill_with_field_elements(this->pb, vec_data);
    para_.fill_with_field_elements(this->pb, vec_para);
    ip_->generate_r1cs_witness();
  }

 private:
  libsnark::pb_variable_array<Fr> data_;
  libsnark::pb_variable_array<Fr> para_;
  libsnark::pb_variable<Fr> ret_;
  std::unique_ptr<libsnark::inner_product_gadget<Fr>> ip_;
};

inline bool TestIp() {
  // Tick tick(__FN__);
  // std::vector<bool> rets;
  //{
  //  constexpr size_t DataCol = 4;
  //  constexpr size_t DataRow = 4;
  //  constexpr size_t ConvCol = 3;
  //  constexpr size_t ConvRow = 3;
  //  std::array<double, DataCol* DataRow> double_data = {
  //      0.666667, 0.203922, 0, 0, 0.996078, 0.54902,   0, 0,
  //      0.996078, 0.415686, 0, 0, 0.819608, 0.0705882, 0, 0};
  //  std::array<double, ConvCol* ConvRow + 1> double_para = {
  //      -1.24212,  -1.46405, -1.25206, 0.767129, 0.201014,
  //      -0.873962, 0.776801, 0.977157, 0.266173, 0.187466};
  //  bool ret = ConvGadget<4, 24, DataCol, DataRow, ConvCol, ConvRow>::Test(
  //      double_data, double_para);
  //  rets.push_back(ret);
  //}
  // return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
  return false;
}
};  // namespace circuit::vgg16