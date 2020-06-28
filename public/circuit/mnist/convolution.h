#pragma once

#include "../fixed_point/fixed_point.h"

// combine conv, relu, max_pooling, flatten

namespace circuit::cnn {

template <size_t D, size_t N, size_t DataCol, size_t DataRow, size_t ConvCol,
          size_t ConvRow>
class ConvGadget : public libsnark::gadget<Fr> {
 public:
  ConvGadget(libsnark::protoboard<Fr>& pb,
             const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    data_.allocate(this->pb, DataCol * DataRow,
                   FMT(this->annotation_prefix, " data"));
    para_.allocate(this->pb, ConvCol * ConvRow + 1,
                   FMT(this->annotation_prefix, " para"));

    std::array<libsnark::pb_linear_combination_array<Fr>, kIpCount> conv_datas;
    for (size_t i = 0; i < DataRow + 1 - ConvRow; ++i) {
      for (size_t j = 0; j < DataCol + 1 - ConvCol; ++j) {
        auto const* topleft_data = &data_[i * DataCol + j];
        auto& conv_data = conv_datas[i * (DataCol + 1 - ConvCol) + j];
        conv_data.resize(ConvCol * ConvRow);
        for (size_t ii = 0; ii < ConvCol; ++ii) {
          for (size_t jj = 0; jj < ConvRow; ++jj) {
            conv_data[ii * ConvCol + jj] = topleft_data[ii * DataCol + jj];
          }
        }
      }
    }

    relu_rets_.resize(kIpCount);
    for (size_t i = 0; i < kIpCount; ++i) {
      auto& ip_gadget = ip_gadgets_[i];
      ip_gadget.reset(new fp::IpGadget<D, N>(
          this->pb, conv_datas[i], para_,
          FMT(this->annotation_prefix, " ip_gadgets_%zu", i)));
      auto& relu_gadget = relu_gadgets_[i];
      auto ip_ret = ip_gadget->ret();
      auto ip_sign = ip_gadget->sign();
      relu_gadget.reset(new fp::ReluGadget<D, N>(
          this->pb, ip_ret, &ip_sign,
          FMT(this->annotation_prefix, " relu_gadgets_%zu", i)));
      relu_rets_[i] = relu_gadget->ret();
    }

    max_gadget_.reset(new fp::MaxGadget<D, N>(
        this->pb, relu_rets_, FMT(this->annotation_prefix, " max_gadgets")));

    generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> ret() const { return max_gadget_->ret(); }

  void Assign(std::array<Fr, DataCol * DataRow> const& data,
              std::array<Fr, ConvCol * ConvRow + 1> const& para) {
    std::vector<Fr> vec_data(data.begin(), data.end());
    std::vector<Fr> vec_para(para.begin(), para.end());
    data_.fill_with_field_elements(this->pb, vec_data);
    para_.fill_with_field_elements(this->pb, vec_para);
    generate_r1cs_witness();

//#ifdef _DEBUG
//    for (size_t i = 0; i < DataCol * DataRow; ++i) {
//      double d = fp::RationalToDouble<4, 24>(data[i]);
//      std::cout << std::right << std::setw(12) << std::setfill(' ') << d;
//      if ((i + 1) % DataCol == 0) std::cout << "\n";
//    }
//    std::cout << "\n";
//    for (size_t i = 0; i < ConvCol * ConvRow + 1; ++i) {
//      double d = fp::RationalToDouble<4, 24>(para[i]);
//      std::cout << std::right << std::setw(12) << std::setfill(' ') << d;
//      if ((i + 1) % ConvCol == 0) std::cout << "\n";
//    }
//    std::cout << "\n";
//
//    double max = fp::RationalToDouble<4, 24>(this->pb.val(ret()));
//    std::cout << max << "\n";
//    std::cout << "------------------------\n";
//#endif
  }

  static bool Test(std::array<double, DataCol * DataRow> const& data,
                   std::array<double, ConvCol * ConvRow + 1> const& para);

 private:
  void generate_r1cs_constraints() {
    for (size_t i = 0; i < kIpCount; ++i) {
      ip_gadgets_[i]->generate_r1cs_constraints();
      relu_gadgets_[i]->generate_r1cs_constraints();
    }
    max_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    for (size_t i = 0; i < kIpCount; ++i) {
      ip_gadgets_[i]->generate_r1cs_witness();
      relu_gadgets_[i]->generate_r1cs_witness();
    }
    max_gadget_->generate_r1cs_witness();
  }

 private:
  static constexpr size_t kIpCount =
      (DataCol + 1 - ConvCol) * (DataRow + 1 - ConvRow);
  libsnark::pb_variable_array<Fr> data_;
  libsnark::pb_variable_array<Fr> para_;
  std::array<std::unique_ptr<fp::IpGadget<D, N>>, kIpCount> ip_gadgets_;
  std::array<std::unique_ptr<fp::ReluGadget<D, N>>, kIpCount> relu_gadgets_;
  libsnark::pb_linear_combination_array<Fr> relu_rets_;
  std::unique_ptr<fp::MaxGadget<D, N>> max_gadget_;
};

template <size_t D, size_t N, size_t DataCol, size_t DataRow, size_t ConvCol,
          size_t ConvRow>
bool ConvGadget<D, N, DataCol, DataRow, ConvCol, ConvRow>::Test(
    std::array<double, DataCol * DataRow> const& double_data,
    std::array<double, ConvCol * ConvRow + 1> const& double_para) {
  std::array<std::array<double, ConvCol * ConvRow>, kIpCount> conv_datas;
  for (size_t i = 0; i < DataRow + 1 - ConvRow; ++i) {
    for (size_t j = 0; j < DataCol + 1 - ConvCol; ++j) {
      auto const* topleft_data = &double_data[i * DataCol + j];
      auto& conv_data = conv_datas[i * (DataCol + 1 - ConvCol) + j];
      for (size_t ii = 0; ii < ConvCol; ++ii) {
        for (size_t jj = 0; jj < ConvRow; ++jj) {
          conv_data[ii * ConvCol + jj] = topleft_data[ii * DataCol + jj];
        }
      }
    }
  }

  std::array<double, kIpCount> double_relu;
  for (size_t i = 0; i < kIpCount; ++i) {
    auto const& conv_data = conv_datas[i];
    double ip = std::inner_product(conv_data.begin(), conv_data.end(),
                                   double_para.begin(), double_para.back());
    ip = ip >= 0 ? ip : 0;
    double_relu[i] = ip;
  }

  double double_ret = *std::max_element(double_relu.begin(), double_relu.end());

  std::array<Fr, DataCol * DataRow> data;
  for (size_t i = 0; i < double_data.size(); ++i) {
    data[i] = fp::DoubleToRational<D, N>(double_data[i]);
  }
  std::array<Fr, ConvCol * ConvRow + 1> para;
  for (size_t i = 0; i < double_para.size(); ++i) {
    para[i] = fp::DoubleToRational<D, N>(double_para[i]);
  }

  libsnark::protoboard<Fr> pb;
  ConvGadget<D, N, DataCol, DataRow, ConvCol, ConvRow> gadget(pb, "TestConv");
  gadget.Assign(data, para);
  assert(pb.is_satisfied());
  if (!pb.is_satisfied()) return false;

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << "fr_ret: " << fr_ret << "\t"
            << fp::RationalToDouble<D, N>(fr_ret) << "\n";
  double diff = fp::RationalToDouble<D, N>(fr_ret) - double_ret;
  assert(std::abs(diff) < 0.001);
  if (std::abs(diff) >= 0.001) return false;

  std::cout << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << "num_variables: " << pb.num_variables() << "\n";

  return true;
}

inline bool TestConvolution() {
  Tick tick(__FN__);
  std::vector<bool> rets;
  {
    constexpr size_t DataCol = 4;
    constexpr size_t DataRow = 4;
    constexpr size_t ConvCol = 3;
    constexpr size_t ConvRow = 3;
    std::array<double, DataCol* DataRow> double_data = {
        0.666667, 0.203922, 0, 0, 0.996078, 0.54902,   0, 0,
        0.996078, 0.415686, 0, 0, 0.819608, 0.0705882, 0, 0};
    std::array<double, ConvCol* ConvRow + 1> double_para = {
        -1.24212,  -1.46405, -1.25206, 0.767129, 0.201014,
        -0.873962, 0.776801, 0.977157, 0.266173, 0.187466};
    bool ret = ConvGadget<4, 24, DataCol, DataRow, ConvCol, ConvRow>::Test(
        double_data, double_para);
    rets.push_back(ret);
  }
  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
};  // namespace circuit::cnn