#pragma once

#include "./para_fr.h"
#include "ecc/ecc.h"
#include "ecc/pc_base.h"
#include "misc/misc.h"
#include "public.h"
#include "utils/fst.h"

namespace clink::vgg16 {

class AuxiPub {
 public:
  AuxiPub() {
    Tick tick(__FN__);
    InitPtr();

    ComputeParaBn(*para_u_bn_);    
    
    ComputeParaConv<32, 64, 64>(*para_u_conv1_);    
    ComputeParaConv<16, 128, 128>(*para_u_conv3_);    
    ComputeParaConv<8, 256, 256>(*para_u_conv6_);    
    ComputeParaConv<4, 512, 512>(*para_u_conv9_);    
    ComputeParaConv<2, 512, 512>(*para_u_conv12_);

    ComputeDataConv<32, 3, 64>(*data_u_conv0_);
    ComputeDataConv<32, 64, 64>(*data_u_conv1_);
    ComputeDataConv<16, 64, 128>(*data_u_conv2_);
    ComputeDataConv<16, 128, 128>(*data_u_conv3_);
    ComputeDataConv<8, 128, 256>(*data_u_conv4_);
    ComputeDataConv<8, 256, 256>(*data_u_conv6_);
    ComputeDataConv<4, 256, 512>(*data_u_conv7_);
    ComputeDataConv<4, 512, 512>(*data_u_conv9_);
    ComputeDataConv<2, 512, 512>(*data_u_conv12_);
  }

  AuxiPub(std::string const& file) {
    Tick tick(__FN__);
    InitPtr();

    if (!Load(file)) {
      throw std::invalid_argument("invalid auxi file: " + file);
    }
  }

  bool Save(std::string const& file) const {
    Tick tick(__FN__);
    try {
      boost::system::error_code dummy;
      fs::remove(file);
      yas::file_ostream os(file.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      oa.serialize(*this);
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }

#ifdef _DEBUG_CHECK
    try {
      AuxiPub check(file);
      if (check != *this) {
        std::cout << "oops\n";
        return false;
      }
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
#endif

    return true;
  }

  bool operator==(AuxiPub const& b) const {
    return *para_u_bn_ == *b.para_u_bn_ && *para_u_conv1_ == *b.para_u_conv1_ &&
           *para_u_conv3_ == *b.para_u_conv3_ &&
           *para_u_conv6_ == *b.para_u_conv6_ &&
           *para_u_conv9_ == *b.para_u_conv9_ &&
           *para_u_conv12_ == *b.para_u_conv12_ &&
           *data_u_conv0_ == *b.data_u_conv0_ &&
           *data_u_conv1_ == *b.data_u_conv1_ &&
           *data_u_conv2_ == *b.data_u_conv2_ &&
           *data_u_conv3_ == *b.data_u_conv3_ &&
           *data_u_conv4_ == *b.data_u_conv4_ &&
           *data_u_conv6_ == *b.data_u_conv6_ &&
           *data_u_conv7_ == *b.data_u_conv7_ &&
           *data_u_conv9_ == *b.data_u_conv9_ &&
           *data_u_conv12_ == *b.data_u_conv12_;
  }

  bool operator!=(AuxiPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP(
        "vgg16.auxi", ("b", *para_u_bn_), ("c1", *para_u_conv1_),
        ("c3", *para_u_conv3_), ("c6", *para_u_conv6_), ("c9", *para_u_conv9_),
        ("c12", *para_u_conv12_), ("d0", *data_u_conv0_),
        ("d1", *data_u_conv1_), ("d2", *data_u_conv2_), ("d3", *data_u_conv3_),
        ("d4", *data_u_conv4_), ("d6", *data_u_conv6_), ("d7", *data_u_conv7_),
        ("d9", *data_u_conv9_), ("d12", *data_u_conv12_));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP(
        "vgg16.auxi", ("b", *para_u_bn_), ("c1", *para_u_conv1_),
        ("c3", *para_u_conv3_), ("c6", *para_u_conv6_), ("c9", *para_u_conv9_),
        ("c12", *para_u_conv12_), ("d0", *data_u_conv0_),
        ("d1", *data_u_conv1_), ("d2", *data_u_conv2_), ("d3", *data_u_conv3_),
        ("d4", *data_u_conv4_), ("d6", *data_u_conv6_), ("d7", *data_u_conv7_),
        ("d9", *data_u_conv9_), ("d12", *data_u_conv12_));
  }

  // para bn
  std::pair<G1 const*, G1 const*> para_u_bn() const {
    return std::make_pair(para_u_bn_->data(), para_u_bn_->data() + para_u_bn_->size());
  }

  // para conv coef
  std::pair<G1 const*, G1 const*> para_u_conv_coef(size_t order) const {
    return std::make_pair(para_u_conv_ptr_[order],
                          para_u_conv_ptr_[order] + para_u_conv_coef_size_[order]);
  }

  // para conv bias
  std::pair<G1 const*, G1 const*> para_u_conv_bias(size_t order) const {
    return std::make_pair(para_u_conv_ptr_[order],
                          para_u_conv_ptr_[order] + para_u_conv_bias_size_[order]);
  }

  // data conv
  std::pair<G1 const*, G1 const*> data_u_conv(size_t order) const {
    return std::make_pair(data_u_conv_ptr_[order],
                          data_u_conv_ptr_[order] + data_u_conv_size_[order]);
  }

 private:
  std::unique_ptr<std::array<G1, 4736>> para_u_bn_;

  // $u_i=\prod_{j=0}^{DD-1}g_{iDD+j},i\in[0,KC-1]$
  std::unique_ptr<std::array<G1, 64 * 64>> para_u_conv1_;
  std::unique_ptr<std::array<G1, 128 * 128>> para_u_conv3_;
  std::unique_ptr<std::array<G1, 256 * 256>> para_u_conv6_;
  std::unique_ptr<std::array<G1, 512 * 512>> para_u_conv9_;
  std::unique_ptr<std::array<G1, 512 * 512>> para_u_conv12_;
  std::array<G1 const*, 13> para_u_conv_ptr_;
  std::array<size_t, 13> para_u_conv_coef_size_;
  std::array<size_t, 13> para_u_conv_bias_size_;

  // $u_i=\prod_{j=0}^{K-1}g_{i+jCDD},i\in[0,CDD-1]$
  std::unique_ptr<std::array<G1, 3 * 32 * 32>> data_u_conv0_;
  std::unique_ptr<std::array<G1, 64 * 32 * 32>> data_u_conv1_;
  std::unique_ptr<std::array<G1, 64 * 16 * 16>> data_u_conv2_;
  std::unique_ptr<std::array<G1, 128 * 16 * 16>> data_u_conv3_;
  std::unique_ptr<std::array<G1, 128 * 8 * 8>> data_u_conv4_;
  std::unique_ptr<std::array<G1, 256 * 8 * 8>> data_u_conv6_;
  std::unique_ptr<std::array<G1, 256 * 4 * 4>> data_u_conv7_;
  std::unique_ptr<std::array<G1, 512 * 4 * 4>> data_u_conv9_;
  std::unique_ptr<std::array<G1, 512 * 2 * 2>> data_u_conv12_;
  std::array<G1 const*, 13> data_u_conv_ptr_;
  std::array<size_t, 13> data_u_conv_size_;

 private:
  void InitPtr() {
    para_u_bn_.reset(new std::array<G1, 4736>);
    para_u_conv1_.reset(new std::array<G1, 64 * 64>);
    para_u_conv3_.reset(new std::array<G1, 128 * 128>);
    para_u_conv6_.reset(new std::array<G1, 256 * 256>);
    para_u_conv9_.reset(new std::array<G1, 512 * 512>);
    para_u_conv12_.reset(new std::array<G1, 512 * 512>);

    data_u_conv0_.reset(new std::array<G1, 3 * 32 * 32>);
    data_u_conv1_.reset(new std::array<G1, 64 * 32 * 32>);
    data_u_conv2_.reset(new std::array<G1, 64 * 16 * 16>);
    data_u_conv3_.reset(new std::array<G1, 128 * 16 * 16>);
    data_u_conv4_.reset(new std::array<G1, 128 * 8 * 8>);
    data_u_conv6_.reset(new std::array<G1, 256 * 8 * 8>);
    data_u_conv7_.reset(new std::array<G1, 256 * 4 * 4>);
    data_u_conv9_.reset(new std::array<G1, 512 * 4 * 4>);
    data_u_conv12_.reset(new std::array<G1, 512 * 2 * 2>);

    para_u_conv_ptr_ = std::array<G1 const*, 13>{
        {para_u_conv1_->data(), para_u_conv1_->data(), para_u_conv3_->data(),
         para_u_conv3_->data(), para_u_conv6_->data(), para_u_conv6_->data(),
         para_u_conv6_->data(), para_u_conv9_->data(), para_u_conv9_->data(),
         para_u_conv9_->data(), para_u_conv12_->data(), para_u_conv12_->data(),
         para_u_conv12_->data()}};
    para_u_conv_coef_size_ = std::array<size_t, 13>{
        {192, 4096, 8192, 16384, 32768, 65536, 65536, 131072, 262144, 262144,
         262144, 262144, 262144}};
    para_u_conv_bias_size_ = std::array<size_t, 13>{
        {64, 64, 128, 128, 256, 256, 256, 512, 512, 512, 512, 512, 512}};

    data_u_conv_ptr_ = std::array<G1 const*, 13>{
        {data_u_conv0_->data(), data_u_conv1_->data(), data_u_conv2_->data(),
         data_u_conv3_->data(), data_u_conv4_->data(), data_u_conv6_->data(),
         data_u_conv6_->data(), data_u_conv7_->data(), data_u_conv9_->data(),
         data_u_conv9_->data(), data_u_conv12_->data(), data_u_conv12_->data(),
         data_u_conv12_->data()}};
    data_u_conv_size_ = std::array<size_t, 13>{
        {data_u_conv0_->size(), data_u_conv1_->size(), data_u_conv2_->size(),
         data_u_conv3_->size(), data_u_conv4_->size(), data_u_conv6_->size(),
         data_u_conv6_->size(), data_u_conv7_->size(), data_u_conv9_->size(),
         data_u_conv9_->size(), data_u_conv12_->size(), data_u_conv12_->size(),
         data_u_conv12_->size()}};
  }

  void ComputeParaBn(std::array<G1, 4736>& u) {
    size_t constexpr kCount = 14;
    std::array<size_t, kCount> g_offsets = {
        {0, 65536, 131072, 163840, 196608, 212992, 229376, 245760, 253952,
         262144, 270336, 272384, 274432, 276480}};
    std::array<size_t, kCount> u_offsets = {{0, 64, 128, 256, 384, 640, 896,
                                             1152, 1664, 2176, 2688, 3200, 3712,
                                             4224}};
    std::array<size_t, kCount> range_i = {
        {64, 64, 128, 128, 256, 256, 256, 512, 512, 512, 512, 512, 512, 512}};
    std::array<size_t, kCount> range_j = {
        {1024, 1024, 256, 256, 64, 64, 64, 16, 16, 16, 4, 4, 4, 1}};

    for (size_t k = 0; k < kCount; ++k) {
      for (size_t i = 0; i < range_i[k]; ++i) {
        u[i + u_offsets[k]] =
            pc::PcComputeSigmaG(g_offsets[k] + i * range_j[k], range_j[k]);
      }
    }
  }

  template <size_t D, size_t C, size_t K>
  void ComputeParaConv(std::array<G1, K * C>& u) {
    auto DD = D * D;
    for (size_t i = 0; i < K * C; ++i) {
      u[i] = pc::PcComputeSigmaG(i * DD, DD);
    }
  }

  template <size_t D, size_t C, size_t K>
  void ComputeDataConv(std::array<G1, C * D * D>& u) {
    for (size_t i = 0; i < C * D * D; ++i) {
      u[i] = G1Zero();
      for (size_t j = 0; j < K; ++j) {
        u[i] += pc::PcG(i + j * C * D * D);
      }
    }
  }

  bool Load(std::string const& file) {
    try {
      yas::file_istream is(file.c_str());
      yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
      ia.serialize(*this);
      return true;
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
  }
};
};  // namespace clink::vgg16