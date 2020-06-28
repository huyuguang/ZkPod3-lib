#pragma once

#include "./para_fr.h"
#include "./auxi_pub.h"
#include "ecc/ecc.h"
#include "ecc/pc_base.h"
#include "misc/misc.h"
#include "public.h"
#include "utils/fst.h"

namespace clink::vgg16 {

struct BnCommitmentPub {  
  G1 mu;
  G1 alpha;
  G1 beta;

  bool operator==(BnCommitmentPub const& b) const {
    return mu == b.mu && alpha == b.alpha && beta == b.beta;
  }

  bool operator!=(BnCommitmentPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.bn", ("m", mu), ("a", alpha),
                       ("b", beta));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.bn", ("m", mu), ("a", alpha),
                       ("b", beta));
  }
};

struct BnCommitmentSec {
  Fr mu_r;
  Fr alpha_r;
  Fr beta_r;
  bool operator==(BnCommitmentSec const& b) const {
    return mu_r == b.mu_r && alpha_r == b.alpha_r && beta_r == b.beta_r;
  }

  bool operator!=(BnCommitmentSec const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.bn", ("m", mu_r), ("a", alpha_r),
                       ("b", beta_r));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.bn", ("m", mu_r), ("a", alpha_r),
                       ("b", beta_r));
  }
};

struct ConvCommitmentPub {
  std::array<std::array<G1, 9>, 13> coef;
  std::array<G1, 13> bias;

  bool operator==(ConvCommitmentPub const& b) const {
    return coef == b.coef && bias == b.bias;
  }

  bool operator!=(ConvCommitmentPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.conv", ("c", coef), ("b", bias));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.conv", ("c", coef), ("b", bias));
  }
};

struct ConvCommitmentSec {
  std::array<std::array<Fr, 9>, 13> coef_r;
  std::array<Fr, 13> bias_r;

  bool operator==(ConvCommitmentSec const& b) const {
    return coef_r == b.coef_r && bias_r == b.bias_r;
  }

  bool operator!=(ConvCommitmentSec const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.conv", ("c", coef_r), ("b", bias_r));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.conv", ("c", coef_r), ("b", bias_r));
  }
};

struct DenseCommitmentPub {
  std::array<G1, 512> d0;
  std::array<G1, 10> d1;

  bool operator==(DenseCommitmentPub const& b) const {
    return d0 == b.d0 && d1 == b.d1;
  }

  bool operator!=(DenseCommitmentPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.dense", ("d0", d0), ("d1", d1));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.compub.dense", ("d0", d0), ("d1", d1));
  }
};

struct DenseCommitmentSec {
  std::array<Fr, 512> d0_r;
  std::array<Fr, 10> d1_r;

  bool operator==(DenseCommitmentSec const& b) const {
    return d0_r == b.d0_r && d1_r == b.d1_r;
  }

  bool operator!=(DenseCommitmentSec const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.dense", ("d0", d0_r), ("d1", d1_r));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec.dense", ("d0", d0_r), ("d1", d1_r));
  }
};

struct ParaCommitmentPub {
  BnCommitmentPub bn;
  ConvCommitmentPub conv;
  DenseCommitmentPub dense;

  ParaCommitmentPub() {}
  ParaCommitmentPub(std::string const& file) {
    if (!Load(file)) {
      throw std::invalid_argument("invalid para commitment pub file: " + file);
    }
  }

  bool operator==(ParaCommitmentPub const& b) const {
    return bn == b.bn && conv == b.conv && dense == b.dense;
  }

  bool operator!=(ParaCommitmentPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.compub", ("b", bn), ("c", conv),
                       ("d", dense));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.compub", ("b", bn), ("c", conv),
                       ("d", dense));
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

  bool Save(std::string const& file) const {
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
    ParaCommitmentPub check;
    if (!check.Load(file)) {
      std::cout << "oops\n";
      return false;
    }
    if (check != *this) {
      std::cout << "oops\n";
      return false;
    }
#endif

    return true;
  }
};

struct ParaCommitmentSec {
  BnCommitmentSec bn;
  ConvCommitmentSec conv;
  DenseCommitmentSec dense;

  ParaCommitmentSec() {}
  ParaCommitmentSec(std::string const& file) {
    if (!Load(file)) {
      throw std::invalid_argument("invalid para commitment sec file: " + file);
    }
  }

  bool operator==(ParaCommitmentSec const& b) const {
    return bn == b.bn && conv == b.conv && dense == b.dense;
  }

  bool operator!=(ParaCommitmentSec const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec", ("b", bn), ("c", conv),
                       ("d", dense));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para.comsec", ("b", bn), ("c", conv),
                       ("d", dense));
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

  bool Save(std::string const& file) const {
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
    ParaCommitmentSec check;
    if (!check.Load(file)) {
      std::cout << "oops\n";
      return false;
    }
    if (check != *this) {
      std::cout << "oops\n";
      return false;
    }
#endif

    return true;
  }
};

inline void ComputeBnCommitment(std::array<Para::BnLayer, 14> const& para,
                                AuxiPub const& auxi, BnCommitmentPub& pub,
                                BnCommitmentSec& sec) {
  Tick tick(__FN__);  
  sec.alpha_r = FrRand();
  sec.beta_r = FrRand();
  sec.mu_r = FrRand();

  // combine alpha,beta,mu
  std::vector<Fr> all_alpha;
  all_alpha.reserve(4736);
  std::vector<Fr> all_beta;
  all_beta.reserve(4736);
  std::vector<Fr> all_mu;
  all_mu.reserve(4736);

  for (auto const& i : para) {    
    all_alpha.insert(all_alpha.end(), i.alpha.begin(), i.alpha.end());
    all_beta.insert(all_beta.end(), i.beta.begin(), i.beta.end());
    all_mu.insert(all_mu.end(), i.mu.begin(), i.mu.end());
  }

#ifdef _DEBUG_CHECK
  if (all_alpha.size() != 4736) throw std::runtime_error("oops");
  if (all_beta.size() != 4736) throw std::runtime_error("oops");
  if (all_mu.size() != 4736) throw std::runtime_error("oops");
#endif

  pub.alpha = pc::PcComputeCommitmentG(all_alpha.size(),auxi.para_u_bn().first,
                                       all_alpha.data(),sec.alpha_r);

  pub.beta = pc::PcComputeCommitmentG(all_beta.size(),auxi.para_u_bn().first,
                                       all_beta.data(),sec.beta_r);

  pub.mu = pc::PcComputeCommitmentG(all_mu.size(),auxi.para_u_bn().first,
                                       all_mu.data(),sec.mu_r);

#ifdef _DEBUG_CHECK
  std::vector<Fr> extended_alpha;
  std::vector<Fr> extended_beta;
  std::vector<Fr> extended_mu;
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first != kReluBn) continue;
    size_t C = kImageInfos[i].channel_count;
    size_t D = kImageInfos[i].dimension;
    auto order = kLayerTypeOrders[i].second;
    auto const& bn_para = para[order];
    if (C != bn_para.alpha.size()) throw std::runtime_error("oops");
    for (size_t j = 0; j < bn_para.alpha.size(); ++j) {
      // repeat DD times
      for (size_t k = 0; k < D*D; ++k) {
        extended_alpha.push_back(bn_para.alpha[j]);
        extended_beta.push_back(bn_para.beta[j]);
        extended_mu.push_back(bn_para.mu[j]);
      }
    }    
  }
  if (pub.alpha != pc::PcComputeCommitmentG(extended_alpha, sec.alpha_r)) {
    throw std::runtime_error("oops");
  }
  if (pub.beta != pc::PcComputeCommitmentG(extended_beta, sec.beta_r)) {
    throw std::runtime_error("oops");
  }
  if (pub.mu != pc::PcComputeCommitmentG(extended_mu, sec.mu_r)) {
    throw std::runtime_error("oops");
  }
#endif
}

inline void ComputeConvCommitment(std::array<Para::ConvLayer, 13> const& para,
                                  AuxiPub const& auxi, ConvCommitmentPub& pub,
                                  ConvCommitmentSec& sec) {
  Tick tick(__FN__);
  auto parallel_f = [&para, &auxi, &pub, &sec](int64_t o) {
    auto C = kConvLayerInfos[o].channel_count;
    auto K = kConvLayerInfos[o].kernel_count;

    auto const& layer = para[o];
    auto& pub_coef = pub.coef[o];
    auto& pub_bias = pub.bias[o];
    auto& sec_coef = sec.coef_r[o];
    auto& sec_bias = sec.bias_r[o];

    FrRand(sec_coef.data(), sec_coef.size());
    sec_bias = FrRand();

    auto get_coef_u = [&auxi, o](int64_t i) -> G1 const& {
      auto range = auxi.para_u_conv_coef(o);
      return i ? range.first[i - 1] : pc::PcH();
    };

    auto parallel_c = [&pub_coef, &get_coef_u, &layer, &sec_coef, K,
                       C](int64_t j) {
      auto const& coefs = layer.coefs;
      auto const& coefs_r = sec_coef[j];
      auto get_coef = [&coefs, &coefs_r, j, K, C](int64_t i) -> Fr const& {
        return i ? coefs[(i - 1) / C][(i - 1) % C][j / 3][j % 3] : coefs_r;
      };
      pub_coef[j] = MultiExpBdlo12<G1>(get_coef_u, get_coef, K * C + 1);
    };
    parallel::For(9, parallel_c);

    auto get_bias_u = [&auxi, o](int64_t i) -> G1 const& {
      auto range = auxi.para_u_conv_bias(o);
      return i ? range.first[i - 1] : pc::PcH();
    };

    auto const& bias = layer.bias;
    auto get_bias = [&bias, &sec_bias](int64_t i) -> Fr const& {
      return i ? bias[i - 1] : sec_bias;
    };
    pub_bias = MultiExpBdlo12<G1>(get_bias_u, get_bias, K + 1);

#ifdef _DEBUG_CHECK
    // TODO: check com u
#endif
  };

  parallel::For((int64_t)para.size(), parallel_f);
}

inline void ComputeDenseCommitment(std::array<Para::DenseLayer, 2> const& para,
                                   DenseCommitmentPub& pub,
                                   DenseCommitmentSec& sec) {
  Tick tick(__FN__);
  FrRand(sec.d0_r.data(), sec.d0_r.size());
  FrRand(sec.d1_r.data(), sec.d1_r.size());
  
  assert(para[0].weight.size() == pub.d0.size());
  for (size_t i = 0; i < para[0].weight.size(); ++i) {
    auto const& w = para[0].weight[i];
    pub.d0[i] = pc::PcComputeCommitmentG(w, sec.d0_r[i]);
  }

  assert(para[1].weight.size() == pub.d1.size());
  for (size_t i = 0; i < para[1].weight.size(); ++i) {
    auto const& w = para[1].weight[i];
    pub.d1[i] = pc::PcComputeCommitmentG(w, sec.d1_r[i]);
  }
}

inline void ComputeParaCommitment(Para const& para, AuxiPub const& auxi,
                                  ParaCommitmentPub& pub,
                                  ParaCommitmentSec& sec) {
  Tick tick(__FN__);
  ComputeBnCommitment(para.bn_layers(), auxi, pub.bn, sec.bn);
  ComputeConvCommitment(para.conv_layers(), auxi, pub.conv, sec.conv);
  ComputeDenseCommitment(para.dense_layers(), pub.dense, sec.dense);
}
}  // namespace clink::vgg16