#pragma once

#include <memory>

#include "./adapt.h"
#include "./r1cs.h"
#include "./context.h"
#include "./image_com.h"
#include "./policy.h"
#include "./r1cs_pub.h"
#include "circuit/vgg16/vgg16.h"
#include "clink/equal_ip.h"
#include "clink/equality2.h"

namespace clink::vgg16 {
inline void OneConvUpdateSeed(h256_t& seed, std::array<G1, 9> const& c) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  for (auto const& i : c) {
    HashUpdate(hash, i);
  }
  hash.Final(seed.data());
}

inline void OneConvUpdateSeed(h256_t& seed, G1 const& c) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, c);
  hash.Final(seed.data());
}

struct OneConvInputPub {
  std::array<G1, 9> cb;

  bool operator==(OneConvInputPub const& b) const { return cb == b.cb; }

  bool operator!=(OneConvInputPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.OneConvInputPub", ("cb", cb));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.OneConvInputPub", ("cb", cb));
  }
};

struct OneConvInputSec {
  std::vector<std::array<Fr const*, 9>> b;
  std::array<Fr, 9> rb;
};

struct OneConvR1csPub {
  std::vector<G1> com_w;
  size_t r1cs_ret_index;

  bool operator==(OneConvR1csPub const& b) const {
    return com_w == b.com_w && r1cs_ret_index == b.r1cs_ret_index;
  }

  bool operator!=(OneConvR1csPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.OneConvR1csPub", ("c", com_w),
                       ("r", r1cs_ret_index));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.OneConvR1csPub", ("c", com_w),
                       ("r", r1cs_ret_index));
  }
};

struct OneConvR1csSec : public BaseR1csSec {
  std::vector<Fr> y;
  Fr ry;
};

struct OneConvOutputPub {
  G1 cy;

  bool operator==(OneConvOutputPub const& b) const { return cy == b.cy; }

  bool operator!=(OneConvOutputPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.OneConvOutputPub", ("c", cy));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.OneConvOutputPub", ("c", cy));
  }
};

struct OneConvProof {
  OneConvInputPub input_pub;
  OneConvR1csPub r1cs_pub;
  OneConvOutputPub output_pub;

  bool operator==(OneConvProof const& b) const {
    return input_pub == b.input_pub &&
           r1cs_pub == b.r1cs_pub && output_pub == b.output_pub;
  }

  bool operator!=(OneConvProof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.OneConvProof", ("i", input_pub), ("r", r1cs_pub),
                       ("o", output_pub));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.OneConvProof", ("i", input_pub), ("r", r1cs_pub),
                       ("o", output_pub));
  }
};

inline void OneConvComputeFst(h256_t const& seed, std::string const& prefix,
                              size_t layer, size_t col, std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(layer);
  salt += " ";
  salt += std::to_string(col);
  ComputeFst(seed, salt, r);
}

inline void OneConvComputeFst(h256_t const& seed, std::string const& prefix,
                              std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(r.size());
  ComputeFst(seed, salt, r);
}

inline void OneConvComputeInputR(h256_t const& seed, size_t layer, size_t KCDD,
                                 std::array<std::vector<Fr>, 9>& r) {
  auto parallel_f = [&seed, &r, layer, KCDD](int64_t j) {
    r[j].resize(KCDD);
    OneConvComputeFst(seed, "conv adapt input ", layer, j, r[j]);
  };
  parallel::For(9, parallel_f);
}

inline std::vector<Fr> OneConvComputeOutputR(h256_t const& seed, size_t layer,
                                             size_t KDD) {
  std::vector<Fr> r(KDD);
  std::string salt = "conv output ";
  salt += std::to_string(layer);
  OneConvComputeFst(seed, salt, r);
  return r;
}

inline std::vector<Fr> OneConvOutputR2S(size_t K, size_t C, size_t D,
                                        std::vector<Fr> const& r) {
  assert(r.size() == K * D * D);
  std::vector<Fr> s(K * C * D * D);
  for (size_t i = 0; i < K; ++i) {
    for (size_t j = 0; j < C; ++j) {
      for (size_t k = 0; k < D; ++k) {
        for (size_t l = 0; l < D; ++l) {
          s[i * C * D * D + j * D * D + k * D + l] = r[i * D * D + k * D + l];
        }
      }
    }
  }
  return s;
}

inline std::string ConvR1csTag(size_t layer) {
  static const std::string tag = "conv r1cs ";
  return tag + std::to_string(layer);
}

inline std::string ConvAdaptTag(bool in, size_t layer) {
  static const std::string in_tag = "conv in ";
  static const std::string out_tag = "conv out ";
  return (in? in_tag:out_tag)+std::to_string(layer);
}
}  // namespace clink::vgg16
