#pragma once

#include <memory>

#include "./image_com.h"
#include "./para_com.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

struct PoolingInputImage {
  std::vector<Fr> const* x;
  size_t C;
  size_t D;
};

// 5 pooling layer + a + b + c + d
struct PoolingInputPub {
  std::array<G1, 9> cx;

  bool operator==(PoolingInputPub const& b) const { return cx == b.cx; }

  bool operator!=(PoolingInputPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.PoolingInputPub", ("c", cx));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.PoolingInputPub", ("c", cx));
  }
};

struct PoolingInputSec {
  std::vector<Fr> a;
  std::vector<Fr> b;
  std::vector<Fr> c;
  std::vector<Fr> d;
  Fr r_a;
  Fr r_b;
  Fr r_c;
  Fr r_d;
};

struct PoolingR1csPub {
  std::vector<G1> com_w;
  size_t r1cs_ret_index;

  bool operator==(PoolingR1csPub const& b) const {
    return com_w == b.com_w && r1cs_ret_index == b.r1cs_ret_index;
  }

  bool operator!=(PoolingR1csPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.PoolingR1csPub", ("c", com_w),
                       ("r", r1cs_ret_index));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.PoolingR1csPub", ("c", com_w),
                       ("r", r1cs_ret_index));
  }
};

struct PoolingR1csSec {
  std::vector<Fr> y;
  Fr ry;

  std::unique_ptr<R1csInfo> r1cs_info;
  std::vector<Fr> com_w_r;
  std::vector<std::vector<Fr>> mutable w;
};

struct PoolingOutputPub {
  std::array<G1, 6> cx;

  bool operator==(PoolingOutputPub const& b) const { return cx == b.cx; }

  bool operator!=(PoolingOutputPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.PoolingOutputPub", ("c", cx));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.PoolingOutputPub", ("c", cx));
  }
};

struct PoolingProof {
  PoolingInputPub input_pub;
  PoolingR1csPub r1cs_pub;
  PoolingOutputPub output_pub;
  clink::ParallelR1cs<R1cs>::Proof r1cs_proof;

  bool operator==(PoolingProof const& b) const {
    return r1cs_proof == b.r1cs_proof && input_pub == b.input_pub &&
           r1cs_pub == b.r1cs_pub && output_pub == b.output_pub;
  }

  bool operator!=(PoolingProof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.PoolingProof", ("i", input_pub), ("r", r1cs_pub),
                       ("o", output_pub), ("p", r1cs_proof));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.PoolingProof", ("i", input_pub), ("r", r1cs_pub),
                       ("o", output_pub), ("p", r1cs_proof));
  }
};

inline void PoolingUpdateSeed(h256_t& seed, G1 const& a, G1 const& b,
                              G1 const& c, G1 const& d) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, a);
  HashUpdate(hash, b);
  HashUpdate(hash, c);
  HashUpdate(hash, d);
  hash.Final(seed.data());
}

inline void PoolingComputeFst(h256_t const& seed, std::string const& prefix,
                              size_t order, std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(order);
  ComputeFst(seed, salt, r);
}

inline size_t PoolingGetCircuitCount() {
  size_t total_size = 0;
  for (size_t i = 0; i < kPoolingLayers.size(); ++i) {
    auto layer = kPoolingLayers[i];
    auto C = kImageInfos[layer].C;
    auto D = kImageInfos[layer].D;
    total_size += C * D * D;
  }
  return total_size / 4;
}
void PoolingImageToAbcd(std::array<PoolingInputImage, 5> const& images,
                        std::vector<Fr>& a, std::vector<Fr>& b,
                        std::vector<Fr>& c, std::vector<Fr>& d) {
  size_t total_size = 0;
  for (size_t i = 0; i < images.size(); ++i) {
    total_size += images[i].x->size();
  }

  a.resize(total_size / 4);
  b.resize(total_size / 4);
  c.resize(total_size / 4);
  d.resize(total_size / 4);

  size_t offset = 0;
  for (size_t l = 0; l < images.size(); ++l) {
    auto C = images[l].C;
    auto D = images[l].D;
    auto DD = D * D;
    auto D_2 = D / 2;
    auto DD_4 = DD / 4;
    auto const& x = *images[l].x;
    for (size_t i = 0; i < C; ++i) {
      for (size_t j = 0; j < D / 2; ++j) {
        for (size_t k = 0; k < D / 2; ++k) {
          auto idx = offset + i * DD_4 + j * D_2 + k;
          a[idx] = x[i * DD + 2 * j * D + 2 * k];
          b[idx] = x[i * DD + 2 * j * D + 2 * k + 1];
          c[idx] = x[i * DD + (2 * j + 1) * D + 2 * k];
          d[idx] = x[i * DD + (2 * j + 1) * D + 2 * k + 1];
        }
      }
    }
    offset += C * DD_4;
  }
}

inline void SelectInputComIpR(std::array<Fr, 9>& com_ip_r) {
  FrRand(com_ip_r.data(), 5);
  auto sum = std::accumulate(com_ip_r.begin(), com_ip_r.begin() + 5, FrZero());
  auto part = SplitFr(sum, 4);
  for (size_t i = 0; i < 4; ++i) {
    com_ip_r[i + 5] = part[i];
  }
}

inline void SelectOutputComIpR(std::array<Fr, 6>& com_ip_r) {
  com_ip_r[0] = FrRand();
  auto part = SplitFr(com_ip_r[0], 5);
  for (size_t i = 0; i < 5; ++i) {
    com_ip_r[i + 1] = part[i];
  }
}

inline void PoolingBuildInputQ(h256_t const& seed,
                               std::array<std::vector<Fr>, 9>& q) {
  std::array<PoolingInputImage, 5> images;
  for (size_t i = 0; i < 5; ++i) {
    auto layer = kPoolingLayers[i];
    images[i].C = kImageInfos[layer].C;
    images[i].D = kImageInfos[layer].D;
    q[i].resize(images[i].C * images[i].D * images[i].D);
    PoolingComputeFst(seed, "pooling input q", i, q[i]);
    images[i].x = &q[i];
  }

  PoolingImageToAbcd(images, q[5], q[6], q[7], q[8]);
}

inline void PoolingBuildAbcd(std::array<std::vector<Fr>, 9>& x) {
  std::array<PoolingInputImage, 5> images;
  for (size_t i = 0; i < 5; ++i) {
    auto layer = kPoolingLayers[i];
    images[i].x = &x[i];
    images[i].C = kImageInfos[layer].C;
    images[i].D = kImageInfos[layer].D;
  }
  PoolingImageToAbcd(images, x[5], x[6], x[7], x[8]);
}

inline void PoolingUpdateSeed(h256_t& seed, PoolingR1csPub const& pub) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(pub);
  auto buf = os.get_shared_buffer();
  HashUpdate(hash, buf.data.get(), buf.size);
  hash.Final(seed.data());
}

inline void PoolingBuildOutputQ(h256_t const& seed,
                                std::array<std::vector<Fr>, 6>& q) {
  std::vector<size_t> size;
  size_t total_size = 0;
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first != kPooling) continue;
    auto const& info = kImageInfos[i + 1];
    size.push_back(info.C * info.D * info.D);
    total_size += size.back();
  }
  q[0].resize(total_size);
  for (size_t i = 0; i < q.size() - 1; ++i) {
    q[i + 1].resize(size[i]);
  }

  ComputeFst(seed, "pooling output q", q[0]);
  std::vector<Fr>::const_iterator q_cur = q[0].begin();
  for (size_t i = 1; i < q.size(); ++i) {
    q[i].assign(q_cur, q_cur + q[i].size());
    q_cur += q[i].size();
  }
}

}  // namespace clink::vgg16
