#pragma once

#include <memory>

#include "./context.h"
#include "./image_com.h"
#include "./r1cs_pub.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

enum ReluBnImageType { kCombinedInput, kCombinedOutput, kInput, kOutput };

struct ReluBnImage {
  ReluBnImageType type;
  std::vector<Fr> x;
  G1 com_x;
  Fr com_x_r;
  std::vector<Fr> a;
};

struct ReluBnInOutPub {
  std::vector<G1> cx;
  bool operator==(ReluBnInOutPub const& b) const { return cx == b.cx; }

  bool operator!=(ReluBnInOutPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnInOutPub", ("c", cx));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnInOutPub", ("c", cx));
  }
};

struct ReluBnInOutSec {
  std::vector<Fr> input;   // combined inputs
  std::vector<Fr> output;  // combined outputs
  Fr com_in_r;
  Fr com_out_r;
};

struct ReluBnR1csPub {
  std::vector<G1> com_w;

  bool operator==(ReluBnR1csPub const& b) const { return com_w == b.com_w; }

  bool operator!=(ReluBnR1csPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnR1csPub", ("c", com_w));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnR1csPub", ("c", com_w));
  }
};

struct ReluBnR1csSec : public BaseR1csSec {
  // nothing
};

struct ReluBnProof {
  ReluBnInOutPub io_pub;
  ReluBnR1csPub r1cs_pub;

  bool operator==(ReluBnProof const& b) const {
    return io_pub == b.io_pub && r1cs_pub == b.r1cs_pub;
  }

  bool operator!=(ReluBnProof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnProof", ("i", io_pub), ("r", r1cs_pub));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.ReluBnProof", ("i", io_pub), ("r", r1cs_pub));
  }
};

inline void ReluBnUpdateSeed(h256_t& seed, G1 const& x, G1 const& y) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, x);
  HashUpdate(hash, y);
  hash.Final(seed.data());
}

inline std::vector<Fr> ReluBnComputeFst(h256_t const& seed,
                                        std::string const& salt, size_t size) {
  std::vector<Fr> r(size);
  ComputeFst(seed, salt + std::to_string(size), r);
  return r;
}

inline void ReluBnBuildChallenge(h256_t const& seed,
                                 std::vector<ReluBnImage>& images) {
  Tick tick(__FN__);
  assert(images[0].type == ReluBnImageType::kCombinedInput);
  assert(images[1].type == ReluBnImageType::kCombinedOutput);

  images[0].a = ReluBnComputeFst(seed, "relubn input", images[0].x.size());
  images[1].a = ReluBnComputeFst(seed, "relubn output", images[1].x.size());

  auto it_i = images[0].a.begin();
  auto it_o = images[1].a.begin();
  for (size_t i = 1; i < images.size() / 2; ++i) {
    auto& input_image = images[2 * i];
    auto& output_image = images[2 * i + 1];
    assert(input_image.type == ReluBnImageType::kInput);
    assert(output_image.type == ReluBnImageType::kOutput);
    auto size = input_image.x.size();
    input_image.a.resize(size);
    output_image.a.resize(size);
    std::copy(it_i, it_i + size, input_image.a.begin());
    std::copy(it_o, it_o + size, output_image.a.begin());
    it_i += size;
    it_o += size;
  }
}

// 0: combined_in
// 1: combined_out
// even: in
// odd: out
inline void ReluBnBuildImages(ProveContext const& context,
                              std::vector<ReluBnImage>& images) {
  Tick tick(__FN__);
  auto const& const_images = context.const_images();
  auto const& com_images = context.image_com_pub().c;
  auto const& com_images_r = context.image_com_sec().r;

  std::vector<Fr> combined_in_x;
  std::vector<Fr> combined_out_x;
  images.resize(kReluBnLayers.size() * 2 + 2);
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    ReluBnImage& in = images[2 + order * 2];
    in.x = const_images[layer]->data;
    in.type = ReluBnImageType::kInput;
    in.com_x = com_images[layer];
    in.com_x_r = com_images_r[layer];
    combined_in_x.insert(combined_in_x.end(), in.x.begin(), in.x.end());

    ReluBnImage& out = images[2 + order * 2 + 1];
    out.x = const_images[layer + 1]->data;
    out.type = ReluBnImageType::kOutput;
    out.com_x = com_images[layer + 1];
    out.com_x_r = com_images_r[layer + 1];
    combined_out_x.insert(combined_out_x.end(), out.x.begin(), out.x.end());
  }

  ReluBnImage& combined_in = images[0];
  combined_in.type = ReluBnImageType::kCombinedInput;
  combined_in.x = std::move(combined_in_x);
  combined_in.com_x_r = FrRand();

  ReluBnImage& combined_out = images[1];
  combined_out.type = ReluBnImageType::kCombinedOutput;
  combined_out.x = std::move(combined_out_x);
  combined_out.com_x_r = FrRand();

  std::array<parallel::VoidTask, 2> tasks;
  tasks[0] = [&combined_in]() {
    combined_in.com_x =
        pc::PcComputeCommitmentG(combined_in.x, combined_in.com_x_r);
  };

  tasks[1] = [&combined_out]() {
    combined_out.com_x =
        pc::PcComputeCommitmentG(combined_out.x, combined_out.com_x_r);
  };

  parallel::Invoke(tasks);
}

inline size_t ReluBnGetCircuitCount() {
  size_t size = 0;
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    size_t C = kImageInfos[layer].C;
    size_t D = kImageInfos[layer].D;
    size += C * D * D;
  }
  return size;
}

inline std::string const& ReluBnR1csTag() {
  static const std::string kTag = "relubn r1cs";
  return kTag;
}

inline std::string const& ReluBnAdaptTag(bool in) {
  static const std::string kTagIn = "relubn in";
  static const std::string kTagOut = "relubn out";
  return in ? kTagIn : kTagOut;
}
}  // namespace clink::vgg16