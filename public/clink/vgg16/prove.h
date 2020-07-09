#pragma once

#include "./adapt.h"
#include "./image_com.h"
#include "./infer.h"
#include "./prove_conv.h"
#include "./prove_dense.h"
#include "./prove_pooling.h"
#include "./prove_relubn.h"

namespace clink::vgg16 {

inline bool InferAndCommit(dbl::Image const& test_image,
                           std::string const& working_path) {
  Tick tick(__FN__);
  try {
    Para para(working_path + "/sec/para");
    std::array<std::unique_ptr<Image>, 35> images;
    Infer(para, test_image, images);
    for (size_t i = 0; i < images.size(); ++i) {
      std::string file = working_path + "/sec";
      fs::create_directories(file);
      file += "/image_" + std::to_string(i);
      fs::remove(file);
      yas::file_ostream os(file.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      oa.serialize(*images[i]);
    }

    ImageCommitmentPub image_com_pub;
    ImageCommitmentSec image_com_sec;
    ComputePerImageCommitment(images, image_com_pub, image_com_sec);
    YasSaveBin(working_path + "/pub/image_com_pub", image_com_pub);
    YasSaveBin(working_path + "/sec/image_com_sec", image_com_sec);
    return true;
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return false;
  }
}

struct Proof {
  std::array<OneConvProof, 13> conv;
  ReluBnProof relubn;
  PoolingProof pooling;
  DenseProof dense0;
  DenseProof dense1;

  Proof(std::string const& file) {
    Tick tick(__FN__);
    if (!YasLoadBin(file, *this)) {
      throw std::invalid_argument("invalid proof file: " + file);
    }
  }
  Proof() {}

  bool operator==(Proof const& b) const {
    return conv == b.conv && relubn == b.relubn && pooling == b.pooling &&
           dense0 == b.dense0 && dense1 == b.dense1;
  }

  bool operator!=(Proof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1));
  }
};

inline bool Prove(h256_t seed, dbl::Image const& test_image,
                  std::string const& working_path, Proof& proof) {
  Tick tick(__FN__);
  if (!InferAndCommit(test_image, working_path)) return false;

  ProveContext context(working_path);

  std::vector<std::function<void()>> tasks;

  // conv
  for (size_t i = 0; i < kConvLayers.size(); ++i) {
    tasks.emplace_back([&context, &seed, &proof, i]() {
      OneConvProve(seed, context, kConvLayers[i], proof.conv[i]);
    });
  }

  // relubn
  tasks.emplace_back([&context, &seed, &proof]() {
    ReluBnProve(seed, context, proof.relubn);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof]() {
    PoolingProve(seed, context, proof.pooling);
  });

  // dense0
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<0>(seed, context, proof.dense0);
  });

  // dense1
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<1>(seed, context, proof.dense1);
  });

  auto parallel_f = [&tasks](int64_t i) { tasks[i](); };
  parallel::For(tasks.size(), parallel_f);

  return true;
}

}  // namespace clink::vgg16