#pragma once

#include "./image_com.h"
#include "./infer.h"
#include "./prove_conv.h"
#include "./prove_dense.h"
#include "./prove_pooling.h"
#include "./prove_relubn.h"

namespace clink::vgg16 {

inline bool InferAndCommit(std::string const& test_image_path,
                           std::string const& working_path) {
  Tick tick(__FN__);
  try {
    dbl::Image test_image(kImageInfos[0]);
    dbl::LoadTestImage(test_image_path, test_image);

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
    image_com_pub.Save(working_path + "/pub/image_com_pub");
    image_com_sec.Save(working_path + "/sec/image_com_sec");
    return true;
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return false;
  }
}

struct Proof {
  std::array<std::unique_ptr<OneConvProof>, 13> conv;
  std::unique_ptr<ReluBnProof> relubn;
  std::unique_ptr<PoolingProof> pooling;
  std::unique_ptr<DenseProof> dense0;
  std::unique_ptr<DenseProof> dense1;
};

inline bool Prove(h256_t seed, std::string const& test_image_path,
                  std::string const& working_path, Proof& proof) {
  Tick tick(__FN__);
  if (!InferAndCommit(test_image_path, working_path)) return false;

  ProveContext context(working_path);

  std::vector<std::function<void()>> tasks;

  // conv
  for (size_t i = 0; i < kConvLayers.size(); ++i) {
    tasks.emplace_back([&context, &seed, &proof, i]() {
      proof.conv[i].reset(new OneConvProof);
      OneConvProve(seed, context, kConvLayers[i], *proof.conv[i]);
    });
  }

  // relubn
  tasks.emplace_back([&context, &seed, &proof]() {
    proof.relubn.reset(new ReluBnProof);
    ReluBnProve(seed, context, *proof.relubn);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof]() {
    proof.pooling.reset(new PoolingProof);
    PoolingProve(seed, context, *proof.pooling);
  });

  // dense0
  tasks.emplace_back([&context, &seed, &proof]() {
    proof.dense0.reset(new DenseProof);
    DenseProve<0>(seed, context, *proof.dense0);
  });

  // dense1
  tasks.emplace_back([&context, &seed, &proof]() {
    proof.dense1.reset(new DenseProof);
    DenseProve<1>(seed, context, *proof.dense1);
  });

  auto parallel_f = [&tasks](int64_t i) { tasks[i](); };
  parallel::For(tasks.size(), parallel_f);

  return true;
}

}  // namespace clink::vgg16