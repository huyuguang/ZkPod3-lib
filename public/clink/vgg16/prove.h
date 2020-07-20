#pragma once

#include "./adapt.h"
#include "./image_com.h"
#include "./infer.h"
#include "./conv_prove.h"
#include "./dense_prove.h"
#include "./pooling_prove.h"
#include "./relubn_prove.h"

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
  hyrax::A4::Proof adapt_proof;
  clink::ParallelR1cs<R1cs>::Proof r1cs_proof;

  Proof(std::string const& file) {
    Tick tick(__FN__);
    if (!YasLoadBin(file, *this)) {
      throw std::invalid_argument("invalid proof file: " + file);
    }
  }
  Proof() {}

  bool operator==(Proof const& b) const {
    return conv == b.conv && relubn == b.relubn && pooling == b.pooling &&
           dense0 == b.dense0 && dense1 == b.dense1 &&
           adapt_proof == b.adapt_proof;
  }

  bool operator!=(Proof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1),
                       ("a", adapt_proof));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1),
                       ("a", adapt_proof));
  }
};

inline bool Prove(h256_t seed, dbl::Image const& test_image,
                  std::string const& working_path, Proof& proof) {
  Tick tick(__FN__);
  if (!InferAndCommit(test_image, working_path)) return false;

  ProveContext context(working_path);

  std::unique_ptr<AdaptProveItemMan> padapt_man(new AdaptProveItemMan);
  auto& adapt_man = *padapt_man;
  std::unique_ptr<R1csProveItemMan> pr1cs_man(new R1csProveItemMan);
  auto& r1cs_man = *pr1cs_man;

  //ParallelVoidTaskMan task_man;

  std::vector<parallel::VoidTask> tasks;

  // conv
  size_t conv_count = kConvLayers.size();
  for (size_t i = 0; i < conv_count; ++i) {
    tasks.emplace_back([&context, &seed, &proof, i, &adapt_man, &r1cs_man]() {
      OneConvProvePreprocess(seed, context, kConvLayers[i], proof.conv[i],
                             adapt_man, r1cs_man);
    });
  }
#if 1
  // relubn
  tasks.emplace_back([&context, &seed, &proof, &adapt_man, &r1cs_man]() {
    ReluBnProvePreprocess(seed, context, proof.relubn, adapt_man, r1cs_man);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof, &adapt_man, &r1cs_man]() {
    PoolingProvePreprocess(seed, context, proof.pooling, adapt_man, r1cs_man);
  });

  // dense0
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<0>(seed, context, proof.dense0);
  });

  // dense1
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<1>(seed, context, proof.dense1);
  });
#endif // #if 0

  {
    Tick subtick("preprocess");
    auto f1 = [&tasks](int64_t i) { tasks[i](); };
    parallel::For(tasks.size(), f1);
  }

  AdaptProve(seed, adapt_man, proof.adapt_proof);
  padapt_man.reset();

  R1csProve(seed, r1cs_man, proof.r1cs_proof);

  //std::vector<parallel::VoidTask> void_tasks;
  //task_man.take(void_tasks);
  //void_tasks.emplace_back([&seed, &adapt_man,&proof]() {
   // AdaptProve(seed, adapt_man, proof.adapt_proof);
  //});
  //{
  //  Tick subtick("r1cs prove");
  //  auto f2 = [&void_tasks](int64_t i) { void_tasks[i](); };
  //  parallel::For(void_tasks.size(), f2);
  //}
  return true;
}

}  // namespace clink::vgg16
