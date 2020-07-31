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
  clink::ParallelR1cs<R1cs>::Proof r1cs_proof_conv;
  clink::ParallelR1cs<R1cs>::Proof r1cs_proof_misc;

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
           r1cs_proof_conv == b.r1cs_proof_conv &&
           r1cs_proof_misc == b.r1cs_proof_misc;
  }

  bool operator!=(Proof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1),
                       ("a1", r1cs_proof_conv),("a2", r1cs_proof_misc));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.Proof", ("c", conv), ("r", relubn),
                       ("p", pooling), ("d0", dense0), ("d1", dense1),
                       ("a1", r1cs_proof_conv),("a2", r1cs_proof_misc));
  }
};

inline bool Prove(h256_t seed, dbl::Image const& test_image,
                  std::string const& working_path, Proof& proof) {
  Tick tick(__FN__);
  if (!InferAndCommit(test_image, working_path)) return false;

  ProveContext context(working_path);

  std::unique_ptr<AdaptProveItemMan> padapt_man(new AdaptProveItemMan);
  auto& adapt_man = *padapt_man;

  std::unique_ptr<R1csProveItemMan> pr1cs_man_conv(new R1csProveItemMan);
  auto& r1cs_man_conv = *pr1cs_man_conv;

  std::unique_ptr<R1csProveItemMan> pr1cs_man_misc(new R1csProveItemMan);
  auto& r1cs_man_misc = *pr1cs_man_misc;

  //ParallelVoidTaskMan task_man;

  std::vector<parallel::VoidTask> tasks;

  // conv
  size_t conv_count = kConvLayers.size();
  for (size_t i = 0; i < conv_count; ++i) {
    tasks.emplace_back([&context, &seed, &proof, i, &adapt_man, &r1cs_man_conv]() {
      OneConvProvePreprocess(seed, context, kConvLayers[i], proof.conv[i],
                             adapt_man, r1cs_man_conv);
    });
  }
#if 1
  // relubn
  tasks.emplace_back([&context, &seed, &proof, &adapt_man, &r1cs_man_misc]() {
    ReluBnProvePreprocess(seed, context, proof.relubn, adapt_man, r1cs_man_misc);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof, &adapt_man, &r1cs_man_misc]() {
    PoolingProvePreprocess(seed, context, proof.pooling, adapt_man, r1cs_man_misc);
  });

  // dense0
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<0>(seed, context, proof.dense0);
  });

  // dense1
  tasks.emplace_back([&context, &seed, &proof]() {
    DenseProve<1>(seed, context, proof.dense1);
  });
#endif

  {
    Tick subtick("preprocess");
    auto f1 = [&tasks](int64_t i) { tasks[i](); };
    parallel::For(tasks.size(), f1);
  }  

  std::vector<parallel::VoidTask> void_tasks;
  void_tasks.emplace_back([&seed, &adapt_man, &padapt_man, &proof]() {
    // TODO: for debug AdaptProve(seed, adapt_man, proof.adapt_proof);
    padapt_man.reset();
  });

  void_tasks.emplace_back([&seed, &r1cs_man_conv, &proof]() {
    // TODO: for debug R1csProve(seed, r1cs_man_conv, proof.r1cs_proof_conv);
  });

  void_tasks.emplace_back([&seed, &r1cs_man_misc, &proof]() {
    R1csProve(seed, r1cs_man_misc, proof.r1cs_proof_misc);
  });

  auto f2 = [&void_tasks](int64_t i) { void_tasks[i](); };
  parallel::For(void_tasks.size(), f2, true);  // TODO: true->false

  return true;
}

}  // namespace clink::vgg16
