#pragma once

#include "../details.h"
#include "../parallel_r1cs.h"
#include "./auxi_pub.h"
#include "./context.h"
#include "./prove.h"
#include "./publish.h"
#include "./verify.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

// --vgg16_infer "E:/code/crypto/pod_doc/vgg16_2/test_image e:/vgg16_publish"
bool TestInfer(std::string const& test_image_path,
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

inline bool TestConv(std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  ProveContext prove_context(working_path);
  VerifyContext verify_context(working_path + "/pub");

  bool all_success = false;
  std::vector<size_t> layers{
      0 /*, 2, 5, 7, 10, 12, 14, 17, 19, 21, 24, 26, 28*/};
  auto parallel_f = [&prove_context, &verify_context, &seed,
                     &layers](int64_t i) {
    auto layer = layers[i];
    OneConvProof proof;
    OneConvProve(seed, prove_context, layer, proof);
    return OneConvVerify(seed, verify_context, layer, proof);
  };
  parallel::For(&all_success, layers.size(), parallel_f);

  return all_success;
}

inline bool TestReluBn(std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  ProveContext prove_context(working_path);
  VerifyContext verify_context(working_path + "/pub");

  ReluBnProof proof;
  ReluBnProve(seed, prove_context, proof);
  return ReluBnVerify(seed, verify_context, proof);
}

inline bool TestPooling(std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  ProveContext prove_context(working_path);
  VerifyContext verify_context(working_path + "/pub");

  PoolingProof proof;
  PoolingProve(seed, prove_context, proof);
  return PoolingVerify(seed, verify_context, proof);
}

inline bool TestDense(std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  ProveContext prove_context(working_path);
  VerifyContext verify_context(working_path + "/pub");
  
  DenseProof proof0;
  DenseProve0(seed, prove_context, proof0);
  if (!DenseVerify0(seed, verify_context, proof0)) return false;

  DenseProof proof1;
  DenseProve1(seed, prove_context, proof1);
  if (!DenseVerify1(seed, verify_context, proof1)) return false;
  
  return true;
}
//
//inline bool TestAll(std::string const& working_path) {
//
//}

inline bool Test(std::string const& working_path) {
  bool ret = false;
  //ret = Publish("E:/code/crypto/pod_doc/vgg16_2/features", working_path);
  // ret =TestInfer("E:/code/crypto/pod_doc/vgg16_2/test_image", working_path);
  // ret =TestConv(working_path);
  //ret= TestReluBn(working_path);
  //ret= TestPooling(working_path);
  ret = TestDense(working_path);
  return ret;
}
}  // namespace clink::vgg16