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
inline bool TestInfer(std::string const& test_image_path,
               std::string const& working_path) {
  return InferAndCommit(test_image_path, working_path);
}

inline bool TestConv(std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  ProveContext prove_context(working_path);
  VerifyContext verify_context(working_path + "/pub");

  bool all_success = false;
  auto parallel_f = [&prove_context, &verify_context, &seed](int64_t i) {
    auto layer = kConvLayers[i];
    OneConvProof proof;
    OneConvProve(seed, prove_context, layer, proof);
    return OneConvVerify(seed, verify_context, layer, proof);
  };
  parallel::For(&all_success, kConvLayers.size(), parallel_f);

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
  DenseProve<0>(seed, prove_context, proof0);
  if (!DenseVerify<0>(seed, verify_context, proof0)) return false;

  DenseProof proof1;
  DenseProve<1>(seed, prove_context, proof1);
  if (!DenseVerify<1>(seed, verify_context, proof1)) return false;
  
  return true;
}

inline bool TestAll(std::string const& test_image_path, std::string const& working_path) {
  h256_t seed = misc::RandH256();
  Proof proof;
  if (!Prove(seed, test_image_path, working_path, proof)) return false;
  return Verify(seed, working_path + "/pub", proof);  
}

inline bool Test(std::string const& working_path) {

  bool ret = false;
  //ret = Publish("E:/code/crypto/pod_doc/vgg16_2/features", working_path);
  // ret =TestInfer("E:/code/crypto/pod_doc/vgg16_2/test_image", working_path);
  //ret =TestConv(working_path);
  //ret= TestReluBn(working_path);
  //ret= TestPooling(working_path);
  ret = TestDense(working_path);
  //ret = TestAll("E:/code/crypto/pod_doc/vgg16_2/features", working_path);
  return ret;
}
}  // namespace clink::vgg16