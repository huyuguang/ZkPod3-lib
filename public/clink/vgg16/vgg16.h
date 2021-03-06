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

inline bool TestPublish(std::string const& para_path,
                    std::string const& working_path) {
  return Publish(para_path, working_path);
}

inline bool TestInfer(std::string const& test_image_path,
               std::string const& working_path) {
  dbl::Image test_image(kImageInfos[0]);
  if (!dbl::LoadTestImage(test_image_path, test_image)) {
    return false;
  }

  return InferAndCommit(test_image, working_path);
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

inline bool TestSerialize(std::string const& working_path) {
  try {
    Proof proof(working_path + "/proof");
    std::cout << "proof: " << YasGetBinLen(proof) << "\n";
    for (size_t i = 0; i < proof.conv.size(); ++i) {
      std::cout << "proof conv " << i << ": " << YasGetBinLen(proof.conv[i])
                << "\n";
      std::cout << "\t input pub:" << YasGetBinLen(proof.conv[i].input_pub) << "\n";
      std::cout << "\t r1cs pub:" << YasGetBinLen(proof.conv[i].r1cs_pub) << "\n";
      std::cout << "\t output pub:" << YasGetBinLen(proof.conv[i].output_pub) << "\n";
    }
    std::cout << "\n";
    std::cout << "proof relubn: " << YasGetBinLen(proof.relubn) << "\n";
    std::cout << "\t io pub: " << YasGetBinLen(proof.relubn.io_pub) << "\n";
    std::cout << "\t r1cs pub: " << YasGetBinLen(proof.relubn.r1cs_pub) << "\n";
    std::cout << "\n";

    std::cout << "proof pooling: " << YasGetBinLen(proof.pooling) << "\n";
    std::cout << "\t input pub: " << YasGetBinLen(proof.pooling.input_pub) << "\n";
    std::cout << "\t r1cs pub: " << YasGetBinLen(proof.pooling.r1cs_pub) << "\n";
    std::cout << "\t output pub: " << YasGetBinLen(proof.pooling.output_pub) << "\n";
    std::cout << "\n";

    std::cout << "proof dense0: " << YasGetBinLen(proof.dense0) << "\n";
    std::cout << "proof dense1: " << YasGetBinLen(proof.dense1) << "\n";
    std::cout << "\n";

    std::cout << "proof adapt: " << YasGetBinLen(proof.adapt_proof) << "\n";
    std::cout << "proof r1cs conv: " << YasGetBinLen(proof.r1cs_proof) << "\n";
    std::cout << "\n";
    return true;
  } catch (std::exception& e) {
    std::cerr << __FN__ << ": " << __LINE__ << " " << e.what() << "\n";
    return false;
  }
}

inline bool TestProve(std::string const& test_image_path,
                    std::string const& working_path) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
   
  dbl::Image test_image(kImageInfos[0]);
  if (!dbl::LoadTestImage(test_image_path, test_image)) {
    return false;
  }

  Proof proof;
  if (!Prove(seed, test_image, working_path, proof)) return false;
  YasSaveBin(working_path + "/proof", proof);
  if (!TestSerialize(working_path)) return false;
  return Verify(seed, working_path + "/pub", test_image, proof);
}
//
//// "E:/code/crypto/pod_doc/vgg16_2/test_image"
//inline bool TestProve(std::string const& test_image_path,
//                      std::string const& working_path) {
//  Tick tick(__FN__);
//
//  bool ret = false;
//
//  //ret = Publish("E:/code/crypto/pod_doc/vgg16_2/features", working_path);
//  // ret =TestInfer("E:/code/crypto/pod_doc/vgg16_2/test_image", working_path);
//  //ret =TestConv(working_path);
//  // ret= TestReluBn(working_path);
//  //ret= TestPooling(working_path);
//  //ret = TestDense(working_path);
//  ret = TestAll(test_image_path, working_path);
//  //ret = TestSerialize(working_path);
//  return ret;
//}

inline bool Test() {
#ifdef _WIN32
  std::string test_image_path = "../../../../data/vgg16/test_image";
  std::string working_path = "../../../../data/vgg16/working";
  std::string features_path = "../../../../data/vgg16/features";
#else
  std::string test_image_path = "../../data/vgg16/test_image";
  std::string working_path = "../../data/vgg16/working";
  std::string features_path = "../../data/vgg16/features";
#endif

  if (!fs::is_directory(working_path) &&
      !fs::create_directories(working_path)) {
    std::cerr << "Create " << working_path << " failed\n";
    return false;
  }

  boost::system::error_code ec;
  if (!fs::is_directory(working_path + "/pub", ec) ||
      !fs::is_directory(working_path + "/sec", ec)) {
    Publish(features_path, working_path);
  }

  //TestConv(working_path); // TODO
  return TestProve(test_image_path, working_path);
}

}  // namespace clink::vgg16