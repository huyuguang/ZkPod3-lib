#pragma once

#include "./image_com.h"
#include "./prove.h"
#include "./verify_conv.h"
#include "./verify_dense.h"
#include "./verify_pooling.h"
#include "./verify_relubn.h"

namespace clink::vgg16 {
//
//inline bool CheckProofFormat(VerifyContext context, Proof const proof) {
//
//}

inline bool Verify(h256_t seed, std::string const& pub_path,
                   Proof const& proof) {
  Tick tick(__FN__);
  VerifyContext context(pub_path);

  std::vector<std::function<bool()>> tasks;

  // conv
  for (size_t i = 0; i < kConvLayers.size(); ++i) {
    tasks.emplace_back([&context, &seed, &proof, i]() {
      return OneConvVerify(seed, context, kConvLayers[i], proof.conv[i]);
    });
  }

  // relubn
  tasks.emplace_back([&context, &seed, &proof]() {
    return ReluBnVerify(seed, context, proof.relubn);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof]() {
    return PoolingVerify(seed, context, proof.pooling);
  });

  // dense0
  tasks.emplace_back([&context, &seed, &proof]() {
    return DenseVerify<0>(seed, context, proof.dense0);
  });

  // dense1
  tasks.emplace_back([&context, &seed, &proof]() {
    return DenseVerify<1>(seed, context, proof.dense1);
  });

  bool all_success = false;
  auto parallel_f = [&tasks](int64_t i) { return tasks[i](); };
  parallel::For(&all_success, tasks.size(), parallel_f);

  if (all_success) {
    std::cout << "Verify success\n";
  } else {
    std::cout << "Verify failed\n";
  }
  return all_success;
}
}  // namespace clink::vgg16