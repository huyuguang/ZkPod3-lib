#pragma once

#include "./adapt.h"
#include "./image_com.h"
#include "./prove.h"
#include "./verify_conv.h"
#include "./verify_dense.h"
#include "./verify_pooling.h"
#include "./verify_relubn.h"

namespace clink::vgg16 {

inline bool Verify(h256_t seed, std::string const& pub_path,
                   dbl::Image const& test_image, Proof const& proof) {
  Tick tick(__FN__);
  VerifyContext context(pub_path);
  AdaptVerifyItemMan item_man;
  ParallelBoolTaskMan task_man;

  std::vector<parallel::BoolTask> tasks;

  // check test image
  tasks.emplace_back([&context, &test_image]() {
    Image image(test_image);
    auto c = pc::PcComputeCommitmentG(image.data, FrZero());
    return c == context.image_com_pub().c[0];
  });

  // conv
  size_t conv_count = kConvLayers.size();
  for (size_t i = 0; i < conv_count; ++i) {
    tasks.emplace_back([&context, &seed, &proof, i, &item_man, &task_man]() {
      return OneConvVerifyPreprocess(seed, context, kConvLayers[i],
                                     proof.conv[i], item_man, task_man);
    });
  }

  // relubn
  tasks.emplace_back([&context, &seed, &proof, &item_man, &task_man]() {
    return ReluBnVerifyPreprocess(seed, context, proof.relubn, item_man,
                                  task_man);
  });

  // pooling
  tasks.emplace_back([&context, &seed, &proof, &item_man, &task_man]() {
    return PoolingVerifyPreprocess(seed, context, proof.pooling, item_man,
                                   task_man);
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
  auto f1 = [&tasks](int64_t i) { return tasks[i](); };
  parallel::For(&all_success, tasks.size(), f1);
  if (!all_success) {
    std::cout << __FN__ << " " << __LINE__<< " Verify failed\n";
    return false;
  }

  std::vector<parallel::BoolTask> bool_tasks;
  task_man.take(bool_tasks);
  bool_tasks.emplace_back([&seed, &item_man,&proof]() {
    return AdaptVerify(seed, item_man, proof.adapt_proof);
  });

  auto f2 = [&bool_tasks](int64_t i) { 
    return bool_tasks[i]();
  };
  parallel::For(&all_success,bool_tasks.size(), f2);

  if (!all_success) {
    std::cout << __FN__ << " " << __LINE__<< " Verify failed\n";
    return false;
  }

  return all_success;
}
}  // namespace clink::vgg16
