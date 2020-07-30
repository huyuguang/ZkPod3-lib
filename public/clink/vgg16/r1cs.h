#pragma once

#include <memory>

#include "./policy.h"
#include "./r1cs_pub.h"
#include "./safevec.h"
#include "clink/batch_r1cs.h"

namespace clink::vgg16 {
struct R1csProveItem {
  using BaseR1cs = ParallelR1cs<Policy>;
  using ProveInput = typename BaseR1cs::ProveInput;
  std::shared_ptr<BaseR1csSec> r1cs_sec;
  std::shared_ptr<ProveInput> r1cs_input;
};

struct R1csVerifyItem {
  using BaseR1cs = ParallelR1cs<Policy>;
  using VerifyInput = typename BaseR1cs::VerifyInput;
  std::shared_ptr<R1csInfo> r1cs_info;
  std::shared_ptr<std::vector<std::vector<Fr>>> public_w;
  std::shared_ptr<VerifyInput> r1cs_input;
};

using R1csProveItemMan = SafeVec<R1csProveItem>;
using R1csVerifyItemMan = SafeVec<R1csVerifyItem>;

inline void R1csProve(h256_t seed, R1csProveItemMan& item_man,
                      clink::ParallelR1cs<R1cs>::Proof& proof) {
  Tick tick(__FN__);
  std::vector<R1csProveItem> items;
  item_man.take(items);

  std::vector<BatchR1cs<Policy>::ProveInput*> inputs(items.size());
  for (size_t i = 0; i < inputs.size(); ++i) {
    inputs[i] = items[i].r1cs_input.get();
  }
  BatchR1cs<Policy>::Prove(proof, seed, std::move(inputs));
}

inline bool R1csVerify(h256_t seed, R1csVerifyItemMan& item_man,
                       clink::ParallelR1cs<R1cs>::Proof const& proof) {
  Tick tick(__FN__);
  std::vector<R1csVerifyItem> items;
  item_man.take(items);

  std::vector<BatchR1cs<Policy>::VerifyInput*> inputs(items.size());
  for (size_t i = 0; i < inputs.size(); ++i) {
    inputs[i] = items[i].r1cs_input.get();
    std::cout << inputs[i]->unique_tag
              << ", public_w.size(): " << inputs[i]->public_w.size() << "\n";
  }

  return BatchR1cs<Policy>::Verify(proof, seed, std::move(inputs));
}
}  // namespace clink::vgg16