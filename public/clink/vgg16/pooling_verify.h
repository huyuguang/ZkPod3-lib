#pragma once

#include "./pooling_prove.h"

namespace clink::vgg16 {
inline bool PoolingInputVerifyPreprocess(h256_t seed,
                                         VerifyContext const& context,
                                         PoolingProof const& proof,
                                         SafeVec<AdaptVerifyItem>& adapt_man) {
  Tick tick(__FN__);
  auto const& input_pub = proof.input_pub;

  for (size_t l = 0; l < kPoolingLayers.size(); ++l) {
    auto layer = kPoolingLayers[l];
    if (context.image_com_pub().c[layer] != input_pub.cx[l]) {
      std::cout << Tick::GetIndentString() << ": proof invalid\n";
      return false;
    }
  }

  PoolingUpdateSeed(seed, input_pub.cx[5], input_pub.cx[6], input_pub.cx[7],
                    input_pub.cx[8]);

  std::array<std::vector<Fr>, 9> q;
  PoolingBuildInputQ(seed, q);

  AdaptVerifyItem adapt_item;
  adapt_item.Init(9, PoolingAdaptTag(true), FrZero());
  for (size_t j = 0; j < 9; ++j) {
    adapt_item.a[j] = std::move(q[j]);
    adapt_item.cx[j] = input_pub.cx[j];
    if (j < 5) adapt_item.a[j] = -adapt_item.a[j];
  }
  adapt_man.emplace(std::move(adapt_item));
  return true;
}

inline bool PoolingR1csVerifyPreprocess(h256_t seed,
                                        VerifyContext const& /*context*/,
                                        PoolingProof const& proof,
                                        SafeVec<R1csVerifyItem>& r1cs_man) {
  Tick tick(__FN__);
  (void)seed;
  if (proof.r1cs_pub.com_w[0] != proof.input_pub.cx[5]) {  // a
    std::cout << Tick::GetIndentString() << ": proof invalid 1\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[1] != proof.input_pub.cx[6]) {  // b
    std::cout << Tick::GetIndentString() << ": proof invalid 2\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[2] != proof.input_pub.cx[7]) {  // c
    std::cout << Tick::GetIndentString() << ": proof invalid 3\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[3] != proof.input_pub.cx[8]) {  // d
    std::cout << Tick::GetIndentString() << ": proof invalid 4\n";
    return false;
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::PoolingGadget<8, 24> gadget(pb, "vgg16 pooling gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  // see protoboard<FieldT>::val
  auto r1cs_ret_index = gadget.ret().index - 1;
  if (proof.r1cs_pub.r1cs_ret_index != r1cs_ret_index) {
    std::cout << Tick::GetIndentString() << ": proof invalid 5\n";
    return false;
  }

  auto n = PoolingGetCircuitCount();
  R1csVerifyItem item;
  item.public_w.reset(new std::vector<std::vector<Fr>>);
  item.r1cs_info.reset(new R1csInfo(pb));
  item.r1cs_input.reset(new R1cs::VerifyInput(
      n, *item.r1cs_info, PoolingR1csTag(), proof.r1cs_pub.com_w,
      *item.public_w, pc::kGetRefG1));

  r1cs_man.emplace(std::move(item));

  return true;
}

inline bool PoolingOutputVerifyPreprocess(h256_t seed,
                                          VerifyContext const& context,
                                          PoolingProof const& proof,
                                          SafeVec<AdaptVerifyItem>& adapt_man) {
  Tick tick(__FN__);
  auto const& output_pub = proof.output_pub;

  for (size_t l = 0; l < kPoolingLayers.size(); ++l) {
    auto layer = kPoolingLayers[l] + 1;
    if (context.image_com_pub().c[layer] != output_pub.cx[l + 1]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  if (proof.r1cs_pub.com_w[proof.r1cs_pub.r1cs_ret_index] != output_pub.cx[0]) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  PoolingUpdateSeed(seed, proof.r1cs_pub);

  std::array<std::vector<Fr>, 6> q;
  PoolingBuildOutputQ(seed, q);

  AdaptVerifyItem adapt_item;
  adapt_item.Init(6, PoolingAdaptTag(false), FrZero());
  for (size_t i = 0; i < 6; ++i) {
    adapt_item.a[i] = std::move(q[i]);
    adapt_item.cx[i] = output_pub.cx[i];
    if (i == 0) {
      adapt_item.a[i] = -adapt_item.a[i];
    }
  }
  adapt_man.emplace(std::move(adapt_item));
  return true;
}

inline bool PoolingVerifyPreprocess(h256_t seed, VerifyContext const& context,
                                    PoolingProof const& proof,
                                    SafeVec<AdaptVerifyItem>& adapt_man,
                                    SafeVec<R1csVerifyItem>& r1cs_man) {
  Tick tick(__FN__);

  std::array<parallel::VoidTask, 3> tasks;

  tasks[0] = [&seed, &context, &proof, &adapt_man]() {
    CHECK(PoolingInputVerifyPreprocess(seed, context, proof, adapt_man), "");
  };

  tasks[1] = [&seed, &context, &proof, &r1cs_man]() {
    CHECK(PoolingR1csVerifyPreprocess(seed, context, proof, r1cs_man), "");
  };

  tasks[2] = [&seed, &context, &proof, &adapt_man]() {
    CHECK(PoolingOutputVerifyPreprocess(seed, context, proof, adapt_man), "");
  };

  parallel::Invoke(tasks);

  return true;
}
}  // namespace clink::vgg16