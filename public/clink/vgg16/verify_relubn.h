#pragma once

#include "./prove_relubn.h"

namespace clink::vgg16 {

inline bool ReluBnInOutVerifyPreprocess(
    h256_t seed, VerifyContext const& context, ReluBnProof const& proof,
    AdaptVerifyItemMan& item_man) {
  auto const& io_pub = proof.io_pub;

  if (io_pub.cx.size() != kReluBnLayers.size() * 2 + 2) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    if (context.image_com_pub().c[layer] != io_pub.cx[2 + order * 2]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
    if (context.image_com_pub().c[layer + 1] !=
        io_pub.cx[2 + order * 2 + 1]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  std::vector<std::vector<Fr>> a(io_pub.cx.size());
  size_t total_size = 0;
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    auto const& info = kImageInfos[layer];
    auto size = info.C * info.D * info.D;
    a[2 + order * 2].resize(size);
    a[2 + order * 2 + 1].resize(size);
    total_size += size;
  }
  a[0].resize(total_size);
  a[1].resize(total_size);

  ReluBnUpdateSeed(seed, io_pub.cx[0], io_pub.cx[1]);
  a[0] = ReluBnComputeFst(seed, "relubn input", total_size);
  a[1] = ReluBnComputeFst(seed, "relubn output", total_size);

  auto it_i = a[0].begin();
  auto it_o = a[1].begin();
  for (size_t i = 1; i < a.size() / 2; ++i) {
    auto& a_i = a[2 * i];
    auto& a_o = a[2 * i + 1];
    auto size = a_i.size();
    std::copy(it_i, it_i + size, a_i.begin());
    std::copy(it_o, it_o + size, a_o.begin());
    it_i += size;
    it_o += size;
  }

  AdaptVerifyItem adapt_item_in;
  AdaptVerifyItem adapt_item_out;
  adapt_item_in.Init(io_pub.cx.size() / 2, "relubn verify in");
  adapt_item_out.Init(io_pub.cx.size() / 2, "relubn verify out");
  for (size_t j = 0; j < io_pub.cx.size() / 2; ++j) {
    adapt_item_in.a[j] = std::move(a[j * 2]);
    adapt_item_in.cx[j] = io_pub.cx[j * 2];
    if (j == 0) adapt_item_in.a[j] = -adapt_item_in.a[j];

    adapt_item_out.a[j] = std::move(a[j * 2 + 1]);
    adapt_item_out.cx[j] = io_pub.cx[j * 2 + 1];
    if (j == 0) adapt_item_out.a[j] = -adapt_item_out.a[j];
  }
  item_man.emplace(std::move(adapt_item_in));
  item_man.emplace(std::move(adapt_item_out));
  return true;
}

inline bool ReluBnR1csVerifyPreprocess(
    h256_t seed, VerifyContext const& context, ReluBnProof const& proof,
    ParallelBoolTaskMan& task_man) {
  Tick tick(__FN__);
  if (proof.r1cs_pub.com_w[0] != proof.io_pub.cx[0]) {  // in
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[1] != proof.io_pub.cx[1]) {  // out
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[2] != context.para_com_pub().bn.alpha) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[3] != context.para_com_pub().bn.beta) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs_pub.com_w[4] != context.para_com_pub().bn.mu) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::ReluBnGadget<8, 24> gadget(pb, "vgg16 relubn gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  std::shared_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  auto n = ReluBnGetCircuitCount();

  parallel::BoolTask task = [seed, &proof, r1cs_info, n]() {
    std::vector<std::vector<Fr>> public_w;  // empty
    R1cs::VerifyInput input(n, *r1cs_info, proof.r1cs_pub.com_w, public_w,
                            pc::kGetRefG);
    if (!R1cs::Verify(proof.r1cs_proof, seed, input)) {
      std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
      return false;
    }
    return true;
  };
  task_man.emplace(std::move(task));
  return true;
}

inline bool ReluBnVerifyPreprocess(h256_t seed, VerifyContext const& context,
                                   ReluBnProof const& proof,
                                   AdaptVerifyItemMan& item_man,
                                   ParallelBoolTaskMan& task_man) {
  Tick tick(__FN__);

  // can parallel but need to protect adapt_items and parallel_tasks
  if (!ReluBnInOutVerifyPreprocess(seed, context, proof, item_man)) {
#ifdef _DEBUG_CHECK
    throw std::runtime_error("oops");
#endif
    return false;
  }

  if (!ReluBnR1csVerifyPreprocess(seed, context, proof, task_man)) {
#ifdef _DEBUG_CHECK
    throw std::runtime_error("oops");
#endif
    return false;
  }

  return true;
}
}  // namespace clink::vgg16