#pragma once

#include "./prove_relubn.h"

namespace clink::vgg16 {

inline bool ReluBnInOutVerify(h256_t seed, VerifyContext const& context,
                              ReluBnProof const& proof) {
  auto const& io_proof = proof.inout;
  if (io_proof.ip_com_pubs.size() != io_proof.ip_proofs.size()) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (io_proof.ip_proofs.size() != kBnLayerInfos.size() * 2 + 2) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  for (size_t l = 0; l < kLayerTypeOrders.size(); ++l) {
    auto layer_type = kLayerTypeOrders[l].first;
    auto layer_order = kLayerTypeOrders[l].second;
    if (layer_type != kReluBn) continue;
    if (context.image_com_pub().c[l] !=
        io_proof.ip_com_pubs[2 + layer_order * 2].xi) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
    if (context.image_com_pub().c[l + 1] !=
        io_proof.ip_com_pubs[2 + layer_order * 2 + 1].xi) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  G1 chk_in_com_xa = G1Zero();
  G1 chk_out_com_xa = G1Zero();
  for (size_t i = 2; i < io_proof.ip_com_pubs.size(); ++i) {
    auto const& pub = io_proof.ip_com_pubs[i];
    if (i%2 == 0) {
      chk_in_com_xa += pub.tau;
    } else {
      chk_out_com_xa += pub.tau;
    }
  }
  if (chk_in_com_xa != io_proof.ip_com_pubs[0].tau) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }
  if (chk_out_com_xa != io_proof.ip_com_pubs[1].tau) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  std::vector<std::vector<Fr>> a(io_proof.ip_proofs.size());
  size_t total_size = 0;
  for (size_t l = 0; l < kLayerTypeOrders.size(); ++l) {
    auto layer_type = kLayerTypeOrders[l].first;
    auto layer_order = kLayerTypeOrders[l].second;
    if (layer_type != kReluBn) continue;
    auto const& info = kImageInfos[l];
    auto size = info.channel_count * info.dimension * info.dimension;
    a[2 + layer_order * 2].resize(size);
    a[2 + layer_order * 2 + 1].resize(size);
    total_size += size;
  }  
  a[0].resize(total_size);
  a[1].resize(total_size);

  ReluBnUpdateSeed(seed, io_proof.ip_com_pubs[0].xi,
                   io_proof.ip_com_pubs[1].xi);
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

  bool all_success = false;
  auto parallel_f = [&seed, &a, &io_proof](int64_t i) {
    HyraxA::VerifyInput input(a[i], io_proof.ip_com_pubs[i], pc::kGetRefG,
                             pc::PcG(0));
    return HyraxA::Verify(io_proof.ip_proofs[i],seed, input);
  };
  parallel::For(&all_success, a.size(), parallel_f);
  if (!all_success) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }
  return all_success;
}

inline bool ReluBnR1csVerify(h256_t seed, VerifyContext const& context,
                             ReluBnProof const& proof) {
  Tick tick(__FN__);
  if (proof.r1cs.com_w[0] != proof.inout.ip_com_pubs[0].xi) {  // in
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs.com_w[1] != proof.inout.ip_com_pubs[1].xi) {  // out
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs.com_w[2] != context.para_com_pub().bn.alpha) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs.com_w[3] != context.para_com_pub().bn.beta) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs.com_w[4] != context.para_com_pub().bn.mu) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::ReluBnGadget<8,24> gadget(pb, "vgg16 relubn gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  // auto s = r1cs_info->num_variables;
  auto n = ReluBnGetCircuitCount();
  std::vector<std::vector<Fr>> public_w;  // empty

  R1cs::VerifyInput input(n, *r1cs_info, proof.r1cs.com_w, public_w,
                          pc::kGetRefG);
  if (!R1cs::Verify(proof.r1cs.r1cs_proof, seed, input)) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  return true;
}

inline bool ReluBnVerify(h256_t seed, VerifyContext const& context,
                         ReluBnProof const& proof) {
  Tick tick(__FN__);

  std::array<std::atomic<bool>, 2> rets;
  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&seed, &context, &proof, &rets]() {
    rets[0] = ReluBnInOutVerify(seed, context, proof);
  };
  tasks[1] = [&seed, &context, &proof, &rets]() {
    rets[1] = ReluBnR1csVerify(seed, context, proof);
  };

  parallel::Invoke(tasks);

  if (!rets[0] || !rets[1]) return false;

  std::cout << "ReluBnVerify success\n";
  return true;
}
}  // namespace clink::vgg16