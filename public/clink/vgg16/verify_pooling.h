#pragma once

#include "./prove_pooling.h"

namespace clink::vgg16 {
inline bool PoolingInputVerify(h256_t seed, VerifyContext const& context,
                               PoolingProof const& proof) {
  auto const& in_proof = proof.input;
  if (in_proof.ip_com_pubs.size() != in_proof.ip_proofs.size()) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (in_proof.ip_com_pubs.size() != 5 + 4) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  auto layers = PoolingGetLayers();
  for (size_t l = 0; l < layers.size(); ++l) {
    auto layer = layers[l];
    if (context.image_com_pub().c[layer] != in_proof.ip_com_pubs[l].xi) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  G1 sum_com_5 = G1Zero();
  for (size_t i = 0; i < 5;i++) {
    sum_com_5 += in_proof.ip_com_pubs[i].tau;
  }
  G1 sum_com_4 = G1Zero();
  for (size_t i = 5; i < 9;i++) {
    sum_com_4 += in_proof.ip_com_pubs[i].tau;
  }
  if (sum_com_5 != sum_com_4) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  PoolingUpdateSeed(seed, in_proof.ip_com_pubs[5].xi,
                    in_proof.ip_com_pubs[6].xi, in_proof.ip_com_pubs[7].xi,
                    in_proof.ip_com_pubs[8].xi);

  std::array<std::vector<Fr>, 9> q;
  PoolingBuildInputQ(seed, q);

  bool all_success = false;
  auto parallel_f = [&seed, &q, &in_proof](int64_t i) {
    HyraxA::VerifyInput input(q[i], in_proof.ip_com_pubs[i], pc::kGetRefG,
                              pc::PcG(0));
    return HyraxA::Verify(in_proof.ip_proofs[i], seed, input);
  };
  parallel::For(&all_success, q.size(), parallel_f);
  if (!all_success) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }
  return all_success;
}

inline bool PoolingR1csVerify(h256_t seed, VerifyContext const& /*context*/,
                              PoolingProof const& proof) {
  Tick tick(__FN__);
  if (proof.r1cs.com_w[0] != proof.input.ip_com_pubs[5].xi) { // a
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }
  
  if (proof.r1cs.com_w[1] != proof.input.ip_com_pubs[6].xi) { // b
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }
  
  if (proof.r1cs.com_w[2] != proof.input.ip_com_pubs[7].xi) { // c
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (proof.r1cs.com_w[3] != proof.input.ip_com_pubs[8].xi) { // d
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::PoolingGadget<8,24> gadget(pb, "vgg16 pooling gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  auto r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
  if (proof.r1cs.r1cs_ret_index != r1cs_ret_index) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  // auto s = r1cs_info->num_variables;
  auto n = PoolingGetCircuitCount();
  std::vector<std::vector<Fr>> public_w;  // empty

  R1cs::VerifyInput input(n, *r1cs_info, proof.r1cs.com_w, public_w,
                          pc::kGetRefG);
  if (!R1cs::Verify(proof.r1cs.r1cs_proof, seed, input)) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  return true;
}

inline bool PoolingOutputVerify(h256_t seed, VerifyContext const& context,
                                PoolingProof const& proof) {
  auto const& out_proof = proof.output;
  if (out_proof.ip_com_pubs.size() != out_proof.ip_proofs.size()) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  if (out_proof.ip_com_pubs.size() != 5 + 1) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  auto layers = PoolingGetLayers();
  for (size_t l = 0; l < layers.size(); ++l) {
    auto layer = layers[l] + 1;
    if (context.image_com_pub().c[layer] != out_proof.ip_com_pubs[l + 1].xi) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }
  
  if (proof.r1cs.com_w[proof.r1cs.r1cs_ret_index] !=
      out_proof.ip_com_pubs[0].xi) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  G1 sum_com = G1Zero();
  for (size_t i = 1; i < 6; i++) {
    sum_com += out_proof.ip_com_pubs[i].tau;
  }

  if (sum_com != out_proof.ip_com_pubs[0].tau) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  PoolingUpdateSeed(seed, proof.r1cs);

  std::array<std::vector<Fr>, 6> q;
  PoolingBuildOutputQ(seed, q);

  bool all_success = false;
  auto parallel_f = [&seed, &q, &out_proof](int64_t i) {
    HyraxA::VerifyInput input(q[i], out_proof.ip_com_pubs[i], pc::kGetRefG,
                              pc::PcG(0));
    return HyraxA::Verify(out_proof.ip_proofs[i], seed, input);
  };
  parallel::For(&all_success, q.size(), parallel_f);
  if (!all_success) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }
  return all_success;
}

inline bool PoolingVerify(h256_t seed, VerifyContext const& context,
                          PoolingProof const& proof) {
  Tick tick(__FN__);

  std::array<std::atomic<bool>, 3> rets;
  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&seed, &context, &proof, &rets]() {
    rets[0] = PoolingInputVerify(seed, context, proof);
  };
  tasks[1] = [&seed, &context, &proof, &rets]() {
    rets[1] = PoolingR1csVerify(seed, context, proof);
  };
  tasks[2] = [&seed, &context, &proof, &rets]() {
    rets[2] = PoolingOutputVerify(seed, context, proof);
  };

  parallel::Invoke(tasks);

  if (!rets[0] || !rets[1] || !rets[2]) return false;

  std::cout << "PoolingVerify success\n";
  return true;
}
}  // namespace clink::vgg16