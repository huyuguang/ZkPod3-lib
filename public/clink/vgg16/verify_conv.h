#pragma once

#include "./prove_conv.h"

namespace clink::vgg16 {

inline bool OneConvInputVerify(h256_t seed, VerifyContext const& context,
                               size_t layer, OneConvProof const& proof) {
  Tick tick(__FN__);
  size_t const order = kLayerTypeOrders[layer].second;
  // auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  G1 const& cx = context.image_com_pub().c[layer];

  struct Ctx {
    std::array<std::vector<Fr>, 9> r;
    std::vector<std::array<int64_t, 9>> B;
    std::vector<Fr> q;
  } ctx;

  auto get_col_b_u = [&context, order](int64_t i) -> G1 const& {
    auto range = context.auxi().data_u_conv(order);
    if (i >= (range.second - range.first)) throw std::runtime_error("oops");
    return range.first[i];
  };

  OneConvUpdateSeed(seed, proof.input.cb);
  OneConvComputeInputR(seed, layer, C * D * D, ctx.r);

  // build B
  ctx.B.resize(C * D * D);
  for (size_t i = 0; i < C * D * D; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      size_t m = j / 3;
      size_t n = j % 3;
      size_t r = i / (D * D);
      size_t p = i % (D * D);
      size_t q = p / D;
      size_t o = p % D;
      size_t ii = q + m;
      size_t jj = o + n;
      if (ii == 0 || jj == 0 || ii == (D + 1) || jj == (D + 1)) {
        ctx.B[i][j] = -1;
      } else {
        ctx.B[i][j] = r * D * D + (ii - 1) * D + (jj - 1);
      }
    }
  }

  // build q base B and r
  ctx.q.resize(C * D * D, FrZero());
  //std::fill(ctx.q.begin(), ctx.q.end(), FrZero());
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < C * D * D; ++i) {
      auto const& Bij = ctx.B[i][j];
      if (Bij != -1) {
        ctx.q[Bij] += ctx.r[j][i];
      }
    }
  }

  bool all_success = false;
  auto parallel_f = [&proof, &ctx, &seed, &cx, &get_col_b_u](int64_t j) {
    if (j < 9) {
      HyraxA::CommitmentPub com_pub(proof.input.cb[j], proof.input.a3_cy[j]);
      HyraxA::VerifyInput input(ctx.r[j], com_pub, get_col_b_u, pc::PcG(0));
      return HyraxA::Verify(proof.input.rb_proofs[j], seed, input);
    } else {
      G1 a3_cz = std::accumulate(proof.input.a3_cy.begin(),
                                 proof.input.a3_cy.end(), G1Zero());
      HyraxA::CommitmentPub com_pub(cx, a3_cz);
      HyraxA::VerifyInput input(ctx.q, com_pub, pc::kGetRefG, pc::PcG(0));
      return HyraxA::Verify(proof.input.xq_proof, seed, input);
    }
  };
  parallel::For(&all_success, 10, parallel_f);

  if (!all_success) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  return true;
}

inline bool OneConvR1csVerify(h256_t seed, VerifyContext const& context,
                              size_t layer, OneConvProof const& proof) {
  Tick tick(__FN__);
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto order = kLayerTypeOrders[layer].second;

  for (size_t i = 0; i < 9; ++i) {
    if (proof.r1cs.com_w[i] != proof.input.cb[i]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
    if (proof.r1cs.com_w[i + 9] != context.para_com_pub().conv.coef[order][i]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::IpGadget gadget(pb, "vgg16 conv gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  auto r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
  if (proof.r1cs.r1cs_ret_index != r1cs_ret_index) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  // auto s = r1cs_info->num_variables;
  auto n = K * C * D * D;
  std::vector<std::vector<Fr>> public_w;  // empty

  R1cs::VerifyInput input(n, *r1cs_info, proof.r1cs.com_w, public_w,
                          pc::kGetRefG);
  if (!R1cs::Verify(proof.r1cs.r1cs_proof, seed, input)) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  return true;
}
inline bool OneConvOutputVerify(h256_t seed, VerifyContext const& context,
                                size_t layer, OneConvProof const& proof) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto order = kLayerTypeOrders[layer].second;

  G1 const& cx = proof.r1cs.com_w[proof.r1cs.r1cs_ret_index];
  G1 const& cz = context.image_com_pub().c[layer + 1];
  G1 cb =
      context.para_com_pub().conv.bias[order] * fp::RationalConst<8, 24>().kFrN;

  if (cz != proof.output.cy + cb) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  OneConvUpdateSeed(seed, proof.output.cy);

  // r.size = KDD
  std::vector<Fr> r = OneConvComputeOutputR(seed, layer, K * D * D);

  // s.size = KCDD
  std::vector<Fr> s = OneConvOutputR2S(K, C, D, r);

  // <x,s>==<y,r>
  EqualIp<HyraxA>::VerifyInput input(s, cx, pc::kGetRefG, r, proof.output.cy,
                                     pc::kGetRefG);
  if (!EqualIp<HyraxA>::Verify(seed, proof.output.eip_proof, input)) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  return true;
}

inline bool OneConvVerify(h256_t seed, VerifyContext const& context,
                          size_t layer, OneConvProof const& proof) {
  Tick tick(__FN__);
  // if (!OneConvInputVerify(seed, context, layer, proof)) return false;
  // if (!OneConvR1csVerify(seed, context, layer, proof)) return false;
  // if (!OneConvOutputVerify(seed, context, layer, proof)) return false;
  // std::cout << "OneConvVerify " << layer << " success\n";
  // return true;

  std::array<std::atomic<bool>, 3> rets;
  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&seed, &context, layer, &proof, &rets]() {
    rets[0] = OneConvInputVerify(seed, context, layer, proof);
  };
  tasks[1] = [&seed, &context, layer, &proof, &rets]() {
    rets[1] = OneConvR1csVerify(seed, context, layer, proof);
  };
  tasks[2] = [&seed, &context, layer, &proof, &rets]() {
    rets[2] = OneConvOutputVerify(seed, context, layer, proof);
  };

  parallel::Invoke(tasks);

  if (!rets[0] || !rets[1] || !rets[2]) return false;

  std::cout << "OneConvVerify " << layer << " success\n";
  return true;
}
}  // namespace clink::vgg16