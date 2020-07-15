#pragma once

#include "./prove_conv.h"

namespace clink::vgg16 {

inline bool OneConvInputVerifyPreprocess(
    h256_t seed, VerifyContext const& context, size_t layer,
    OneConvProof const& proof, AdaptVerifyItemMan& item_man) {
  Tick tick(__FN__);
 
  struct Ctx {
    std::array<std::vector<Fr>, 9> r;
    std::vector<std::array<int64_t, 9>> B;
    std::vector<Fr> q;
  } ctx;
  
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto DD = D * D;
  auto CDD = C * DD;
  auto KCDD = K * CDD;
  
  // build B
  ctx.B.resize(KCDD);
  for (size_t i = 0; i < CDD; ++i) {
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

  // extend ctx.B
  for (size_t i = CDD; i < KCDD; ++i) {
    ctx.B[i] = ctx.B[i % (CDD)];
  }

  OneConvUpdateSeed(seed, proof.input_pub.cb);
  OneConvComputeInputR(seed, layer, KCDD, ctx.r);

  // build q base ctx.B and ctx.r
  ctx.q.resize(CDD, FrZero());
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < KCDD; ++i) {
      auto const& Bij = ctx.B[i][j];
      if (Bij != -1) {
        ctx.q[Bij] += ctx.r[j][i];
      }
    }
  }

  AdaptVerifyItem adapt_item;
  adapt_item.Init(10, "conv_input_" + std::to_string(layer));
  for (size_t j = 0; j < 9; ++j) {
    adapt_item.a[j] = ctx.r[j];
    adapt_item.cx[j] = proof.input_pub.cb[j];
  }
  adapt_item.a.back() = -ctx.q;
  adapt_item.cx.back() = context.image_com_pub().c[layer];

  item_man.emplace(std::move(adapt_item));
  return true;
}

inline bool OneConvR1csVerifyPreprocess(
    h256_t seed, VerifyContext const& context, size_t layer,
    OneConvProof const& proof,
    ParallelBoolTaskMan& task_man) {
  Tick tick(__FN__);
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto KCDD = K * C * D * D;
  auto order = kLayerTypeOrders[layer].second;

  for (size_t i = 0; i < 9; ++i) {
    if (proof.r1cs_pub.com_w[i] != proof.input_pub.cb[i]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
    if (proof.r1cs_pub.com_w[i + 9] != context.para_com_pub().conv.coef[order][i]) {
      std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
      return false;
    }
  }

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::IpGadget gadget(pb, "vgg16 conv gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  auto r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
  if (proof.r1cs_pub.r1cs_ret_index != r1cs_ret_index) {
    std::cout << __FN__ << ": " << __LINE__ << ": proof invalid\n";
    return false;
  }

  std::shared_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  parallel::BoolTask task = [seed, &proof, r1cs_info, KCDD]() {
    std::vector<std::vector<Fr>> public_w;  // empty
    R1cs::VerifyInput input(KCDD, *r1cs_info, proof.r1cs_pub.com_w, public_w,
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

inline bool OneConvOutputVerifyPreprocess(
    h256_t seed, VerifyContext const& context, size_t layer,
    OneConvProof const& proof, AdaptVerifyItemMan& item_man) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto order = kLayerTypeOrders[layer].second;

  G1 const& cx = proof.r1cs_pub.com_w[proof.r1cs_pub.r1cs_ret_index];
  G1 const& cz = context.image_com_pub().c[layer + 1];
  G1 cb =
      context.para_com_pub().conv.bias[order] * fp::RationalConst<8, 24>().kFrN;

  if (cz != proof.output_pub.cy + cb) {
    std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
    return false;
  }

  OneConvUpdateSeed(seed, proof.output_pub.cy);

  // r.size = KDD
  std::vector<Fr> r = OneConvComputeOutputR(seed, layer, K * D * D);

  // s.size = KCDD
  std::vector<Fr> s = OneConvOutputR2S(K, C, D, r);

  AdaptVerifyItem adapt_item;
  adapt_item.Init(2, "conv_output_" + std::to_string(layer));
  adapt_item.a[0] = std::move(s);
  adapt_item.cx[0] = cx;
  adapt_item.a[1] = std::move(r);
  adapt_item.a[1] = -adapt_item.a[1];
  adapt_item.cx[1] = proof.output_pub.cy;
  item_man.emplace(std::move(adapt_item));

  return true;

  //// <x,s>==<y,r>
  //EqualIp<HyraxA>::VerifyInput input(s, cx, pc::kGetRefG, r, proof.output.cy,
  //                                   pc::kGetRefG);
  //if (!EqualIp<HyraxA>::Verify(seed, proof.output.eip_proof, input)) {
  //  std::cout << __FN__ << ": " << __LINE__ << ": verify failed\n";
  //  return false;
  //}

  //return true;
}

inline bool OneConvVerifyPreprocess(
    h256_t seed, VerifyContext const& context, size_t layer,
    OneConvProof const& proof, AdaptVerifyItemMan& item_man,
    ParallelBoolTaskMan& task_man) {
  Tick tick(__FN__);

  // can parallel but need to protect adapt_items and parallel_tasks
  if (!OneConvInputVerifyPreprocess(seed, context, layer, proof, item_man)) {
#ifdef _DEBUG_CHECK
    throw std::runtime_error("oops");
#endif
    return false;
  }

  if (!OneConvR1csVerifyPreprocess(seed, context, layer, proof, task_man)) {
#ifdef _DEBUG_CHECK
    throw std::runtime_error("oops");
#endif
    return false;
  }

  if (!OneConvOutputVerifyPreprocess(seed, context, layer, proof, item_man)) {
#ifdef _DEBUG_CHECK
    throw std::runtime_error("oops");
#endif
    return false;
  }

  return true;
}
}  // namespace clink::vgg16
