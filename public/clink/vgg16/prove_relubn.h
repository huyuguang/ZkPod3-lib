#pragma once

#include "./context.h"
#include "./image_com.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

inline void ReluBnUpdateSeed(h256_t& seed, G1 const& x, G1 const& y) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, x);
  HashUpdate(hash, y);
  hash.Final(seed.data());
}

inline std::vector<Fr> ReluBnComputeFst(h256_t const& seed,
                                        std::string const& salt, size_t size) {
  std::vector<Fr> r(size);
  ComputeFst(seed, salt + std::to_string(size), r);
  return r;
}

enum ReluBnImageType { kCombinedInput, kCombinedOutput, kInput, kOutput };

struct ReluBnImage {
  ReluBnImageType type;
  std::vector<Fr> x;
  G1 com_x;
  Fr com_x_r;

  std::vector<Fr> a;
  Fr xa;
  Fr com_xa_r;
  G1 com_xa;
};

struct ReluBnInOutPub {
  std::vector<HyraxA::CommitmentPub> ip_com_pubs;
  std::vector<HyraxA::Proof> ip_proofs;
};

struct ReluBnInOutSec {
  std::vector<Fr> input;   // combined inputs
  std::vector<Fr> output;  // combined outputs
  Fr com_in_r;
  Fr com_out_r;
};

inline void ReluBnBuildChallenge(h256_t const& seed,
                                 std::vector<ReluBnImage>& images) {
  Tick tick(__FN__);
  assert(images[0].type == ReluBnImageType::kCombinedInput);
  assert(images[1].type == ReluBnImageType::kCombinedOutput);

  images[0].a = ReluBnComputeFst(seed, "relubn input", images[0].x.size());
  images[1].a = ReluBnComputeFst(seed, "relubn output", images[1].x.size());

  auto it_i = images[0].a.begin();
  auto it_o = images[1].a.begin();
  for (size_t i = 1; i < images.size() / 2; ++i) {
    auto& input_image = images[2 * i];
    auto& output_image = images[2 * i + 1];
    assert(input_image.type == ReluBnImageType::kInput);
    assert(output_image.type == ReluBnImageType::kOutput);
    auto size = input_image.x.size();
    input_image.a.resize(size);
    output_image.a.resize(size);
    std::copy(it_i, it_i + size, input_image.a.begin());
    std::copy(it_o, it_o + size, output_image.a.begin());
    it_i += size;
    it_o += size;
  }
}

inline void ReluBnHandleChallenge(std::vector<ReluBnImage>& images) {
  Tick tick(__FN__);

  // genrate com_xa_r and split
  assert((images.size() - 2) / 2 == kBnLayerInfos.size());
  images[0].com_xa_r = FrRand();
  std::vector<Fr> in_com_xa_r =
      SplitFr(images[0].com_xa_r, kBnLayerInfos.size());
  for (size_t i = 0; i < in_com_xa_r.size(); ++i) {
    images[2 + i * 2].com_xa_r = in_com_xa_r[i];
  }

  images[1].com_xa_r = FrRand();
  std::vector<Fr> out_com_xa_r =
      SplitFr(images[1].com_xa_r, kBnLayerInfos.size());
  for (size_t i = 0; i < out_com_xa_r.size(); ++i) {
    images[2 + i * 2 + 1].com_xa_r = out_com_xa_r[i];
  }

  // parallel compute xa and com_xa
  auto parallel_f = [&images](int64_t i) {
    auto& image = images[i];
    image.xa = InnerProduct(image.a, image.x);
    image.com_xa = pc::PcComputeCommitmentG(image.xa, image.com_xa_r);
  };
  parallel::For((int64_t)images.size(), parallel_f);

#ifdef _DEBUG_CHECK
  Fr chk_in_xa = FrZero();
  Fr chk_out_xa = FrZero();
  G1 chk_in_com_xa = G1Zero();
  G1 chk_out_com_xa = G1Zero();
  for (size_t i = 2; i < images.size(); ++i) {
    if (images[i].type == ReluBnImageType::kInput) {
      chk_in_xa += images[i].xa;
      chk_in_com_xa += images[i].com_xa;
    } else if (images[i].type == ReluBnImageType::kOutput) {
      chk_out_xa += images[i].xa;
      chk_out_com_xa += images[i].com_xa;
    }
  }
  if (images[0].xa != chk_in_xa) throw std::runtime_error("oops 1");
  if (images[1].xa != chk_out_xa) throw std::runtime_error("oops 2");
  if (images[0].com_xa != chk_in_com_xa) throw std::runtime_error("oops 3");
  if (images[1].com_xa != chk_out_com_xa) throw std::runtime_error("oops 4");
#endif
}

inline void ReluBnBuildImages(ProveContext const& context,
                              std::vector<ReluBnImage>& images) {
  Tick tick(__FN__);
  auto const& const_images = context.const_images();
  auto const& com_images = context.image_com_pub().c;
  auto const& com_images_r = context.image_com_sec().r;

  std::vector<Fr> combined_in_x;
  std::vector<Fr> combined_out_x;
  images.resize(kReluBnLayers.size() * 2 + 2);
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    ReluBnImage& in = images[2 + order * 2];
    in.x = const_images[layer]->copy_vec();
    in.type = ReluBnImageType::kInput;
    in.com_x = com_images[layer];
    in.com_x_r = com_images_r[layer];
    combined_in_x.insert(combined_in_x.end(), in.x.begin(), in.x.end());

    ReluBnImage& out = images[2 + order * 2 + 1];
    out.x = const_images[layer + 1]->copy_vec();
    out.type = ReluBnImageType::kOutput;
    out.com_x = com_images[layer + 1];
    out.com_x_r = com_images_r[layer + 1];
    combined_out_x.insert(combined_out_x.end(), out.x.begin(), out.x.end());
  }

  ReluBnImage& combined_in = images[0];
  combined_in.type = ReluBnImageType::kCombinedInput;
  combined_in.x = std::move(combined_in_x);
  combined_in.com_x_r = FrRand();

  ReluBnImage& combined_out = images[1];
  combined_out.type = ReluBnImageType::kCombinedOutput;
  combined_out.x = std::move(combined_out_x);
  combined_out.com_x_r = FrRand();

  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&combined_in]() {
    combined_in.com_x =
        pc::PcComputeCommitmentG(combined_in.x, combined_in.com_x_r);
  };

  tasks[1] = [&combined_out]() {
    combined_out.com_x =
        pc::PcComputeCommitmentG(combined_out.x, combined_out.com_x_r);
  };

  parallel::Invoke(tasks);
}

inline void ReluBnInOutProve(h256_t seed, ProveContext const& context,
                             ReluBnInOutPub& pub, ReluBnInOutSec& sec) {
  Tick tick(__FN__);
  std::vector<ReluBnImage> images;
  ReluBnBuildImages(context, images);

  ReluBnUpdateSeed(seed, images[0].com_x, images[1].com_x);

  ReluBnBuildChallenge(seed, images);

  ReluBnHandleChallenge(images);

  pub.ip_com_pubs.resize(images.size());
  pub.ip_proofs.resize(images.size());

  auto parallel_f = [&seed, &images, &pub](int64_t i) {
    auto const& image = images[i];
    HyraxA::ProveInput input(image.x, image.a, image.xa, pc::kGetRefG,
                             pc::PcG(0));
    pub.ip_com_pubs[i].xi = image.com_x;
    pub.ip_com_pubs[i].tau = image.com_xa;
    HyraxA::CommitmentSec ip_com_sec(image.com_x_r, image.com_xa_r);
    HyraxA::Prove(pub.ip_proofs[i], seed, input, pub.ip_com_pubs[i],
                  ip_com_sec);
  };
  parallel::For(images.size(), parallel_f);

  sec.input = std::move(images[0].x);
  sec.output = std::move(images[1].x);
  sec.com_in_r = images[0].com_x_r;
  sec.com_out_r = images[1].com_x_r;
}

inline void ReluBnBuildPara(ProveContext const& context, std::vector<Fr>& alpha,
                            std::vector<Fr>& beta, std::vector<Fr>& mu) {
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    size_t C = kImageInfos[layer].C;
    size_t D = kImageInfos[layer].D;
    auto const& para = context.para().bn_layer(order);
    if (C != para.alpha.size()) throw std::runtime_error("oops");
    for (size_t j = 0; j < para.alpha.size(); ++j) {
      // repeat DD times
      for (size_t k = 0; k < D * D; ++k) {
        alpha.push_back(para.alpha[j]);
        beta.push_back(para.beta[j]);
        mu.push_back(para.mu[j]);
      }
    }
  }
}

inline size_t ReluBnGetCircuitCount() {
  size_t size = 0;
  for (size_t order = 0; order < kReluBnLayers.size(); ++order) {
    auto layer = kReluBnLayers[order];
    size_t C = kImageInfos[layer].C;
    size_t D = kImageInfos[layer].D;
    size += C * D * D;
  }
  return size;
}

struct ReluBnR1csPub {
  clink::ParallelR1cs<typename R1cs>::Proof r1cs_proof;
  std::vector<G1> com_w;
};

// prove y=<x,para>
inline void ReluBnR1csProve(h256_t seed, ProveContext const& context,
                            ReluBnInOutPub const& inout_pub,
                            ReluBnInOutSec const& inout_sec,
                            ReluBnR1csPub& pub) {
  Tick tick(__FN__);

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::ReluBnGadget<8, 24> gadget(pb, "vgg16 relubn gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  auto r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val  
  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  auto s = r1cs_info->num_variables;
  std::vector<std::vector<Fr>> w(s);
  auto n = inout_sec.input.size();
  for (auto& i : w) i.resize(n);

  std::vector<Fr> alpha;
  alpha.reserve(n);
  std::vector<Fr> beta;
  beta.reserve(n);
  std::vector<Fr> mu;
  mu.reserve(n);
  ReluBnBuildPara(context, alpha, beta, mu);
  
#ifdef _DEBUG_CHECK
  if (alpha.size() != n || beta.size() != n || mu.size() != n) {
    throw std::runtime_error("oops");
  }
  if (n != ReluBnGetCircuitCount()) {
    throw std::runtime_error("oops");
  }

  if (context.para_com_pub().bn.alpha !=
      pc::PcComputeCommitmentG(alpha, context.para_com_sec().bn.alpha_r)) {
    throw std::runtime_error("oops");
  }
  if (context.para_com_pub().bn.beta !=
      pc::PcComputeCommitmentG(beta, context.para_com_sec().bn.beta_r)) {
    throw std::runtime_error("oops");
  }
  if (context.para_com_pub().bn.mu !=
      pc::PcComputeCommitmentG(mu, context.para_com_sec().bn.mu_r)) {
    throw std::runtime_error("oops");
  }

#endif

  for (size_t j = 0; j < n; ++j) {
    gadget.Assign(inout_sec.input[j], alpha[j], beta[j], mu[j]);
    assert(pb.is_satisfied());
    auto v = pb.full_variable_assignment();
    for (int64_t i = 0; i < s; ++i) {
      w[i][j] = v[i];
    }
    if (w[r1cs_ret_index][j] != inout_sec.output[j]) {
      throw std::runtime_error("oops");
    }
  }

  pub.com_w.resize(s);
  std::vector<Fr> com_w_r(s);

  std::cout << "compute conv com(witness)\n";
  assert(r1cs_ret_index == 1);
  com_w_r[0] = inout_sec.com_in_r;
  pub.com_w[0] = inout_pub.ip_com_pubs[0].xi; // in
  
  com_w_r[1] = inout_sec.com_out_r;
  pub.com_w[1] = inout_pub.ip_com_pubs[1].xi; // out

  com_w_r[2] = context.para_com_sec().bn.alpha_r;
  pub.com_w[2] = context.para_com_pub().bn.alpha;

  com_w_r[3] = context.para_com_sec().bn.beta_r;
  pub.com_w[3] = context.para_com_pub().bn.beta;

  com_w_r[4] = context.para_com_sec().bn.mu_r;
  pub.com_w[4] = context.para_com_pub().bn.mu;

  auto parallel_f = [&com_w_r, &pub, &w](int64_t i) {
    std::vector<Fr> const& x = w[i + 5];
    Fr& r = com_w_r[i + 5];
    G1& c = pub.com_w[i + 5];
    r = FrRand();
    c = pc::PcComputeCommitmentG(x, r, true);
  };
  parallel::For<int64_t>(s - 5, parallel_f);

  typename R1cs::ProveInput r1cs_input(*r1cs_info, std::move(w), pub.com_w,
                                       com_w_r, pc::kGetRefG);

  R1cs::Prove(pub.r1cs_proof, seed, std::move(r1cs_input));
}

struct ReluBnProof {
  ReluBnInOutPub inout;
  ReluBnR1csPub r1cs;
};

inline void ReluBnProve(h256_t seed, ProveContext const& context,
                        ReluBnProof& proof) {
  ReluBnInOutPub inout_pub;
  ReluBnInOutSec inout_sec;
  ReluBnInOutProve(seed, context, inout_pub, inout_sec);

  ReluBnR1csPub r1cs_pub;
  ReluBnR1csProve(seed, context, inout_pub, inout_sec, r1cs_pub);

  proof.inout = inout_pub;
  proof.r1cs = r1cs_pub;
}
};  // namespace clink::vgg16
