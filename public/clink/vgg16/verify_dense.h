#pragma once

#include "./prove_dense.h"

namespace clink::vgg16 {

template <size_t M, size_t N>
struct VerifyDenseInput {
  VerifyDenseInput(G1 const& com_x, std::array<G1, N> const& com_weight)
      : com_x(com_x), com_weight(com_weight) {
    namespace fp = circuit::fp;
    this->com_x += pc::PcG(M) * fp::RationalConst<8, 24>().kFrN;
  }
  G1 com_x;
  std::array<G1, N> const& com_weight;
};

template <size_t M, size_t N>
static bool VerifyDense(DenseProof const& proof, h256_t seed,
                        VerifyDenseInput<M, N> const& input) {
  Tick tick(__FN__);
  std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
  G1 com_e = G1Zero();
  for (size_t i = 0; i < N; ++i) {
    com_e += input.com_weight[i] * x[i];
  }

  HyraxA::CommitmentPub com_pub_hy(proof.com_y, proof.com_z);
  HyraxA::VerifyInput input_hy(x, com_pub_hy, pc::kGetRefG, pc::PcG(0));
  if (!HyraxA::Verify(proof.proof_hy, seed, input_hy)) {
    assert(false);
    return false;
  }

  // std::cout << "verify, seed: " << misc::HexToStr(seed) << "\n";
  // std::cout << "verify, com_e: " << com_e << "\n";
  // std::cout << "verify, com_y: " << input.com_y << "\n";
  // std::cout << "verify, com_z: " << proof.com_z << "\n";
  groth09::Sec51a::CommitmentPub com_pub_51(com_e, input.com_x, proof.com_z);
  groth09::Sec51a::VerifyInput input_51(com_pub_51, pc::kGetRefG, pc::kGetRefG,
                                        pc::PcG(0));
  if (!groth09::Sec51a::Verify(proof.proof_51, seed, input_51)) {
    assert(false);
    return false;
  }

  return true;
}

inline bool DenseVerify0(h256_t seed, VerifyContext const& context,
                        DenseProof const& proof) {
  Tick tick(__FN__);

  constexpr size_t kOrder = 0;
  constexpr size_t kLayer = 31;
  constexpr size_t M = 512;
  constexpr size_t N = 512;

#ifdef _DEBUG_CHECK
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first == kDense &&
        kLayerTypeOrders[i].second == kOrder) {
      if (kLayer != i) throw std::runtime_error("oops");
      break;
    }
  }

  auto const& info = kImageInfos[kLayer];
  if (info.channel_count * info.dimension * info.dimension != M)
    throw std::runtime_error("oops");
  auto const& info2 = kImageInfos[kLayer+1];
  if (info2.channel_count * info2.dimension * info2.dimension != N)
    throw std::runtime_error("oops");
#endif

  auto const& com_weight = context.para_com_pub().dense.d0;
  auto const& com_x = context.image_com_pub().c[kLayer];
  VerifyDenseInput<M, N> input(com_x,com_weight);
  return VerifyDense<M, N>(proof, seed, input);
}


inline bool DenseVerify1(h256_t seed, VerifyContext const& context,
                        DenseProof const& proof) {
  Tick tick(__FN__);

  constexpr size_t kOrder = 1;
  constexpr size_t kLayer = 33;
  constexpr size_t M = 512;
  constexpr size_t N = 10;

#ifdef _DEBUG_CHECK
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first == kDense &&
        kLayerTypeOrders[i].second == kOrder) {
      if (kLayer != i) throw std::runtime_error("oops");
      break;
    }
  }

  auto const& info = kImageInfos[kLayer];
  if (info.channel_count * info.dimension * info.dimension != M)
    throw std::runtime_error("oops");
  auto const& info2 = kImageInfos[kLayer+1];
  if (info2.channel_count * info2.dimension * info2.dimension != N)
    throw std::runtime_error("oops");
#endif

  auto const& com_weight = context.para_com_pub().dense.d1;
  auto const& com_x = context.image_com_pub().c[kLayer];
  VerifyDenseInput<M, N> input(com_x,com_weight);
  return VerifyDense<M, N>(proof, seed, input);
}

}  // namespace clink::vgg16