#pragma once

#include "./dense_prove.h"

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
  HyraxA::VerifyInput input_hy("dense", x, com_pub_hy, pc::kGetRefG1,
                               pc::PcG(0));
  CHECK(HyraxA::Verify(proof.proof_hy, seed, input_hy), "");

  // std::cout << "verify, seed: " << misc::HexToStr(seed) << "\n";
  // std::cout << "verify, com_e: " << com_e << "\n";
  // std::cout << "verify, com_y: " << input.com_y << "\n";
  // std::cout << "verify, com_z: " << proof.com_z << "\n";
  Sec51::CommitmentPub com_pub_51(com_e, input.com_x, proof.com_z);
  std::vector<Fr> t(M + 1, FrOne());
  Sec51::VerifyInput input_51(t, com_pub_51, pc::kGetRefG1, pc::kGetRefG1,
                              pc::PcG(0));
  CHECK(Sec51::Verify(proof.proof_51, seed, input_51), "");

  return true;
}

template <size_t kOrder>
bool DenseVerify(h256_t seed, VerifyContext const& context,
                 DenseProof const& proof) {
  Tick tick(__FN__);

  constexpr size_t kLayer = kDenseLayers[kOrder];
  constexpr ImageInfo const& info_in = kImageInfos[kLayer];
  constexpr ImageInfo const& info_out = kImageInfos[kLayer + 1];
  constexpr size_t M = info_in.C * info_in.D * info_in.D;
  constexpr size_t N = info_out.C * info_out.D * info_out.D;

  auto const& com_weight = context.para_com_pub().dense.get<kOrder>();
  auto const& com_x = context.image_com_pub().c[kLayer];
  VerifyDenseInput<M, N> input(com_x, com_weight);
  return VerifyDense<M, N>(proof, seed, input);
}

}  // namespace clink::vgg16