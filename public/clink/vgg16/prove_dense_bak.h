#pragma once

#include "./context.h"
#include "./image_com.h"

namespace clink::vgg16 {

template <size_t M, size_t N>
struct ProveDenseInput {
  ProveDenseInput(std::vector<std::vector<Fr>> const& weight,
                  std::array<G1, N> const& com_weight,
                  std::array<Fr, N> const& com_weight_r,
                  std::vector<Fr> const& x,
                  G1 const& com_x,
                  Fr const& com_x_r)
      : weight(weight),
        com_weight(com_weight),
        com_weight_r(com_weight_r),
        x(x),
        com_x_r(com_x_r),
        com_x(com_x) {
    namespace fp = circuit::fp;

#ifdef _DEBUG_CHECK
    if (x.size() != M) throw std::runtime_error("invalid x");
    if (weight.size() != N) throw std::runtime_error("invalid weight");
    for (auto const& i : weight) {
      if (i.size() != M+1) throw std::runtime_error("invalid weight");
    }
#endif

    this->x.push_back(fp::RationalConst<8, 24>().kFrN);
    this->com_x += pc::PcG(M) * data.back();
  }
  std::vector<std::vector<Fr>> const& weight;
  std::array<G1, N> const& com_weight;
  std::array<Fr, N> const& com_weight_r;

  std::vector<Fr> x;
  Fr com_x_r;
  G1 com_x;
};

struct DenseProof {
  G1 com;
  G1 com_z;
  groth09::Sec51a::Proof proof_51;
  hyrax::A2::Proof proof_hy;

  bool operator==(DenseProof const& b) const {
    return com == b.com && com_z == b.com_z && proof_51 == b.proof_51 &&
           proof_hy == b.proof_hy;
  }

  bool operator!=(DenseProof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("d.p", ("cd", com), ("cz", com_z), ("51", proof_51),
                       ("hy", proof_hy));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("d.p", ("cd", com), ("cz", com_z), ("51", proof_51),
                       ("hy", proof_hy));
  }
};

template <size_t M, size_t N>
static std::vector<Fr> ComputeDenseFst(h256_t const& seed) {
  std::vector<Fr> x(N);
  std::string str = "dense";
  str += std::to_string(M);
  str += "_";
  str += std::to_string(N);
  ComputeFst(seed, str, x);
  return x;
}

struct ProveOutput {
  Fr com_r;
  G1 com;
  std::vector<Fr> data;
};

template <size_t M, size_t N>
static void ProveDense(DenseProof& proof, ProveOutput& output, h256_t seed,
                       ProveDenseInput<M, N> const& input) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  assert(input.weight.size() == N);
  assert(input.data.size() == M + 1);

  if (M == 10) {
    std::cout << "data:\n";
    for (auto const& i : input.data) {
      double ret = fp::RationalToDouble<8, 24>(i);
      std::cout << std::right << std::setw(12) << std::setfill(' ') << ret;
    }
    std::cout << "\n";
  }

  // build output
  output.data.resize(N);
  auto parallel_f = [&input, &output](int64_t i) {
    assert(input.weight[i].size() == M + 1);
    output.data[i] = std::inner_product(input.data.begin(), input.data.end(),
                                        input.para_dense[i].begin(), FrZero());
  };
  parallel::For(N, parallel_f);

  //#ifdef _DEBUG
  std::cout << "dense before relu:\n";
  for (size_t i = 0; i < N; ++i) {
    double ret = fp::RationalToDouble<8, 24 + 24>(output.data[i]);
    std::cout << std::right << std::setw(12) << std::setfill(' ') << ret;
  }
  std::cout << "\n";
  //#endif

  output.com_r = FrRand();
  output.com = pc::PcComputeCommitmentG(output.data, output.com_r);

  // prove
  std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
  G1 com_e = G1Zero();
  Fr com_e_r = FrZero();
  // TODO: change to multiexp and inner_product
  for (size_t i = 0; i < N; ++i) {
    com_e += input.com_weight[i] * x[i];
    com_e_r += input.com_weight_r[i] * x[i];
  }
  std::vector<Fr> e(M + 1);
  for (size_t i = 0; i < M + 1; ++i) {
    e[i] = FrZero();
    for (size_t j = 0; j < N; ++j) {
      e[i] += input.weight[j][i] * x[j];
    }
  }

  Fr z = std::inner_product(x.begin(), x.end(), output.data.begin(), FrZero());
  Fr com_z_r = FrRand();
  G1 com_z = pc::PcComputeCommitmentG(z, com_z_r);
  proof.com_z = com_z;
  proof.com = output.com;

#ifdef _DEBUG
  assert(com_e == pc::PcComputeCommitmentG(e, com_e_r));
  assert(z ==
         std::inner_product(e.begin(), e.end(), input.data.begin(), FrZero()));
#endif

  // prove left
  hyrax::A2::ProveInput input_hy(output.data, x, z, pc::kGetRefG, pc::PcG(0));
  hyrax::A2::CommitmentPub com_pub_hy(output.com, com_z);
  hyrax::A2::CommitmentSec com_sec_hy(output.com_r, com_z_r);
  hyrax::A2::Prove(proof.proof_hy, seed, input_hy, com_pub_hy, com_sec_hy);

  // prove right
  // std::cout << "prove, seed: " << misc::HexToStr(seed) << "\n";
  // std::cout << "prove, com_e: " << com_e << "\n";
  // std::cout << "prove, com_data: " << input.com_data << "\n";
  // std::cout << "prove, com_z: " << com_z << "\n";
  groth09::Sec51a::ProveInput input_51(e, input.data, z, pc::kGetRefG,
                                       pc::kGetRefG, pc::PcG(0));
  groth09::Sec51a::CommitmentPub com_pub_51(com_e, input.com_data, com_z);
  groth09::Sec51a::CommitmentSec com_sec_51(com_e_r, input.com_data_r, com_z_r);
  groth09::Sec51a::Prove(proof.proof_51, seed, input_51, com_pub_51,
                         com_sec_51);
}

template <size_t M, size_t N>
struct VerifyDenseInput {
  VerifyDenseInput(G1 const& com_data, std::array<G1, N> const& com_weight)
      : com_data(com_data), com_weight(com_weight) {
    namespace fp = circuit::fp;
    this->com_data += pc::PcG(M) * fp::RationalConst<8, 24>().kFrN;
  }
  G1 com_data;
  std::array<G1, N> const& com_weight;
};

template <size_t M, size_t N>
static bool VerifyDense(DenseProof const& proof, h256_t seed,
                        VerifyDenseInput<M, N> const& input) {
  Tick tick(__FN__);
  std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
  G1 com_e = G1Zero();
  // TODO: change to multiexp
  for (size_t i = 0; i < N; ++i) {
    com_e += input.com_weight[i] * x[i];
  }

  hyrax::A2::CommitmentPub com_pub_hy(proof.com, proof.com_z);
  hyrax::A2::VerifyInput input_hy(x, com_pub_hy, pc::kGetRefG, pc::PcG(0));
  if (!hyrax::A2::Verify(proof.proof_hy, seed, input_hy)) {
    assert(false);
    return false;
  }

  // std::cout << "verify, seed: " << misc::HexToStr(seed) << "\n";
  // std::cout << "verify, com_e: " << com_e << "\n";
  // std::cout << "verify, com: " << input.com << "\n";
  // std::cout << "verify, com_z: " << proof.com_z << "\n";
  groth09::Sec51a::CommitmentPub com_pub_51(com_e, input.com_data, proof.com_z);
  groth09::Sec51a::VerifyInput input_51(com_pub_51, pc::kGetRefG, pc::kGetRefG,
                                        pc::PcG(0));
  if (!groth09::Sec51a::Verify(proof.proof_51, seed, input_51)) {
    assert(false);
    return false;
  }

  return true;
}
}  // namespace clink::vgg16