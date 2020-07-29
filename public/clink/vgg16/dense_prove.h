#pragma once

#include "./context.h"
#include "./image_com.h"

namespace clink::vgg16 {

template <size_t M, size_t N>
struct ProveDenseInput {
  ProveDenseInput(std::vector<std::vector<Fr>> const& weight,
                  std::array<G1, N> const& com_weight,
                  std::array<Fr, N> const& com_weight_r, std::vector<Fr>&& x,
                  G1 const& com_x, Fr const& com_x_r, std::vector<Fr>&& y,
                  G1 const& com_y, Fr const& com_y_r)
      : weight(weight),
        com_weight(com_weight),
        com_weight_r(com_weight_r),
        x(std::move(x)),
        com_x_r(com_x_r),
        com_x(com_x),
        y(std::move(y)),
        com_y_r(com_y_r),
        com_y(com_y) {
    namespace fp = circuit::fp;
    this->x.push_back(fp::RationalConst<8, 24>().kFrN);
    this->com_x += pc::PcG(M) * this->x.back();

#ifdef _DEBUG_CHECK
    if (this->x.size() != M + 1) throw std::runtime_error("invalid x");
    if (this->y.size() != N) throw std::runtime_error("invalid y");
    if (weight.size() != N) throw std::runtime_error("invalid weight");
    for (auto const& i : weight) {
      if (i.size() != M + 1) throw std::runtime_error("invalid weight");
    }

    bool all_success = false;
    auto parallel_f = [this](int64_t i) {
      auto check_y = InnerProduct(this->x, this->weight[i]);
      return check_y == this->y[i];
    };
    parallel::For(&all_success, N, parallel_f);
    if (!all_success) throw std::runtime_error("invalid y");

    if (pc::ComputeCom(this->y, this->com_y_r) != this->com_y) {
      std::runtime_error("invalid com_y");
    }

#endif
  }
  std::vector<std::vector<Fr>> const& weight;
  std::array<G1, N> const& com_weight;
  std::array<Fr, N> const& com_weight_r;

  std::vector<Fr> x;
  Fr com_x_r;
  G1 com_x;
  std::vector<Fr> y;
  Fr com_y_r;
  G1 com_y;
};

struct DenseProof {
  G1 com_y;
  G1 com_z;
  Sec51::Proof proof_51;
  HyraxA::Proof proof_hy;

  bool operator==(DenseProof const& b) const {
    return com_y == b.com_y && com_z == b.com_z && proof_51 == b.proof_51 &&
           proof_hy == b.proof_hy;
  }

  bool operator!=(DenseProof const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("d.p", ("cd", com_y), ("cz", com_z), ("51", proof_51),
                       ("hy", proof_hy));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("d.p", ("cd", com_y), ("cz", com_z), ("51", proof_51),
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

template <size_t M, size_t N>
static void ProveDense(h256_t seed, ProveDenseInput<M, N> const& input,
                       DenseProof& proof) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;

  // prove
  std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
  G1 com_e = G1Zero();
  Fr com_e_r = FrZero();

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

  Fr z = InnerProduct(x, input.y);
  Fr com_z_r = FrRand();
  proof.com_z = pc::ComputeCom(z, com_z_r);
  proof.com_y = input.com_y;

#ifdef _DEBUG_CHECK
  assert(com_e == pc::ComputeCom(e, com_e_r));
  assert(z == InnerProduct(e, input.x));
#endif

  // prove left
  HyraxA::ProveInput input_hy(input.y, x, z, pc::kGetRefG1, pc::PcG(0));
  HyraxA::CommitmentPub com_pub_hy(input.com_y, proof.com_z);
  HyraxA::CommitmentSec com_sec_hy(input.com_y_r, com_z_r);
  HyraxA::Prove(proof.proof_hy, seed, input_hy, com_pub_hy, com_sec_hy);

  // prove right
  std::vector<Fr> t(e.size(), FrOne());
  Sec51::ProveInput input_51(e, input.x, t, input.x, z, pc::kGetRefG1,
                             pc::kGetRefG1, pc::PcG(0));
  Sec51::CommitmentPub com_pub_51(com_e, input.com_x, proof.com_z);
  Sec51::CommitmentSec com_sec_51(com_e_r, input.com_x_r, com_z_r);
  Sec51::Prove(proof.proof_51, seed, input_51, com_pub_51, com_sec_51);
}

template <size_t kOrder>
void DenseProve(h256_t seed, ProveContext const& context, DenseProof& proof) {
  Tick tick(__FN__);

  constexpr size_t kLayer = kDenseLayers[kOrder];
  constexpr ImageInfo const& info_in = kImageInfos[kLayer];
  constexpr ImageInfo const& info_out = kImageInfos[kLayer + 1];
  constexpr size_t M = info_in.C * info_in.D * info_in.D;
  constexpr size_t N = info_out.C * info_out.D * info_out.D;

  auto const& para = context.para().dense_layer(kOrder);
  auto const& weight = para.weight;
  auto const& com_weight = context.para_com_pub().dense.get<kOrder>();
  auto const& com_weight_r = context.para_com_sec().dense.get<kOrder>();
  std::vector<Fr> x = context.const_images()[kLayer]->data;
  auto const& com_x = context.image_com_pub().c[kLayer];
  auto const& com_x_r = context.image_com_sec().r[kLayer];
  std::vector<Fr> y = context.const_images()[kLayer + 1]->data;
  auto const& com_y = context.image_com_pub().c[kLayer + 1];
  auto const& com_y_r = context.image_com_sec().r[kLayer + 1];

  ProveDenseInput<M, N> input(weight, com_weight, com_weight_r, std::move(x),
                              com_x, com_x_r, std::move(y), com_y, com_y_r);

  ProveDense<M, N>(seed, input, proof);
}
}  // namespace clink::vgg16