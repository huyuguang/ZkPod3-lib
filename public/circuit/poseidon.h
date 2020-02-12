#pragma once

#include <cryptopp/keccak.h>

#include <vector>

#include "ecc/ecc.h"

namespace circuit {
struct PoseidonConstants {
  std::vector<Fr> C;  // `t` constants
  std::vector<Fr> M;  // `t * t` matrix of constants
  bool operator==(PoseidonConstants const& right) const {
    return C == right.C && M == right.M;
  }
  bool operator!=(PoseidonConstants const& right) const {
    return !(*this == right);
  }
};

inline void poseidon_constants_fill(const std::string& seed,
                                    unsigned n_constants,
                                    std::vector<Fr>& result) {
  result.resize(n_constants);
  for (size_t i = 0; i < n_constants; ++i) {
    h256_t digest;
    CryptoPP::Keccak_256 hash;
    std::string s = seed + std::to_string(i);
    hash.Update((uint8_t const*)s.data(), s.size());
    hash.Final(digest.data());

    bool success;
    result[i].setArray(&success, digest.data(), digest.size(), mcl::fp::Mod);
    assert(success);
  }
}

inline const std::vector<Fr> poseidon_constants(const std::string& seed,
                                                unsigned n_constants) {
  std::vector<Fr> result;
  poseidon_constants_fill(seed, n_constants, result);
  return result;
}

inline void poseidon_matrix_fill(const std::string& seed, unsigned t,
                                 std::vector<Fr>& result) {
  const std::vector<Fr> c = poseidon_constants(seed, t * 2);

  result.reserve(t * 2);

  for (unsigned i = 0; i < t; i++) {
    for (unsigned j = 0; j < t; j++) {
      result.emplace_back((c[i] - c[t + j]).inverse());
    }
  }
}

inline const std::vector<Fr> poseidon_matrix(const std::string& seed,
                                             unsigned t) {
  std::vector<Fr> result;
  poseidon_matrix_fill(seed, t, result);
  return result;
}

template <unsigned param_t, unsigned param_F, unsigned param_P>
const PoseidonConstants& poseidon_params() {
  static PoseidonConstants constants;
  static std::once_flag flag;

  std::call_once(flag, []() {
    poseidon_constants_fill("poseidon_constants", param_F + param_P,
                            constants.C);
    poseidon_matrix_fill("poseidon_matrix_0000", param_t, constants.M);
  });

  return constants;
}

template <unsigned param_t, unsigned nSBox, unsigned nInputs, unsigned nOutputs>
std::array<Fr, nOutputs> PoseidonRound(const Fr& C_i, const std::vector<Fr>& M,
                              const std::array<Fr, nInputs>& state) {
  std::array<Fr, nSBox> sbox;
  for (unsigned h = 0; h < nSBox; h++) {
    auto value = C_i;
    if (h < nInputs) {
      value += state[h];
    }
    auto value2 = value * value;
    sbox[h] = value2 * value2 * value;
  }

  std::array<Fr, nOutputs> ret;

  for (unsigned i = 0; i < nOutputs; i++) {
    const unsigned M_offset = i * param_t;

    // Any element which isn't passed through an sbox
    // Can be accumulated separately as part of the constant term
    Fr& lc = ret[i];
    lc = FrZero();
    if (nSBox < param_t) {
      for (unsigned j = nSBox; j < param_t; j++) {
        lc += C_i * M[M_offset + j];
      }
    }

    // Add S-Boxes to the output row
    for (unsigned s = 0; s < nSBox; s++) {
      lc += sbox[s] * M[M_offset + s];
    }

    // Then add inputs (from the state) multiplied by the matrix element
    for (unsigned k = nSBox; k < nInputs; k++) {
      lc += state[k] * M[M_offset + k];
    }
  }
  return ret;
}

template <unsigned param_t, unsigned param_c, unsigned param_F,
          unsigned param_P, unsigned nInputs, unsigned nOutputs>
std::array<Fr, nOutputs> Poseidon(std::array<Fr, nInputs> const& inputs) {
  static constexpr unsigned partial_begin = (param_F / 2);
  static constexpr unsigned partial_end = (partial_begin + param_P);
  static constexpr unsigned total_rounds = param_F + param_P;

  auto const& constants = poseidon_params<param_t, param_F, param_P>();

  // first round
  auto first_round = PoseidonRound<param_t, param_t, nInputs, param_t>(
      constants.C[0], constants.M, inputs);
  //misc::PrintArray(first_round);

  // prefix_full_round
  auto prefix_full_round = PoseidonRound<param_t, param_t, param_t, param_t>(
      constants.C[1], constants.M, first_round);

  for (size_t i = 1 + 1; i < partial_begin; ++i) {
    prefix_full_round = PoseidonRound<param_t, param_t, param_t, param_t>(
        constants.C[i], constants.M, prefix_full_round);
  }
  //misc::PrintArray(prefix_full_round);

  // partial_round
  auto partial_round = PoseidonRound<param_t, param_c, param_t, param_t>(
      constants.C[partial_begin], constants.M, prefix_full_round);

  for (size_t i = partial_begin + 1; i < partial_end; ++i) {
    partial_round = PoseidonRound<param_t, param_c, param_t, param_t>(
        constants.C[i], constants.M, partial_round);
  }
  //misc::PrintArray(partial_round);

  // suffix_full_round
  auto suffix_full_round = PoseidonRound<param_t, param_t, param_t, param_t>(
      constants.C[partial_end], constants.M, partial_round);
  for (size_t i = partial_end + 1; i < total_rounds - 1; ++i) {
    suffix_full_round = PoseidonRound<param_t, param_t, param_t, param_t>(
        constants.C[i], constants.M, suffix_full_round);
  }
  //misc::PrintArray(suffix_full_round);

  // last round
  auto last_round = PoseidonRound<param_t, param_t, param_t, nOutputs>(
      constants.C.back(), constants.M, suffix_full_round);

  return last_round;
}
}  // namespace circuit
