#pragma once

#include "./types.h"
#include "ecc/ecc.h"
#include "parallel/parallel.h"
#include "utils/fst.h"

namespace pod::details {
inline bool CheckCommitedData(CommitedData const& data) {
  if (!data.n || !data.s || data.s > PcBase::kGSize) return false;

  bool all_success = false;
  auto parallel_f = [&data](int64_t i) {
    auto x = [&data, i](int64_t j) -> Fr const& { return data.get_m(i, j); };
    auto const& r = data.get_r(i);
    return PcComputeCommitment(data.s, x, r) == data.get_com(i);
  };
  parallel::For(&all_success, data.n, parallel_f);
  return all_success;
}

inline void UpdateSeed(h256_t& seed, std::vector<G1> const& k) {
  // Tick tick(__FUNCTION__);
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, k);
  hash.Final(seed.data());
}

inline bool CheckVrs(int64_t n, int64_t s, h256_t const& seed,
                     ProvedData const& proved_data, VerifyOutput& output,
                     vrs::VerifyOutput& vrs_output) {
  output.plain.resize((n + 1) * (s + 1));
  auto parallel_f = [&output, &proved_data](int64_t i) {
    output.plain[i] = vrs::GeneratePlain(proved_data.vrs_plain_seed, i);
  };
  parallel::For((int64_t)output.plain.size(), parallel_f);

  vrs::PublicInput public_input(
      output.plain.size(), [&output](int64_t i) { return output.plain[i]; });

  vrs::LargeVerifier verifier(public_input);
  auto get_w = [&output, s](int64_t i) { return output.w[i / (s + 1)]; };
  if (!verifier.Verify(seed, get_w, proved_data.vrs_proofs, vrs_output)) {
    assert(false);
    return false;
  }

  auto const& vw = proved_data.vw;
  output.sigma_vw = parallel::Accumulate(vw.begin(), vw.end(), FrZero());
  G1 check_vw_com = PcComputeCommitment(output.sigma_vw, proved_data.vw_com_r);
  if (verifier.com_vw() != check_vw_com) {
    assert(false);
    return false;
  }
  return true;
}
}  // namespace pod::details