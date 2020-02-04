#pragma once

#include <cryptopp/keccak.h>

#include <boost/endian/conversion.hpp>
#include <string>
#include <vector>

#include "./mimc5.h"
#include "ecc/ecc.h"

namespace vrs {

inline static const std::string kFstHpCom = "vrs_fst_hp_com";

inline void GeneratePlain(Fr* out, h256_t const& plain_seed, int64_t position) {
  CryptoPP::Keccak_256 hash;
  h256_t plain;
  hash.Update(plain_seed.data(), plain_seed.size());
  auto offset_big = boost::endian::native_to_big(position);
  hash.Update((uint8_t const*)&offset_big, sizeof(offset_big));
  hash.Final(plain.data());
  *out = H256ToFr(plain);
}

inline Fr GeneratePlain(h256_t const& plain_seed, int64_t position) {
  Fr ret;
  GeneratePlain(&ret, plain_seed, position);
  return ret;
}

// Since the plain data is random, we do not need to use cbc mode.
inline void GenerateV(int64_t offset, Fr const& key, h256_t const& plain_seed,
                      Fr* out) {
  auto plain = GeneratePlain(plain_seed, offset);
  *out = Mimc5Enc(plain, key);
}

inline Fr GenerateV(int64_t offset, Fr const& key, h256_t const& plain_seed) {
  Fr ret;
  GenerateV(offset, key, plain_seed, &ret);
  return ret;
}

inline std::vector<Fr> SplitFr(Fr const& fr, int64_t count) {
  std::vector<Fr> ret(count);
  if (count == 1) {
    ret[0] = fr;
  } else {
    for (int64_t i = 0; i < count - 1; ++i) {
      ret[i] = FrRand();
    }

    auto sum_exclude_last =
        std::accumulate(ret.begin(), ret.begin() + count - 1, FrZero());
    ret.back() = fr - sum_exclude_last;

    assert(fr == std::accumulate(ret.begin(), ret.end(), FrZero()));
  }
  return ret;
}

template <typename Output>
void MergeOutputs(Output& output, std::vector<Output> const& outputs) {
  output.h = outputs[0].h;
  output.g = parallel::Accumulate(
      outputs.begin(), outputs.end(), G1Zero(),
      [](G1 const& a, Output const& b) { return a + b.g; });
  output.key_com = parallel::Accumulate(
      outputs.begin(), outputs.end(), G1Zero(),
      [](G1 const& a, Output const& b) { return a + b.key_com; });
}

}  // namespace vrs