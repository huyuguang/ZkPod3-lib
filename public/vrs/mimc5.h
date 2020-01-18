#pragma once

#include <cryptopp/keccak.h>
#include <vector>
#include "ecc/ecc.h"

namespace vrs {

// P = Fr::getModulo(), that is the order of the  G1, equal
// 21888242871839275222246405745257275088548364400416034343698204186575808495617
// gcd(P-1, 5) = 1, gcd(P-1, 7) = 1, gcd(P-1, 3)=0
// so use 5 or 7 is safe.

namespace details {
inline std::vector<Fr> ComputeConst(std::string const& prefix, size_t count) {
  std::vector<Fr> ret(count);
  for (size_t i = 0; i < count; ++i) {
    h256_t digest;
    CryptoPP::Keccak_256 hash;
    std::string s = prefix + std::to_string(i);
    hash.Update((uint8_t const*)s.data(), s.size());
    hash.Final(digest.data());

    bool success;
    ret[i].setArray(&success, digest.data(), digest.size(), mcl::fp::Mod);
    assert(success);
  }
  return ret;
}
}  // namespace details

inline static const int64_t kMimc5Round = 40;

inline std::vector<Fr> const& Mimc5Const() {
  static auto instance = details::ComputeConst("mimc_5_const", kMimc5Round);
  return instance;
}

// TODO: wrapped by OneWayGadget
inline Fr Mimc5Enc(Fr const& plain, Fr const& key) {
  auto const& kMimc5Const = Mimc5Const();
  auto box = [](Fr const& v) {
    auto vv = v * v;
    return vv * vv * v;
  };
  Fr enc = plain;
  for (size_t i = 0; i < kMimc5Const.size(); ++i) {
    enc = box(enc + key + kMimc5Const[i]);
  }
  return enc + key;
}
}  // namespace vrs
