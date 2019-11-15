#pragma once

#include <stdint.h>
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "../ecc.h"
#include "../ecc_pub.h"
#include "../misc.h"
#include "../multiexp.h"
#include "../parallel.h"
#include "../pds_pub.h"
#include "../tick.h"

namespace hyrax::details {
inline void HashUpdate(CryptoPP::Keccak_256& hash, h256_t const& d) {
  hash.Update((uint8_t const*)d.data(), d.size());
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, G1 const& d) {
  HashUpdate(hash, G1ToBin(d));
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, Fr const& d) {
  HashUpdate(hash, FrToBin(d));
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, std::vector<G1> const& d) {
  for (auto const& i : d) {
    HashUpdate(hash, G1ToBin(i));
  }
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, std::vector<Fr> const& d) {
  for (auto const& i : d) {
    HashUpdate(hash, FrToBin(i));
  }
}

inline G1 ComputeCommitment(std::vector<Fr> const& x, Fr const& r) {
  auto const& pds_pub = GetPdsPub();
  // Tick tick(__FUNCTION__, std::to_string(x.size()));
  assert(PdsPub::kGSize >= x.size());
  auto get_g = [&pds_pub](int64_t i) -> G1 const& {
    return i ? pds_pub.g()[i - 1] : pds_pub.h();
  };
  auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x[i - 1] : r; };
  return MultiExpBdlo12Inner<G1>(get_g, get_f, x.size() + 1);
}

inline G1 ComputeCommitment(Fr const& x, Fr const& r) {
  return ComputeCommitment(std::vector<Fr>{x}, r);
}

inline void VectorMul(std::vector<Fr>& c, std::vector<Fr> const& a, Fr b) {
  c.resize(a.size());
  // for (size_t i = 0; i < a.size(); ++i) {
  //  c[i] = a[i] * b;
  //}
  auto parallel_f = [&c, &a, &b](size_t i) mutable { c[i] = a[i] * b; };
  parallel::For(a.size(), parallel_f, a.size() < 16 * 1024);
}

inline void VectorInc(std::vector<Fr>& a, std::vector<Fr> const& b) {
  assert(a.size() == b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] += b[i];
  }
}
}  // namespace hyrax::details