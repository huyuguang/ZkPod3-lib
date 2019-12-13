#pragma once

#include <stdint.h>

#include <algorithm>
#include <boost/endian/conversion.hpp>
#include <functional>
#include <memory>
#include <vector>

#include "ecc.h"
#include "ecc_pub.h"

inline void HashUpdate(CryptoPP::Keccak_256& hash, uint64_t d) {
  auto big_d = boost::endian::native_to_big(d);
  hash.Update((uint8_t const*)&big_d, sizeof(big_d));
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, std::string const& d) {
  hash.Update((uint8_t const*)d.data(), d.size());
}

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

// Fiat-Shamir transform
inline void ComputeFst(h256_t const& seed, std::string const& salt,
                       std::vector<Fr>& c) {
  assert(!c.empty());
  auto parallel_f = [&seed, &c, &salt](int64_t i) {
    CryptoPP::Keccak_256 hash;
    h256_t digest;
    HashUpdate(hash, seed);
    HashUpdate(hash, salt);
    HashUpdate(hash, i);
    hash.Final(digest.data());
    c[i] = H256ToFr(digest);
  };
  parallel::For(c.size(), parallel_f, c.size() < 16 * 1024);
}
