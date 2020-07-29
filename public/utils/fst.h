#pragma once

#include "ecc/ecc.h"
#include "public.h"

inline void HashUpdate(CryptoPP::Keccak_256& hash, uint64_t d) {
  auto big_d = boost::endian::native_to_big(d);
  hash.Update((uint8_t const*)&big_d, sizeof(big_d));
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, std::string const& d) {
  hash.Update((uint8_t const*)d.data(), d.size());
}

inline void HashUpdate(CryptoPP::Keccak_256& hash, void const* p, size_t len) {
  hash.Update((uint8_t const*)p, len);
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

template <typename T>
inline void HashUpdate(CryptoPP::Keccak_256& hash, std::vector<T> const& d) {
  for (auto const& i : d) {
    HashUpdate(hash, i);
  }
}

// Fiat-Shamir transform
// c[i] = hash(seed,salt,i)
inline void ComputeFst1(h256_t const& seed, std::string const& salt,
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

// Fiat-Shamir transform
// c[i] = hash(seed, salt)^i
inline void ComputeFst2(h256_t const& seed, std::string const& salt,
                        std::vector<Fr>& c) {
  assert(c.size() >= 2);
  assert(!c.empty());
  CryptoPP::Keccak_256 hash;
  h256_t digest;
  HashUpdate(hash, seed);
  HashUpdate(hash, salt);
  hash.Final(digest.data());
  c[0] = 1;
  c[1] = H256ToFr(digest);
  H256ToFr(digest);
  for (size_t i = 2; i < c.size(); ++i) {
    c[i] = c[i - 1] * c[1];
  }
}

inline void ComputeFst(h256_t const& seed, std::string const& salt,
                       std::vector<Fr>& c) {
  return ComputeFst1(seed, salt, c);
}