#pragma once

#include <cryptopp/osrng.h>

#include "./mpz.h"
#include "./types.h"

namespace misc {
// NOTE: NonblockingRng should enough for linux & windows
// thread_local CryptoPP::AutoSeededRandomPool rng;
inline thread_local CryptoPP::NonblockingRng tls_rng;

inline void RandomBytes(uint8_t* x, uint64_t xlen) {
  tls_rng.GenerateBlock(x, xlen);
}

inline h256_t RandH256() {
  h256_t ret;
  RandomBytes(ret.data(), ret.size());
  return ret;
}

inline mpz_class RandMpz32() {
  h256_t h = RandH256();
  return MpzFromBE(h.data(), h.size());
}

}  // namespace misc