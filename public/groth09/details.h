#pragma once

#include <stdint.h>

#include <algorithm>
#include <boost/endian/conversion.hpp>
#include <functional>
#include <memory>
#include <vector>

#include "ecc/ecc.h"
#include "hyrax/hyrax.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "parallel/parallel.h"
#include "utils/fst.h"

namespace groth09::details {

inline void ComputePowOfE(Fr const& e, int64_t m, std::vector<Fr>& vec,
                          std::vector<Fr>& rev) {
  // Tick tick(__FUNCTION__, std::to_string(m));
  vec.resize(m * 2 - 1);
  rev.resize(m);

  vec[0] = FrOne();
  for (int64_t i = 1; i < m * 2 - 1; ++i) {
    vec[i] = e * vec[i - 1];
  }
  for (int64_t i = 0; i < m; ++i) {
    rev[i] = vec[m - i - 1];
  }
}

inline void HadamardProduct(std::vector<Fr>& c, std::vector<Fr> const& a,
                            std::vector<Fr> const& b) {
  assert(a.size() == b.size());
  c.resize(a.size());
  auto parallel_f = [&c, &a, &b](size_t i) { c[i] = a[i] * b[i]; };
  parallel::For(a.size(), parallel_f, true);//a.size() < 16 * 1024);
}

inline std::vector<Fr> HadamardProduct(std::vector<Fr> const& a,
                                       std::vector<Fr> const& b) {
  assert(a.size() == b.size());
  std::vector<Fr> c(a.size());
  HadamardProduct(c, a, b);
  return c;
}

inline void PrintVector(std::vector<Fr> const& a) {
  std::cout << "\n";
  for (auto const& i : a) {
    std::cout << i << "\n";
  }
  std::cout << "\n";
}

inline void BuildChallengeVector(std::vector<Fr>& challenge, h256_t const& seed,
                                 std::string const& suffix, int64_t count) {
  auto hash = [&suffix, &seed](int64_t i) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    hash.Update((uint8_t const*)&i, sizeof(i));
    HashUpdate(hash, suffix);
    h256_t digest;
    hash.Final(digest.data());
    return H256ToFr(digest);
  };

  challenge.resize(count);
  for (int64_t i = 0; i < count; ++i) {
    challenge[i] = hash(i);
  }
}

}  // namespace groth09::details