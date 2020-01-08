#pragma once

#include <string>
#include <vector>

#include "ecc/ecc.h"

namespace scheme::atomic_swap {
struct Request {
  h256_t bob_nonce;
  std::vector<Range> demands;
};

struct Response {
  h256_t alice_nonce;
  std::vector<G1> k;   // n+1
  std::vector<Fr> m;   // n*s
  std::vector<Fr> vw;  // s
};

struct Receipt {
  h256_t seed2;
  Fr sigma_vw;
  uint64_t count;
};

struct Secret {
  h256_t seed0;
};
}  // namespace scheme::atomic_swap
