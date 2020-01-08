#pragma once

#include <string>
#include <vector>

#include "ecc/ecc.h"

namespace scheme::complaint {
struct Request {
  h256_t bob_nonce;
  std::vector<Range> demands;
};

struct Response {
  h256_t alice_nonce;
  std::vector<G1> k;
  std::vector<Fr> m;
};

struct Receipt {
  h256_t seed2;
  h256_t k_mkl_root;
  uint64_t count;
};

struct Secret {
  h256_t seed0;
};

struct Claim {
  uint64_t i;
  G1 ki;
  std::vector<h256_t> mkl_path;
};
}  // namespace scheme::complaint
