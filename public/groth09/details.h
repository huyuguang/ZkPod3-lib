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

inline void PrintVector(std::vector<Fr> const& a) {
  std::cout << "\n";
  for (auto const& i : a) {
    std::cout << i << "\n";
  }
  std::cout << "\n";
}

}  // namespace groth09::details