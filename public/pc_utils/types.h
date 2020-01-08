#pragma once

#include <functional>

#include "ecc/ecc.h"

namespace pc_utils {
typedef std::function<Fr const&(int64_t i, int64_t j)> GetMatrix;
typedef std::function<Fr const&(int64_t i)> GetVector;
}  // namespace pc_utils