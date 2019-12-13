#pragma once

#include <memory>
#include <string>

#include "basic_types.h"
#include "scheme_atomic_swap_protocol.h"

namespace scheme::atomic_swap {
bool VerifySecret(uint64_t s, Receipt const& receipt, Secret const& secret);
bool VerifySecret(uint64_t s, uint64_t count, Fr const& sigma_vw,
                  std::vector<Fr> const& v, std::vector<Fr> const& w);
}  // namespace scheme::atomic_swap
