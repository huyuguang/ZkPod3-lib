#pragma once

#include "./types.h"

namespace pod {
inline bool VerifySecret(Receipt const& receipt, Secret const& secret) {
  return vrs::VerifySecret(receipt.h, receipt.g, receipt.seed0_com,
                           secret.seed0_com_r, secret.seed0);
}
}  // namespace pod
