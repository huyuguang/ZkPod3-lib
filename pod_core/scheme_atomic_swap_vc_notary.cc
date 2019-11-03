#include "scheme_atomic_swap_vc_notary.h"
#include "vrs/vrs.h"

namespace scheme::atomic_swap_vc {
bool VerifySecret(Receipt const& receipt, Secret const& secret) {
  return vrs::VerifySecret(receipt.h, receipt.g, receipt.seed0_com,
                           secret.seed0_r, secret.seed0);
}
}  // namespace scheme::atomic_swap_vc
