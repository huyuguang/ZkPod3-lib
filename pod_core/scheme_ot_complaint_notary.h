#pragma once

#include <memory>
#include <string>
#include "basic_types.h"
#include "chain.h"
#include "mkl_tree.h"
#include "scheme_misc.h"

namespace scheme::ot_complaint {
template <typename Receipt, typename Secret, typename Claim>
bool VerifyClaim(uint64_t s, Receipt const& receipt, Secret const& secret,
                 Claim const& claim) {
  if (!VerifyPathOfK(claim.ki, claim.i, receipt.count + 1,
                     receipt.k_mkl_root, claim.mkl_path)) {
    assert(false);
    return false;
  }

  // NOTE: Blockchain vm does not have ecc pub, must call u^v directly
  std::vector<Fr> v(s);
  ChainKeccak256(secret.seed0, claim.i * s, claim.i * s + s, v);

  G1 check_k = MultiExpU1(s, [&v](uint64_t j) -> Fr const& { return v[j]; });

  if (check_k == claim.ki) {
    return false;
  }

  return true;
}
}  // namespace scheme::ot_complaint
