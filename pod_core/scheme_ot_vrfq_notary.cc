#include "scheme_ot_vrfq_notary.h"

#include "chain.h"
#include "mkl_tree.h"
#include "scheme_misc.h"

namespace scheme::table::ot_vrfq {
bool VerifySecret(G1 const& g_exp_r, Secret const& secret) {
  auto const& ecc_pub = GetEccPub();
  return ecc_pub.PowerG1(secret.r) == g_exp_r;
}

bool VerifySecret(Receipt const& receipt, Secret const& secret) {
  return VerifySecret(receipt.g_exp_r, secret);
}
}  // namespace scheme::table::ot_vrfq
