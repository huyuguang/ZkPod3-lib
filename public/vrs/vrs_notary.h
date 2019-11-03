#pragma once

#include "ecc.h"
#include "ecc_pub.h"
#include "vrs_verifier.h"

namespace vrs {
inline bool VerifySecret(G1 const& h, G1 const& g, G1 const& key_com,
                        Fr const& r, Fr const& key) {
  return key_com == h * r + g * key;
}
}