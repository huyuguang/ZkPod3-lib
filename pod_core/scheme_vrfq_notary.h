#pragma once

#include <memory>
#include <string>

#include "misc/misc.h"
#include "scheme_vrfq_protocol.h"

namespace scheme::table::vrfq {
bool VerifySecret(G1 const& g_exp_r, Secret const& secret);
bool VerifySecret(Receipt const& receipt, Secret const& secret);
}  // namespace scheme::table::vrfq