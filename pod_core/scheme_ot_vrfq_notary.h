#pragma once

#include <memory>
#include <string>

#include "basic_types.h"
#include "scheme_ot_vrfq_protocol.h"

namespace scheme::table::ot_vrfq {
bool VerifySecret(G1 const& g_exp_r, Secret const& secret);
bool VerifySecret(Receipt const& receipt, Secret const& secret);
}  // namespace scheme::table::ot_vrfq