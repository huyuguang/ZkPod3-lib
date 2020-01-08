#pragma once

#include <memory>
#include <string>

#include "misc/misc.h"
#include "scheme_atomic_swap_vc_protocol.h"

namespace scheme::atomic_swap_vc {
bool VerifySecret(Receipt const& receipt, Secret const& secret);
}  // namespace scheme::atomic_swap_vc
