#pragma once

#include <stdint.h>

#include <algorithm>
#include <functional>
#include <memory>

#include "ecc/ecc.h"
#include "log/tick.h"

namespace vrf::details {
inline G2 GetU() {
  auto& ecc_pub = GetEccPub();
  return ecc_pub.u2()[1];
}
}  // namespace vrf::details