#pragma once

#include <stdint.h>

#include <algorithm>
#include <boost/endian/conversion.hpp>
#include <functional>
#include <memory>

#include "ecc/ecc.h"
#include "misc/misc.h"

namespace bp::details {
inline size_t PackGCount(size_t count) {
  assert(count);
  auto align_count = misc::Pow2UB(count);
  if (align_count == 1) return 1;
  return align_count / 2;
}

inline size_t PackXCount(size_t count) {
  auto g_count = PackGCount(count);
  return misc::Log2UB(g_count);
}

inline G1 MultiExpGH(G1 const* g, Fr const* a, G1 const* h, Fr const* b,
                     size_t n) {
  auto get_g = [g, h, n](size_t i) -> G1 const& {
    return i < n ? g[i] : h[i - n];
  };
  auto get_f = [a, b, n](size_t i) -> Fr const& {
    return i < n ? a[i] : b[i - n];
  };
  return MultiExpBdlo12<G1>(get_g, get_f, n * 2, true);
}
}  // namespace bp::details
