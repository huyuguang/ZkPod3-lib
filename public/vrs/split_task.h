#pragma once

#include "./misc.h"

namespace vrs {
template <typename Scheme>
std::vector<std::pair<int64_t, int64_t>> SplitLargeTask(int64_t count) {
  auto const kMaxUnitPerZkp = Scheme::kMaxUnitPerZkp;
  std::vector<std::pair<int64_t, int64_t>> items((count + kMaxUnitPerZkp - 1) /
                                                 kMaxUnitPerZkp);
  for (int64_t i = 0; i < (int64_t)items.size(); ++i) {
    auto& item = items[i];
    item.first = i * kMaxUnitPerZkp;
    item.second = item.first + kMaxUnitPerZkp;
    if (item.second > count) {
      item.second = count;
    }
  }
  return items;
}
}  // namespace vrs