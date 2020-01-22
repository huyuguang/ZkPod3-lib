#pragma once

#include "./pc_utils.h"

namespace pc_utils {

inline bool Test() {
  std::vector<bool> rets;
  rets.push_back(matrix::Test(10,3));
  rets.push_back(equal_ip::Test(10,3));
  rets.push_back(overlap::Test());
  rets.push_back(divide::Test());
  rets.push_back(match::Test());
  rets.push_back(substr::Test());
  rets.push_back(pack::Test());
  rets.push_back(substrpack::Test());
  rets.push_back(matchpack::Test());
  auto success = std::all_of(rets.begin(), rets.end(), [](bool r) { return r; });
  assert(success);
  return success;
}
}  // namespace pc_utils