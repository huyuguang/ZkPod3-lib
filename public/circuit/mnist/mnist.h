#pragma once

#include "./convolution.h"

namespace circuit::cnn {
inline bool Test() {
  bool ret;
  std::vector<bool> rets;

  ret = TestConvolution();
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::cnn