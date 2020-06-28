#pragma once 

#include "./ip_gadget.h"
#include "./relubn_gadget.h"
#include "./pooling_gadget.h"

namespace circuit::vgg16 {
inline bool Test() {
  //bool ret;
  std::vector<bool> rets;

  // ret = TestConvolution();
  // rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace circuit::vgg16