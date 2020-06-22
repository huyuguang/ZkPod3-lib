#pragma once

#include "./fixed_point/fixed_point.h"
#include "./has0_gadget.h"
#include "./match_gadget.h"
#include "./mimc5_gadget.h"
#include "./poseidon_gadget.h"
#include "./sha256c_gadget.h"
#include "./substr_gadget.h"
#include "./mnist/mnist.h"

namespace circuit {
inline bool Test() {
  fixed_point::Test();
  //cnn::Test();
  return true;
}
}  // namespace circuit