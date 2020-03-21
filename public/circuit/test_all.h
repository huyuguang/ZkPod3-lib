#pragma once

#include "./fixed_point/fixed_point.h"
#include "./has0_gadget.h"
#include "./logistic_regression.h"
#include "./match_gadget.h"
#include "./mimc5_gadget.h"
#include "./poseidon_gadget.h"
#include "./sha256c_gadget.h"
#include "./substr_gadget.h"

namespace circuit {
inline bool Test() {
  fixed_point::Test();
  TestFuncH();
  return true;
}
}  // namespace circuit