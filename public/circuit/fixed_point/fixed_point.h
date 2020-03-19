#pragma once

#include "./abs_gadget.h"
#include "./div_gadget.h"
#include "./div2_gadget.h"
#include "./exp_gadget.h"
#include "./inv_gadget.h"
#include "./ip_gadget.h"
#include "./mul2_gadget.h"
#include "./mul_gadget.h"
#include "./sign_gadget.h"
#include "./type_gadget.h"

// D: bits of integral part of fixed_point rational
// N: bits of fractional of fixed_point rational
// Fr of data = P/2 + data * 2^N

namespace circuit {

namespace fixed_point {
inline bool Test() {
  //return TestMul();
  //return TestIp();
  // return TestSign();
  // return TestAbs();
  return TestDiv();
  // return TestInv();
  // return TestMul2();
  // return TestExp();
  return TestDiv2();
}
}  // namespace

namespace fp = fixed_point;
}  // namespace circuit::fixed_point
