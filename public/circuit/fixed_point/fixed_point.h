#pragma once

#include "./abs_gadget.h"
#include "./add_gadget.h"
#include "./div_gadget.h"
#include "./exp_gadget.h"
#include "./ip_gadget.h"
#include "./mul_gadget.h"
#include "./sign_gadget.h"

// D: bits of integral part of fixed_point rational
// N: bits of fractional of fixed_point rational
// Fr of data = P/2 + data * 2^N

namespace circuit {

namespace fixed_point {
inline bool Test() {
  bool ret;
  std::vector<bool> rets;

  ret = TestIp();
  rets.push_back(ret);

  ret = TestSign();
  rets.push_back(ret);

  ret = TestAbs();
  rets.push_back(ret);

  ret = TestInv();
  rets.push_back(ret);

  ret = TestMul();
  rets.push_back(ret);

  ret = TestDiv();
  rets.push_back(ret);

  ret = TestExp();
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}
}  // namespace fixed_point

namespace fp = fixed_point;
}  // namespace circuit
