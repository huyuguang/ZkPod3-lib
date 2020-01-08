#pragma once

#include "./vrf.h"

namespace vrf {
inline bool Test() {
  auto& ecc_pub = GetEccPub();
  Pk<> pk;
  Sk<> sk;
  Generate<>(pk, sk);
  std::array<uint8_t, 32> x;
  for (auto& i : x) {
    i = (uint8_t)rand();  // NOTE: use rand() for test
  }

  Fsk fsk = Vrf<>(sk, x.data());

  Psk<> psk;
  Prove<>(sk, x.data(), psk);

  bool ret = Verify<>(pk, x.data(), fsk, psk);
  assert(ret);
  if (!ret) return false;

  Fr r = FrRand();
  Psk<> psk_exp_r;
  ProveWithR<>(sk, x.data(), r, psk_exp_r);

  G1 g1_exp_r = ecc_pub.PowerG1(r);
  ret = VerifyWithR<>(pk, x.data(), psk_exp_r, g1_exp_r);
  assert(ret);
  if (!ret) return false;

  Fsk fsk2;
  GetFskFromPskExpR(psk_exp_r.back(), r, fsk2);

  assert(fsk2 == fsk);
  if (fsk2 != fsk) return false;
  return true;
}
}  // namespace vrf