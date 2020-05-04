#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>
#include <mutex>
#include <thread>
#include <vector>

#include "ecc/ecc.h"

namespace circuit::fixed_point {
template <size_t D, size_t N>
class RationalConst {
 public:
  RationalConst() {
    std::call_once(once_flag, []() {
      kMpzD3N = mpz_class(1) << (D + 3 * N);
      kFrD3N.setMpz(kMpzD3N);

      kMpzD2N = mpz_class(1) << (D + 2 * N);
      kFrD2N.setMpz(kMpzD2N);

      kMpzDN = mpz_class(1) << (D + N);
      kFrDN.setMpz(kMpzDN);

      kMpzN = mpz_class(1) << N;
      kFrN.setMpz(kMpzN);

      kMpz2N = mpz_class(1) << 2 * N;
      kFr2N.setMpz(kMpz2N);

      kMpz3N = mpz_class(1) << 3 * N;
      kFr3N.setMpz(kMpz3N);
    });
  }

  Fr GetFrDxN(size_t x) const {
    mpz_class mpz_dxn = mpz_class(1) << (D + x * N);
    return SignedMpzToFr(mpz_dxn);
  }

  bool IsOverflow(Fr const& fr) const {
    bool neg = fr.isNegative();
    mpz_class abs_mpz = neg ? (-fr).getMpz() : fr.getMpz();
    return IsOverflow(abs_mpz, neg);
  }

  bool IsOverflow(mpz_class const& abs_mpz, bool neg) const {
    if (neg) {
      return abs_mpz > kMpzDN;
    } else {
      return abs_mpz >= kMpzDN;
    }
  }

  inline static mpz_class kMpzD3N;
  inline static Fr kFrD3N;

  inline static mpz_class kMpzD2N;
  inline static Fr kFrD2N;

  inline static mpz_class kMpzDN;
  inline static Fr kFrDN;

  inline static mpz_class kMpzN;
  inline static Fr kFrN;

  inline static mpz_class kMpz2N;
  inline static Fr kFr2N;

  inline static mpz_class kMpz3N;
  inline static Fr kFr3N;

 private:
  inline static std::once_flag once_flag;
};

template <size_t D, size_t N>
inline double RationalToDouble(Fr const& fr_x) {
  static_assert(N < 64, "N too large");
  auto constances = RationalConst<D, N>();
  bool neg = fr_x.isNegative();
  mpz_class mpz_x = neg ? (-fr_x).getMpz() : fr_x.getMpz();
  if (constances.IsOverflow(mpz_x, neg)) {
    throw std::runtime_error(__FN__);
  }
  double double_x = mpz_x.get_d();
  double_x /= (double)(1ULL << N);
  return neg ? -double_x : double_x;
}

template <size_t D, size_t N>
inline Fr DoubleToRational(double double_x) {
  static_assert(N < 64, "N too large");
  auto constances = RationalConst<D, N>();
  bool neg = double_x < 0;
  mpz_class mpz_x = std::abs(double_x) * (1ULL << N);
  if (constances.IsOverflow(mpz_x, neg)) {
    throw std::runtime_error(__FN__);
  }
  Fr fr_x;
  fr_x.setMpz(mpz_x);
  if (neg) fr_x = -fr_x;
  return fr_x;
}
}  // namespace circuit::fixed_point