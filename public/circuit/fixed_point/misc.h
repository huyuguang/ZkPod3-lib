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
    bool ret;
    if (neg) {
      ret = abs_mpz > kMpzDN;
    } else {
      ret = abs_mpz >= kMpzDN;
    }

    if (DEBUG_CHECK) {
      if (ret) {
        double double_x = abs_mpz.get_d();
        uint64_t n = N;
        for (;;) {
          if (n > 63) {
            double_x /= (double)(1ULL << 63);
            n -= 63;
          } else {
            double_x /= (double)(1ULL << n);
            break;
          }
        }
        double_x = neg ? -double_x : double_x;
        std::cout << Tick::GetIndentString() << "fixpoint overflow: "
                  << "<" << D << "," << N << ">, double_x=" << double_x << "\n";
      }
    }
    return ret;
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
  auto constances = RationalConst<D, N>();
  bool neg = fr_x.isNegative();
  mpz_class mpz_x = neg ? (-fr_x).getMpz() : fr_x.getMpz();
  CHECK(!constances.IsOverflow(mpz_x, neg), "");

  double double_x = mpz_x.get_d();
  uint64_t n = N;
  for (;;) {
    if (n > 63) {
      double_x /= (double)(1ULL << 63);
      n -= 63;
    } else {
      double_x /= (double)(1ULL << n);
      break;
    }
  }
  return neg ? -double_x : double_x;
}

template <size_t D, size_t N>
Fr DoubleToRational(double double_x) {
  auto constances = RationalConst<D, N>();
  bool neg = double_x < 0;
  double_x = std::abs(double_x);
  uint64_t n = N;
  for (;;) {
    if (n > 63) {
      double_x *= (double)(1ULL << 63);
      n -= 63;
    } else {
      double_x *= (double)(1ULL << n);
      break;
    }
  }
  mpz_class mpz_x = double_x;
  CHECK(!constances.IsOverflow(mpz_x, neg), "");
  Fr fr_x;
  fr_x.setMpz(mpz_x);
  if (neg) fr_x = -fr_x;
  return fr_x;
}

// NOTE
// if a>=0
//  ReducePrecision = (a+2^(D+N))/2^(N-M)-2^M = a/2^(N-M)+2^M-2^M = a/2^(N-M)
// if a<0
//  ReducePrecision = P - ((2^(D+N)-|a|)/2^(N-M)-2^M) = P - |a|/2^(N-M)
// Must same with PrecisionGadget
template <size_t D, size_t N, size_t M>
Fr ReducePrecision(Fr const& a) {
  static_assert(D + N < 253, "invalid D or N");
  static_assert(N > M, "invalid N or M");

  CHECK((!RationalConst<D, N>().IsOverflow(a)), "");

  auto a_off = a + RationalConst<D, N>().kFrDN;

  mpz_class mpz_a_off = a_off.getMpz();
  mpz_a_off /= (mpz_class(1) << (N - M));
  Fr fr_x;
  fr_x.setMpz(mpz_a_off);
  fr_x -= RationalConst<D, M>().kFrDN;

  return fr_x;
}

}  // namespace circuit::fixed_point