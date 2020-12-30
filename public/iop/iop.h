#pragma once

#include "log/tick.h"
#include "snark/snark.h"

namespace iop {

inline Fr get_root_of_unity(const size_t nth) {
  Fr omega(
      "191032190679217139442913928276920700361456519573292863153056420048214621"
      "61904");
  constexpr size_t s = 28;
  for (size_t i = s; i > nth; --i) {
    omega *= omega;
  }

  return omega;
}

inline void _basic_serial_radix2_FFT(std::vector<Fr> &a, const Fr &omega) {
  const size_t n = a.size(), logn = libff::log2(n);
  assert(n == ((size_t)1 << logn));

  /* swapping in place (from Storer's book) */
  for (size_t k = 0; k < n; ++k) {
    const size_t rk = libff::bitreverse(k, logn);
    if (k < rk) std::swap(a[k], a[rk]);
  }

  size_t m = 1;  // invariant: m = 2^{s-1}
  size_t mul_count = 0;
  for (size_t s = 1; s <= logn; ++s) {
    // w_m is 2^s-th root of unity now
    const Fr w_m = FrPower(omega, (n / (2 * m)));  // omega ^ (n / (2 * m));
#ifndef _MSC_VER
    asm volatile("/* pre-inner */");
#endif
    for (size_t k = 0; k < n; k += 2 * m) {
      Fr w = FrOne();
      for (size_t j = 0; j < m; ++j) {
        const Fr t = w * a[k + j + m];
        a[k + j + m] = a[k + j] - t;
        a[k + j] += t;
        w *= w_m;
        mul_count += 2;
      }
    }
#ifndef _MSC_VER
    asm volatile("/* post-inner */");
#endif
    m *= 2;
  }
  std::cout << "mul_count :" << mul_count << "\n";
}

inline bool Test() {
  snark::InitZkp(true);
  constexpr size_t n = 1ULL << 20;
  auto omega_n = get_root_of_unity(20);
  std::vector<Fr> u(n);
  FrRand(u);

  {
    // lagrange point to coeff
    Tick tick("lagrange point to coeff");
    _basic_serial_radix2_FFT(u, omega_n.inverse());
  }

  constexpr size_t m = 1ULL << 25;
  Fr omega_m = get_root_of_unity(25);
  auto v = u;
  v.resize(m, FrZero());

  {
    Tick tick("coeff to point");
    _basic_serial_radix2_FFT(v, omega_m);
  }
  return true;
}
}  // namespace iop