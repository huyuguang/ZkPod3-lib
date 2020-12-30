#pragma once

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

#include <gmp.h>
#include <gmpxx.h>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

namespace misc {
inline bool MpzIsUint256(mpz_class const& v) {
  auto count = (mpz_sizeinbase(v.get_mpz_t(), 2) + 7) / 8;
  return count <= 32;
}

inline void MpzToBE(mpz_class const& v, uint8_t* buf, size_t len) {
  auto count = (mpz_sizeinbase(v.get_mpz_t(), 2) + 7) / 8;
  CHECK(count <= len, "");
  mpz_export(buf + (len - count), nullptr, 1, 1, 0, 0, v.get_mpz_t());
  if (len > count) {
    memset(buf, 0, len - count);
  }
}

inline void MpzToLE(mpz_class const& v, uint8_t* buf, size_t len) {
  auto count = (mpz_sizeinbase(v.get_mpz_t(), 2) + 7) / 8;
  CHECK(count <= len, "");
  mpz_export(buf, nullptr, -1, 1, 0, 0, v.get_mpz_t());
  if (len > count) {
    memset(buf + count, 0, len - count);
  }
}

inline mpz_class MpzFromBE(uint8_t const* buf, size_t len) {
  mpz_class z;
  mpz_import(z.get_mpz_t(), len, 1, 1, 0, 0, buf);
  return z;
}

inline mpz_class MpzFromLE(uint8_t const* buf, size_t len) {
  mpz_class z;
  mpz_import(z.get_mpz_t(), len, -1, 1, 0, 0, buf);
  return z;
}

inline bool MpzFromStr(std::string const& s, mpz_class* mpz) {
  if (mpz->set_str(s, 10) != 0) return false;
  return true;
}

inline std::string MpzToStr(mpz_class const& mpz) { return mpz.get_str(); }

}  // namespace misc
