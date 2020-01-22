#pragma once

#include <cryptopp/keccak.h>

#include <vector>

#include "ecc/ecc.h"

namespace vrs {

namespace details {
inline uint32_t ReadBE32(const uint8_t* ptr) {
  return ((uint32_t) * ((ptr) + 3)) | ((uint32_t) * ((ptr) + 2) << 8) |
         ((uint32_t) * ((ptr) + 1) << 16) | ((uint32_t) * ((ptr) + 0) << 24);
}

inline void WriteBE32(uint8_t* ptr, uint32_t x) {
  *((ptr) + 3) = (uint8_t)((x));
  *((ptr) + 2) = (uint8_t)((x) >> 8);
  *((ptr) + 1) = (uint8_t)((x) >> 16);
  *((ptr) + 0) = (uint8_t)((x) >> 24);
}

inline uint32_t Ch(uint32_t x, uint32_t y, uint32_t z) {
  return z ^ (x & (y ^ z));
}

inline uint32_t Maj(uint32_t x, uint32_t y, uint32_t z) {
  return (x & y) | (z & (x | y));
}

inline uint32_t Sigma0(uint32_t x) {
  return (x >> 2 | x << 30) ^ (x >> 13 | x << 19) ^ (x >> 22 | x << 10);
}

inline uint32_t Sigma1(uint32_t x) {
  return (x >> 6 | x << 26) ^ (x >> 11 | x << 21) ^ (x >> 25 | x << 7);
}

inline uint32_t sigma0(uint32_t x) {
  return (x >> 7 | x << 25) ^ (x >> 18 | x << 14) ^ (x >> 3);
}

inline uint32_t sigma1(uint32_t x) {
  return (x >> 17 | x << 15) ^ (x >> 19 | x << 13) ^ (x >> 10);
}

/** One round of SHA-256. */
inline void Round(uint32_t a, uint32_t b, uint32_t c, uint32_t& d, uint32_t e,
                  uint32_t f, uint32_t g, uint32_t& h, uint32_t k, uint32_t w) {
  uint32_t t1 = h + Sigma1(e) + Ch(e, f, g) + k + w;
  uint32_t t2 = Sigma0(a) + Maj(a, b, c);
  d += t1;
  h = t1 + t2;
}

/** Initialize SHA-256 state. */
inline void Initialize(uint32_t* s) {
  s[0] = 0x6a09e667ul;
  s[1] = 0xbb67ae85ul;
  s[2] = 0x3c6ef372ul;
  s[3] = 0xa54ff53aul;
  s[4] = 0x510e527ful;
  s[5] = 0x9b05688cul;
  s[6] = 0x1f83d9abul;
  s[7] = 0x5be0cd19ul;
}

/** Perform one SHA-256 transformation, processing a 64-byte chunk. */
inline void Transform(uint32_t* s, const uint8_t* chunk) {
  uint32_t a = s[0], b = s[1], c = s[2], d = s[3], e = s[4], f = s[5], g = s[6],
           h = s[7];
  uint32_t w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;

  Round(a, b, c, d, e, f, g, h, 0x428a2f98, w0 = ReadBE32(chunk + 0));
  Round(h, a, b, c, d, e, f, g, 0x71374491, w1 = ReadBE32(chunk + 4));
  Round(g, h, a, b, c, d, e, f, 0xb5c0fbcf, w2 = ReadBE32(chunk + 8));
  Round(f, g, h, a, b, c, d, e, 0xe9b5dba5, w3 = ReadBE32(chunk + 12));
  Round(e, f, g, h, a, b, c, d, 0x3956c25b, w4 = ReadBE32(chunk + 16));
  Round(d, e, f, g, h, a, b, c, 0x59f111f1, w5 = ReadBE32(chunk + 20));
  Round(c, d, e, f, g, h, a, b, 0x923f82a4, w6 = ReadBE32(chunk + 24));
  Round(b, c, d, e, f, g, h, a, 0xab1c5ed5, w7 = ReadBE32(chunk + 28));
  Round(a, b, c, d, e, f, g, h, 0xd807aa98, w8 = ReadBE32(chunk + 32));
  Round(h, a, b, c, d, e, f, g, 0x12835b01, w9 = ReadBE32(chunk + 36));
  Round(g, h, a, b, c, d, e, f, 0x243185be, w10 = ReadBE32(chunk + 40));
  Round(f, g, h, a, b, c, d, e, 0x550c7dc3, w11 = ReadBE32(chunk + 44));
  Round(e, f, g, h, a, b, c, d, 0x72be5d74, w12 = ReadBE32(chunk + 48));
  Round(d, e, f, g, h, a, b, c, 0x80deb1fe, w13 = ReadBE32(chunk + 52));
  Round(c, d, e, f, g, h, a, b, 0x9bdc06a7, w14 = ReadBE32(chunk + 56));
  Round(b, c, d, e, f, g, h, a, 0xc19bf174, w15 = ReadBE32(chunk + 60));

  Round(a, b, c, d, e, f, g, h, 0xe49b69c1,
        w0 += sigma1(w14) + w9 + sigma0(w1));
  Round(h, a, b, c, d, e, f, g, 0xefbe4786,
        w1 += sigma1(w15) + w10 + sigma0(w2));
  Round(g, h, a, b, c, d, e, f, 0x0fc19dc6,
        w2 += sigma1(w0) + w11 + sigma0(w3));
  Round(f, g, h, a, b, c, d, e, 0x240ca1cc,
        w3 += sigma1(w1) + w12 + sigma0(w4));
  Round(e, f, g, h, a, b, c, d, 0x2de92c6f,
        w4 += sigma1(w2) + w13 + sigma0(w5));
  Round(d, e, f, g, h, a, b, c, 0x4a7484aa,
        w5 += sigma1(w3) + w14 + sigma0(w6));
  Round(c, d, e, f, g, h, a, b, 0x5cb0a9dc,
        w6 += sigma1(w4) + w15 + sigma0(w7));
  Round(b, c, d, e, f, g, h, a, 0x76f988da, w7 += sigma1(w5) + w0 + sigma0(w8));
  Round(a, b, c, d, e, f, g, h, 0x983e5152, w8 += sigma1(w6) + w1 + sigma0(w9));
  Round(h, a, b, c, d, e, f, g, 0xa831c66d,
        w9 += sigma1(w7) + w2 + sigma0(w10));
  Round(g, h, a, b, c, d, e, f, 0xb00327c8,
        w10 += sigma1(w8) + w3 + sigma0(w11));
  Round(f, g, h, a, b, c, d, e, 0xbf597fc7,
        w11 += sigma1(w9) + w4 + sigma0(w12));
  Round(e, f, g, h, a, b, c, d, 0xc6e00bf3,
        w12 += sigma1(w10) + w5 + sigma0(w13));
  Round(d, e, f, g, h, a, b, c, 0xd5a79147,
        w13 += sigma1(w11) + w6 + sigma0(w14));
  Round(c, d, e, f, g, h, a, b, 0x06ca6351,
        w14 += sigma1(w12) + w7 + sigma0(w15));
  Round(b, c, d, e, f, g, h, a, 0x14292967,
        w15 += sigma1(w13) + w8 + sigma0(w0));

  Round(a, b, c, d, e, f, g, h, 0x27b70a85,
        w0 += sigma1(w14) + w9 + sigma0(w1));
  Round(h, a, b, c, d, e, f, g, 0x2e1b2138,
        w1 += sigma1(w15) + w10 + sigma0(w2));
  Round(g, h, a, b, c, d, e, f, 0x4d2c6dfc,
        w2 += sigma1(w0) + w11 + sigma0(w3));
  Round(f, g, h, a, b, c, d, e, 0x53380d13,
        w3 += sigma1(w1) + w12 + sigma0(w4));
  Round(e, f, g, h, a, b, c, d, 0x650a7354,
        w4 += sigma1(w2) + w13 + sigma0(w5));
  Round(d, e, f, g, h, a, b, c, 0x766a0abb,
        w5 += sigma1(w3) + w14 + sigma0(w6));
  Round(c, d, e, f, g, h, a, b, 0x81c2c92e,
        w6 += sigma1(w4) + w15 + sigma0(w7));
  Round(b, c, d, e, f, g, h, a, 0x92722c85, w7 += sigma1(w5) + w0 + sigma0(w8));
  Round(a, b, c, d, e, f, g, h, 0xa2bfe8a1, w8 += sigma1(w6) + w1 + sigma0(w9));
  Round(h, a, b, c, d, e, f, g, 0xa81a664b,
        w9 += sigma1(w7) + w2 + sigma0(w10));
  Round(g, h, a, b, c, d, e, f, 0xc24b8b70,
        w10 += sigma1(w8) + w3 + sigma0(w11));
  Round(f, g, h, a, b, c, d, e, 0xc76c51a3,
        w11 += sigma1(w9) + w4 + sigma0(w12));
  Round(e, f, g, h, a, b, c, d, 0xd192e819,
        w12 += sigma1(w10) + w5 + sigma0(w13));
  Round(d, e, f, g, h, a, b, c, 0xd6990624,
        w13 += sigma1(w11) + w6 + sigma0(w14));
  Round(c, d, e, f, g, h, a, b, 0xf40e3585,
        w14 += sigma1(w12) + w7 + sigma0(w15));
  Round(b, c, d, e, f, g, h, a, 0x106aa070,
        w15 += sigma1(w13) + w8 + sigma0(w0));

  Round(a, b, c, d, e, f, g, h, 0x19a4c116,
        w0 += sigma1(w14) + w9 + sigma0(w1));
  Round(h, a, b, c, d, e, f, g, 0x1e376c08,
        w1 += sigma1(w15) + w10 + sigma0(w2));
  Round(g, h, a, b, c, d, e, f, 0x2748774c,
        w2 += sigma1(w0) + w11 + sigma0(w3));
  Round(f, g, h, a, b, c, d, e, 0x34b0bcb5,
        w3 += sigma1(w1) + w12 + sigma0(w4));
  Round(e, f, g, h, a, b, c, d, 0x391c0cb3,
        w4 += sigma1(w2) + w13 + sigma0(w5));
  Round(d, e, f, g, h, a, b, c, 0x4ed8aa4a,
        w5 += sigma1(w3) + w14 + sigma0(w6));
  Round(c, d, e, f, g, h, a, b, 0x5b9cca4f,
        w6 += sigma1(w4) + w15 + sigma0(w7));
  Round(b, c, d, e, f, g, h, a, 0x682e6ff3, w7 += sigma1(w5) + w0 + sigma0(w8));
  Round(a, b, c, d, e, f, g, h, 0x748f82ee, w8 += sigma1(w6) + w1 + sigma0(w9));
  Round(h, a, b, c, d, e, f, g, 0x78a5636f,
        w9 += sigma1(w7) + w2 + sigma0(w10));
  Round(g, h, a, b, c, d, e, f, 0x84c87814,
        w10 += sigma1(w8) + w3 + sigma0(w11));
  Round(f, g, h, a, b, c, d, e, 0x8cc70208,
        w11 += sigma1(w9) + w4 + sigma0(w12));
  Round(e, f, g, h, a, b, c, d, 0x90befffa,
        w12 += sigma1(w10) + w5 + sigma0(w13));
  Round(d, e, f, g, h, a, b, c, 0xa4506ceb,
        w13 += sigma1(w11) + w6 + sigma0(w14));
  Round(c, d, e, f, g, h, a, b, 0xbef9a3f7,
        w14 + sigma1(w12) + w7 + sigma0(w15));
  Round(b, c, d, e, f, g, h, a, 0xc67178f2,
        w15 + sigma1(w13) + w8 + sigma0(w0));

  s[0] += a;
  s[1] += b;
  s[2] += c;
  s[3] += d;
  s[4] += e;
  s[5] += f;
  s[6] += g;
  s[7] += h;
}

inline void Sha256Compress(const uint8_t data[64], uint8_t hash[32]) {
  uint32_t s[8];
  Initialize(s);
  Transform(s, data);
  WriteBE32(hash, s[0]);
  WriteBE32(hash + 4, s[1]);
  WriteBE32(hash + 8, s[2]);
  WriteBE32(hash + 12, s[3]);
  WriteBE32(hash + 16, s[4]);
  WriteBE32(hash + 20, s[5]);
  WriteBE32(hash + 24, s[6]);
  WriteBE32(hash + 28, s[7]);
}

inline bool test_sha256_compress() {
  uint8_t data[64] = {183, 231, 178, 111, 197, 66,  169, 241, 210, 48,  239,
                      205, 118, 75,  152, 233, 23,  244, 68,  121, 155, 134,
                      181, 131, 32,  157, 253, 177, 49,  186, 62,  132, 0x80,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                      0,   0,   0,   0,   0,   0,   0,   1,   0};

  uint8_t sig[32] = {78,  144, 206, 42,  80,  100, 176, 75,  200, 232, 113,
                     98,  19,  218, 162, 124, 58,  186, 16,  209, 143, 237,
                     155, 247, 76,  51,  189, 234, 207, 145, 110, 196};

  uint8_t out[32];
  Sha256Compress(data, out);
  if (memcmp(out, sig, 32) != 0) {
    assert(false);
    return false;
  }

  // copy from test_sha256_gadget.cpp
  uint8_t data2[64] = {
      0x42, 0x6b, 0xc2, 0xd8,  // 0x426bc2d8
      0x4d, 0xc8, 0x67, 0x82,  // 0x4dc86782
      0x81, 0xe8, 0x95, 0x7a,  // 0x81e8957a
      0x40, 0x9e, 0xc1, 0x48,  // 0x409ec148
      0xe6, 0xcf, 0xfb, 0xe8,  // 0xe6cffbe8
      0xaf, 0xe6, 0xba, 0x4f,  // 0xafe6ba4f
      0x9c, 0x6f, 0x19, 0x78,  // 0x9c6f1978
      0xdd, 0x7a, 0xf7, 0xe9,  // 0xdd7af7e9

      0x03, 0x8c, 0xce, 0x42,  // 0x038cce42
      0xab, 0xd3, 0x66, 0xb8,  // 0xabd366b8
      0x3e, 0xde, 0x7e, 0x00,  // 0x3ede7e00
      0x91, 0x30, 0xde, 0x53,  // 0x9130de53
      0x72, 0xcd, 0xf7, 0x3d,  // 0x72cdf73d
      0xee, 0x82, 0x51, 0x14,  // 0xee825114
      0x8c, 0xb4, 0x8d, 0x1b,  // 0x8cb48d1b
      0x9a, 0xf6, 0x8a, 0xd0   // 0x9af68ad0
  };

  uint8_t sig2[32] = {
      0xef, 0xfd, 0x0b, 0x7f,  // 0xeffd0b7f
      0x1c, 0xcb, 0xa1, 0x16,  // 0x1ccba116
      0x2e, 0xe8, 0x16, 0xf7,  // 0x2ee816f7
      0x31, 0xc6, 0x2b, 0x48,  // 0x31c62b48
      0x59, 0x30, 0x51, 0x41,  // 0x59305141
      0x99, 0x0e, 0x5c, 0x0a,  // 0x990e5c0a
      0xce, 0x40, 0xd3, 0x3d,  // 0xce40d33d
      0x0b, 0x11, 0x67, 0xd1   // 0x0b1167d1
  };

  Sha256Compress(data2, out);
  if (memcmp(out, sig2, 32) != 0) {
    assert(false);
    return false;
  }

  uint8_t data3[64] = {
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,

      0x80, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00,
  };

  uint8_t sig3[32] = {
      0xa3,0xf4,0x65,0x0a,
      0x88,0x4a,0x0d,0x87,
      0x81,0x81,0xf5,0x59,
      0xeb,0xea,0x9b,0x82,
      0x04,0x29,0x94,0x94,
      0x12,0x25,0x8e,0x30,
      0xe2,0xd5,0xea,0x86,
      0xa0,0x1b,0xde,0x57
  };
  Sha256Compress(data3, out);
  if (memcmp(out, sig3, 32) != 0) {
    assert(false);
    return false;
  }
  return true;
}

inline h256_t Sha256Compress(h256_t const& a, h256_t const& b) {
  uint8_t data[64];
  memcpy(data, a.data(), 32);
  memcpy(data + 32, b.data(), 32);
  h256_t c;
  Sha256Compress(data, c.data());
  return c;
}

inline uint8_t ReverseBits(uint8_t num) {
  uint8_t reverse_num = 0;
  for (int i = 0; i < 8; i++) {
    if ((num & (1 << i))) {
      reverse_num |= 1 << ((8 - 1) - i);
    }
  }
  return reverse_num;
}
}  // namespace details

inline Fr Sha256Enc(Fr const& plain, Fr const& key) {
  h256_t p = FrToBin(plain);
  h256_t k = FrToBin(key);

  for (auto& i : p) i = details::ReverseBits(i);
  for (auto& i : k) i = details::ReverseBits(i);

  h256_t v = details::Sha256Compress(p, k);
  for (auto& i : v) i = details::ReverseBits(i);
  return H256ToFr(v);
}
}  // namespace vrs