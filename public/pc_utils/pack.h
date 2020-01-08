#pragma once

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"

// x: vector<Fr>, (size+252)/253 = n, {0,1}
// y: vector<Fr>, size = n
// open com(x), com(y)
// prove binary of y == x
// let b: vector<Fr>, size = n, b = fst(com(x), com(y))
// let c: vector<Fr>, size = 253, a = {2^0, 2^1.... 2^252}
// let a = {{c*b0}, {c*b1}, {c*b2}.... {c*b252}}, a.size = 253*n
// prove <x, a> == <y, b>

namespace pc_utils::pack::details {
inline void UpdateSeed(h256_t& seed, G1 const& c1, G1 const& c2, int64_t n) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, c1);
  HashUpdate(hash, c2);
  HashUpdate(hash, n);
  hash.Final(seed.data());
}
}  // namespace pc_utils::pack::details

namespace pc_utils::pack {

using Proof = equal_ip::Proof;

inline void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                  G1 const& com_x, Fr const& com_x_r, std::vector<Fr> const& y,
                  G1 const& com_y, Fr const& com_y_r) {
#ifdef _DEBUG
  assert(FrBitsToFrs(x) == y);
  auto check_x = FrsToFrBits(y);
  assert(check_x.size() >= x.size());
  for (auto i = x.size(); i < check_x.size(); ++i) {
    assert(check_x[i] == 0);
  }
  check_x.resize(x.size());
  assert(check_x == x);
#endif

  int64_t xn = (int64_t)x.size();
  int64_t yn = (xn + 252) / 253;
  assert(yn == (int64_t)y.size());

  details::UpdateSeed(seed, com_x, com_y, xn);

  // b
  std::vector<Fr> b(yn);
  ComputeFst(seed, "pack:b", b);

  // c
  std::vector<Fr> c(253);
  c[0] = FrOne();
  for (auto i = 1; i < c.size(); ++i) {
    c[i] = c[i - 1] + c[i - 1];
  }

  // a
  std::vector<Fr> a;
  a.reserve(yn * 253);
  for (int64_t i = 0; i < yn; ++i) {
    auto cb = c * b[i];
    a.insert(a.end(), cb.begin(), cb.end());
  }
  a.resize(xn);

  auto z = InnerProduct(y, b);

  assert(z == InnerProduct(x, a));

  equal_ip::Prove(proof, seed, x, a, com_x, com_x_r, y, b, com_y, com_y_r, z);
}

// NOET: n is x.size()
inline bool Verify(Proof const& proof, h256_t seed, int64_t xn, G1 const& com_x,
                   G1 const& com_y) {
  int64_t yn = (xn + 252) / 253;

  details::UpdateSeed(seed, com_x, com_y, xn);

  // b
  std::vector<Fr> b(yn);
  ComputeFst(seed, "pack:b", b);

  // c
  std::vector<Fr> c(253);
  c[0] = FrOne();
  for (auto i = 1; i < c.size(); ++i) {
    c[i] = c[i - 1] + c[i - 1];
  }

  // a
  std::vector<Fr> a;
  a.reserve(yn * 253);
  for (int64_t i = 0; i < yn; ++i) {
    auto cb = c * b[i];
    a.insert(a.end(), cb.begin(), cb.end());
  }
  a.resize(xn);

  return equal_ip::Verify(seed, a, com_x, b, com_y, proof);
}

inline bool Test() {
  auto seed = misc::RandH256();
  int64_t xn = 1000;
  std::vector<Fr> x(xn);
  for (auto& i : x) i = rand() % 2;
  auto y = FrBitsToFrs(x);
  auto com_x_r = FrRand();
  auto com_x = PcComputeCommitment(x, com_x_r);
  auto com_y_r = FrRand();
  auto com_y = PcComputeCommitment(y, com_y_r);

  Proof proof;
  Prove(proof, seed, x, com_x, com_x_r, y, com_y, com_y_r);

  return Verify(proof, seed, xn, com_x, com_y);
}
}  // namespace pc_utils::pack