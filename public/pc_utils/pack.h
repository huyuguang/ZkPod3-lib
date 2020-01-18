#pragma once

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"

// x: vector<Fr>, (size+252)/253 = n, {0,1}
// y: vector<Fr>, size = n
// open com(gx,x), com(gy,y)
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

struct ProverInput {
  ProverInput(std::vector<Fr> const& x, G1 const& com_x, Fr const& com_x_r,
              int64_t x_g_offset, std::vector<Fr> const& y, G1 const& com_y,
              Fr const& com_y_r, int64_t y_g_offset)
      : x(x),
        com_x(com_x),
        com_x_r(com_x_r),
        x_g_offset(x_g_offset),
        y(y),
        com_y(com_y),
        com_y_r(com_y_r),
        y_g_offset(y_g_offset) {
#ifdef _DEBUG
    assert((xn() + 252) / 253 == yn());
    assert(FrBitsToFrs(x) == y);
    auto check_x = FrsToFrBits(y);
    assert(check_x.size() >= x.size());
    for (auto i = x.size(); i < check_x.size(); ++i) {
      assert(check_x[i] == 0);
    }
    check_x.resize(x.size());
    assert(check_x == x);
#endif
  }
  int64_t xn() const { return (int64_t)x.size(); }
  int64_t yn() const { return (int64_t)y.size(); }
  std::vector<Fr> const& x;
  G1 const& com_x;
  Fr const& com_x_r;
  int64_t const x_g_offset;
  std::vector<Fr> const& y;
  G1 const& com_y;
  Fr const& com_y_r;
  int64_t const y_g_offset;
};

inline void Prove(Proof& proof, h256_t seed, ProverInput const& input) {
  int64_t xn = input.xn();  
  int64_t yn = input.yn();

  details::UpdateSeed(seed, input.com_x, input.com_y, xn);

  // b
  std::vector<Fr> b(yn);
  ComputeFst(seed, "pack:b", b);

  // c
  std::vector<Fr> c(253);
  c[0] = FrOne();
  for (size_t i = 1; i < c.size(); ++i) {
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

  auto z = InnerProduct(input.y, b);

  assert(z == InnerProduct(input.x, a));

  equal_ip::ProverInput eip_input(input.x, a, input.com_x, input.com_x_r,
                                  input.x_g_offset, input.y, b, input.com_y,
                                  input.com_y_r, input.y_g_offset, z);
  equal_ip::Prove(proof, seed, eip_input);
}

struct VerifierInput {
  VerifierInput(int64_t xn, G1 const& com_x, int64_t x_g_offset,
                G1 const& com_y, int64_t y_g_offset)
      : xn(xn),
        yn((xn + 252) / 253),
        com_x(com_x),
        x_g_offset(x_g_offset),
        com_y(com_y),
        y_g_offset(y_g_offset) {}
  int64_t const xn;
  int64_t const yn;
  G1 const& com_x;
  int64_t const x_g_offset;
  G1 const& com_y;
  int64_t const y_g_offset;
};

// NOET: n is x.size()
inline bool Verify(Proof const& proof, h256_t seed, VerifierInput const& input) {
  int64_t xn = input.xn;
  int64_t yn = input.yn;

  details::UpdateSeed(seed, input.com_x, input.com_y, xn);

  // b
  std::vector<Fr> b(yn);
  ComputeFst(seed, "pack:b", b);

  // c
  std::vector<Fr> c(253);
  c[0] = FrOne();
  for (size_t i = 1; i < c.size(); ++i) {
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

  equal_ip::VerifierInput eip_input(a, input.com_x, input.x_g_offset, b,
                                    input.com_y, input.y_g_offset);
  return equal_ip::Verify(seed, proof,eip_input);
}

inline bool Test() {
  int64_t x_g_offset = 2990;
  int64_t y_g_offset = 670;
  auto seed = misc::RandH256();
  int64_t xn = 1000;
  std::vector<Fr> x(xn);
  for (auto& i : x) i = rand() % 2;
  auto y = FrBitsToFrs(x);
  auto com_x_r = FrRand();
  auto com_x = PcComputeCommitmentG(x_g_offset,x, com_x_r);
  auto com_y_r = FrRand();
  auto com_y = PcComputeCommitmentG(y_g_offset,y, com_y_r);

  ProverInput prover_input(x, com_x, com_x_r, x_g_offset,y, com_y, com_y_r,
                           y_g_offset);
  Proof proof;
  Prove(proof, seed, prover_input);

  VerifierInput verifier_input(xn, com_x,x_g_offset, com_y,y_g_offset);
  bool success = Verify(proof, seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace pc_utils::pack