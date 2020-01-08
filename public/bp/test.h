#pragma once

#include "./bp.h"

namespace bp {
inline bool Test() {
  size_t count = 1024;  // * 1024;
  std::vector<G1> g(count);
  for (auto& i : g) {
    i = G1Rand();
  }
  std::vector<Fr> f(count);
  for (auto& i : f) {
    i = FrRand();
  }

  Challenge challenge(count, rand(), false);  // NOTE: use rand() for test

  auto get_g = [&g](uint64_t i) {
    if (i < g.size())
      return g[i];
    else
      return G1Zero();
  };

  auto get_f = [&f](uint64_t i) {
    if (i < f.size())
      return f[i];
    else
      return FrZero();
  };

  P1Proof p1_proof = P1Prove(get_g, get_f, count);

  if (!P1Verify(p1_proof, get_g, count)) {
    assert(false);
    return false;
  }

  std::vector<uint8_t> buf(p1_proof.GetBufSize());
  auto size = p1_proof.serialize(buf.data(), buf.size());
  if (size != buf.size()) {
    assert(false);
    return false;
  }

  P1Proof p1_proof_temp;
  size = p1_proof_temp.deserialize(buf.data(), buf.size());
  if (!size) {
    assert(false);
    return false;
  }
  if (p1_proof_temp != p1_proof) {
    assert(false);
    return false;
  }

  return true;
}
}  // namespace bp