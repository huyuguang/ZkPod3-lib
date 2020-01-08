#pragma once

#include "./details.h"
#include "./types.h"

namespace bp {
template <typename GET_G, typename GET_F>
P1Committment P1Commit(GET_G const& get_g, GET_F const& get_f, uint64_t count) {
  Tick tick(__FUNCTION__);
  using namespace details;
  P1Committment p1_committment;
  auto g_count = PackGCount(count);
  p1_committment.c = InnerProduct(
      get_f, [&get_f, g_count](size_t i) { return get_f(i + g_count); },
      g_count);
  p1_committment.p = MultiExpBdlo12<G1>(get_g, get_f, g_count * 2);
  return p1_committment;
}

template <typename GET_G, typename GET_F>
P1Proof P1Prove(GET_G const& get_g, GET_F const& get_f, uint64_t count) {
  P1Proof proof;
  proof.committment = P1Commit(get_g, get_f, count);
  Challenge challenge(count, G1ToBin(proof.committment.p), false);
  proof.p2_proof = P2Prove(proof.committment, get_g, get_f, challenge);

  return proof;
}

template <typename GET_G>
bool P1Verify(P1Committment const& p1_committment, GET_G const& get_g,
              P2Proof const& p2_proof, Challenge const& challenge) {
  if (!P2Verify(p1_committment, get_g, p2_proof, challenge)) return false;
  return p2_proof.q == p1_committment.q(challenge.u());
}

template <typename GET_G>
bool P1Verify(P1Proof const& p1_proof, GET_G const& get_g, uint64_t count) {
  Challenge challenge(count, G1ToBin(p1_proof.committment.p), false);
  if (!P2Verify(p1_proof.committment, get_g, p1_proof.p2_proof, challenge))
    return false;
  return p1_proof.p2_proof.q == p1_proof.committment.q(challenge.u());
}
}  // namespace bp