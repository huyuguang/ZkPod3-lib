#pragma once

#include "./protocol2.h"

// protocol 1
// for a given P, the prover proves that it has vectors a, b for which
// P = g*a + h*b and c= <a,b>

namespace bp {

struct Protocol1Proof {
  G1 p;
  Fr c;
  Protocol2Proof p2;
};

inline bool operator==(Protocol1Proof const& left,
                       Protocol1Proof const& right) {
  return left.p == right.p && left.c == right.c && left.p2 == right.p2;
}

inline bool operator!=(Protocol1Proof const& left,
                       Protocol1Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Protocol1Proof const& t) {
  ar& YAS_OBJECT_NVP("bp.p1proof", ("p", t.p), ("c", t.c), ("p2", t.p2));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Protocol1Proof& t) {
  ar& YAS_OBJECT_NVP("bp.p1proof", ("p", t.p), ("c", t.c), ("p2", t.p2));
}

void UpdateSeed(h256_t& seed, G1 const& p, Fr const& c, size_t size) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, p);
  HashUpdate(hash, c);
  HashUpdate(hash, size);
  hash.Final(seed.data());
}

Protocol1Proof Protocol1Prove(h256_t seed, std::vector<G1>&& g,
                              std::vector<G1>&& h, std::vector<Fr>&& a,
                              std::vector<Fr>&& b, G1 const& p, Fr const& c) {
  assert(g.size() == h.size());
  assert(g.size() == a.size());
  assert(g.size() == b.size());
  assert(InnerProduct(a, b) == c);
  assert(MultiExpBdlo12(g, a) + MultiExpBdlo12(h, b) == p);

  UpdateSeed(seed, p, c, a.size());
  G1 u = MapToG1(seed.data(), seed.size());

  Protocol1Proof proof;
  proof.p = p;
  proof.c = c;

  proof.p2 = Protocol2Prove(seed, p + u * c, u, std::move(g), std::move(h),
                            std::move(a), std::move(b), c);
  return proof;
}

bool Protocol1Verify(h256_t seed, std::vector<G1>&& g, std::vector<G1>&& h,
                     Protocol1Proof const& proof) {
  UpdateSeed(seed, proof.p, proof.c, g.size());
  G1 u = MapToG1(seed.data(), seed.size());
  return Protocol2Verify(seed, proof.p + u * proof.c, u, std::move(g),
                         std::move(h), proof.p2);
}

inline bool TestProtocol1(int64_t n) {
  Tick tick(__FUNCTION__);
  h256_t seed = misc::RandH256();
  std::vector<G1> g(n);
  G1Rand(g.data(), n);

  std::vector<G1> h(n);
  G1Rand(h.data(), n);

  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  std::vector<Fr> b(n);
  FrRand(b.data(), n);

  Fr c = InnerProduct(a, b);
  G1 p = MultiExpBdlo12(g, a) + MultiExpBdlo12(h, b);

  auto g2 = g;
  auto h2 = h;
  auto proof = Protocol1Prove(seed, std::move(g), std::move(h), std::move(a),
                              std::move(b), p, c);

  bool success = Protocol1Verify(seed, std::move(g2), std::move(h2), proof);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace bp