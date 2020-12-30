#pragma once

#include "./protocol2.h"

// protocol 1
// for a given P, the prover proves that it has vectors a, b for which
// P = g*a + h*b and c= <a,b>

namespace bp::p1 {

struct Proof {
  G1 p;
  Fr c;
  p2::Proof p2;
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.p == right.p && left.c == right.c && left.p2 == right.p2;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("bp.p1.sub_proof", ("p", t.p), ("c", t.c), ("p2", t.p2));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("bp.p1.sub_proof", ("p", t.p), ("c", t.c), ("p2", t.p2));
}

void UpdateSeed(h256_t& seed, G1 const& p, Fr const& c, size_t size) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, p);
  HashUpdate(hash, c);
  HashUpdate(hash, size);
  hash.Final(seed.data());
}

inline void Prove(Proof& proof, h256_t seed, std::vector<G1>&& g,
                  std::vector<G1>&& h, std::vector<Fr>&& a, std::vector<Fr>&& b,
                  G1 const& p, Fr const& c) {
  assert(g.size() == h.size());
  assert(g.size() == a.size());
  assert(g.size() == b.size());
  assert(InnerProduct(a, b) == c);
  assert(MultiExpBdlo12(g, a) + MultiExpBdlo12(h, b) == p);

  UpdateSeed(seed, p, c, a.size());
  G1 u = MapToG1(seed.data(), seed.size());

  proof.p = p;
  proof.c = c;

  p2::Prove(proof.p2, seed, p + u * c, u, std::move(g), std::move(h),
            std::move(a), std::move(b), c);
}

bool Verify(h256_t seed, std::vector<G1>&& g, std::vector<G1>&& h,
            Proof const& proof) {
  UpdateSeed(seed, proof.p, proof.c, g.size());
  G1 u = MapToG1(seed.data(), seed.size());
  return p2::Verify(seed, proof.p + u * proof.c, u, std::move(g), std::move(h),
                    proof.p2);
}

inline bool Test(int64_t n) {
  Tick tick(__FN__);
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

  Proof proof;
  Prove(proof, seed, std::move(g), std::move(h), std::move(a), std::move(b), p,
        c);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << Tick::GetIndentString()
            << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  CHECK(proof == proof2, "");
#endif

  bool success = Verify(seed, std::move(g2), std::move(h2), proof);
  std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace bp::p1