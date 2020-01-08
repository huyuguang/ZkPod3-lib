#pragma once

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"

// x: vector<Fr>, size = n, n%s = 0
// s[i]: vector<Fr>, size = sn, i = [0~n/sn), sn can be 1.
// open com(x), com(s[i])
// prove x = {s[i]}

namespace pc_utils::divide::details {
inline void UpdateSeed(h256_t& seed, G1 const& com_x, int64_t sn,
                       std::vector<G1> const& com_s) {
  // update seed
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, sn);
  HashUpdate(hash, com_x);
  HashUpdate(hash, com_s);
  hash.Final(seed.data());
}
}  // namespace pc_utils::divide::details

namespace pc_utils::divide {

using Proof = equal_ip::Proof;

inline void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                  G1 const& com_x, Fr const& com_x_r, int64_t sn,
                  std::vector<G1> const& com_s,
                  std::vector<Fr> const& com_s_r) {
  assert(sn > 0);
  assert(x.size() % sn == 0);
  assert(x.size() / sn == com_s.size());
  assert(com_s.size() == com_s_r.size());

  details::UpdateSeed(seed, com_x, sn, com_s);

  // d, e
  std::vector<Fr> d(com_s.size());
  ComputeFst(seed, "divide:d", d);
  std::vector<Fr> e(sn);
  ComputeFst(seed, "divide:e", e);

  // f
  std::vector<Fr> f;
  f.reserve(x.size());
  for (auto i = 0; i < com_s.size(); ++i) {
    auto de = e * d[i];
    f.insert(f.end(), de.begin(), de.end());
  }

  // y
  std::vector<Fr> y(sn, FrZero());
  auto get_s = [&x, sn](int64_t i) {
    auto begin = x.begin() + i * sn;
    auto end = begin + sn;
    return std::vector<Fr>(begin, end);
  };
  for (auto i = 0; i < com_s.size(); ++i) {
    auto ds = get_s(i) * d[i];
    y = y + ds;
  }

  // com_y_r, com_y
  Fr com_y_r = InnerProduct(com_s_r, d);
  G1 com_y = PcComputeCommitment(y, com_y_r);  // multiexp(sn+1)
#ifdef _DEBUG
  G1 check_com_y = MultiExpBdlo12(com_s, d);  // multiexp(n/sn)
  assert(check_com_y == com_y);
#endif

  // prove <x, f> == <y, e>
  auto z = InnerProduct(x, f);
  assert(z == InnerProduct(y, e));
  equal_ip::Prove(proof, seed, x, f, com_x, com_x_r, y, e, com_y, com_y_r, z);
}

inline bool Verify(h256_t seed, int64_t sn, G1 const& com_x,
                   std::vector<G1> const& com_s, Proof const& proof) {
  details::UpdateSeed(seed, com_x, sn, com_s);
  int64_t n = com_s.size() * sn;

  // d, e
  std::vector<Fr> d(com_s.size());
  ComputeFst(seed, "divide:d", d);
  std::vector<Fr> e(sn);
  ComputeFst(seed, "divide:e", e);

  // f
  std::vector<Fr> f;
  f.reserve(n);
  for (auto i = 0; i < com_s.size(); ++i) {
    auto de = e * d[i];
    f.insert(f.end(), de.begin(), de.end());
  }

  // com_y
  G1 com_y = MultiExpBdlo12(com_s, d);  // multiexp(com_s.size())

  return equal_ip::Verify(seed, f, com_x, e, com_y, proof);
}

inline bool Test() {
  auto seed = misc::RandH256();
  std::vector<Fr> x(15);
  FrRand(x);
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitment(x, com_x_r);

  int64_t sn = 3;
  int64_t m = x.size() / sn;
  std::vector<G1> com_s(m);
  std::vector<Fr> com_s_r(m);
  for (auto i = 0; i < m; ++i) {
    auto begin = x.begin() + i * sn;
    auto end = begin + sn;
    std::vector<Fr> s(begin, end);
    com_s_r[i] = FrRand();
    com_s[i] = PcComputeCommitment(s, com_s_r[i]);
  }

  Proof proof;
  Prove(proof, seed, x, com_x, com_x_r, sn, com_s, com_s_r);

  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);

  Proof proof2;
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream> ia(is);
  ia.serialize(proof2);

  assert(proof == proof2);

  return Verify(seed, sn, com_x, com_s, proof2);
}
}  // namespace pc_utils::divide