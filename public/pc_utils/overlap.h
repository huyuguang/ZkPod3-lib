#pragma once

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"

// for vector x and y, open com(x) and com(y), prove overlap

namespace pc_utils::overlap {

using Proof = equal_ip::Proof;

struct OverlapPosition {
  int64_t x;
  int64_t y;
};

inline void UpdateSeed(h256_t& seed, G1 const& c1, G1 const& c2, int64_t xn,
                       int64_t yn, std::vector<OverlapPosition> const& overlap) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, c1);
  HashUpdate(hash, c2);
  HashUpdate(hash, xn);
  HashUpdate(hash, yn);
  for (auto& i : overlap) {
    HashUpdate(hash, i.x);
    HashUpdate(hash, i.y);
  }
  hash.Final(seed.data());
}

inline void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                  std::vector<Fr> const& y,
                  std::vector<OverlapPosition> const& overlap, G1 const& com_x,
                  G1 const& com_y, Fr const& com_x_r, Fr const& com_y_r) {
  int64_t xn = (int64_t)x.size();
  int64_t yn = (int64_t)y.size();

  UpdateSeed(seed, com_x, com_y, xn, yn, overlap);
  std::vector<Fr> c(overlap.size());
  ComputeFst(seed, "consistency::overlap::c", c);

  std::vector<Fr> a(xn, FrZero());
  std::vector<Fr> b(yn, FrZero());
  for (auto i = 0; i < overlap.size(); ++i) {
    auto const& o = overlap[i];
    assert(x[o.x] == y[o.y]);
    a[o.x] += c[i];  // NOTE: here is +=
    b[o.y] += c[i];
  }

  Fr z = InnerProduct(x, a);
  assert(z == InnerProduct(y, b));
  equal_ip::Prove(proof, seed, x, a, com_x, com_x_r, y, b, com_y, com_y_r, z);
}

inline bool Verify(h256_t seed, int64_t xn, G1 const& com_x, int64_t yn,
                   G1 const& com_y, std::vector<OverlapPosition> const& overlap,
                   Proof const& proof) {
  UpdateSeed(seed, com_x, com_y, xn, yn, overlap);
  std::vector<Fr> c(overlap.size());
  ComputeFst(seed, "consistency::overlap::c", c);

  std::vector<Fr> a(xn, FrZero());
  std::vector<Fr> b(yn, FrZero());
  for (auto i = 0; i < overlap.size(); ++i) {
    auto const& o = overlap[i];
    a[o.x] += c[i];  // NOTE: here is +=
    b[o.y] += c[i];
  }

  return equal_ip::Verify(seed, a, com_x, b, com_y, proof);
}

inline bool Test() {
  auto seed = misc::RandH256();
  std::vector<Fr> x(10);
  FrRand(x);
  std::vector<Fr> y(7);
  FrRand(y);

  y[0] = x[1];
  y[1] = x[3];
  y[2] = x[9];
  y[6] = x[9];
  y[4] = x[7];

  std::vector<OverlapPosition> overlap(5);
  overlap[0].x = 1;
  overlap[0].y = 0;
  overlap[1].x = 3;
  overlap[1].y = 1;
  overlap[2].x = 9;
  overlap[2].y = 2;
  overlap[3].x = 9;
  overlap[3].y = 6;
  overlap[4].x = 7;
  overlap[4].y = 4;

  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitment(x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitment(y, com_y_r);

  Proof proof;
  Prove(proof, seed, x, y, overlap, com_x, com_y, com_x_r, com_y_r);

  auto xn = (int64_t)x.size();
  auto yn = (int64_t)y.size();
  if (!Verify(seed, xn, com_x, yn, com_y, overlap, proof)) {
    assert(false);
    return false;
  }
  return true;
}
}  // namespace pc_utils::overlap