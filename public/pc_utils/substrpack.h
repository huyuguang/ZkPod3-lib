#pragma once

#include "./pack.h"
#include "./substr.h"

namespace pc_utils::substrpack {

struct Proof {
  substr::Proof substr_proof;
  pack::Proof pack_proof;
  G1 com_pack_y;
};

inline bool operator==(Proof const& a, Proof const& b) {
  return a.substr_proof == b.substr_proof && a.pack_proof == b.pack_proof &&
    a.com_pack_y == b.com_pack_y;
}

inline bool operator!=(Proof const& a, Proof const& b) { return !(a == b); }

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("sp.p", ("s", t.substr_proof), ("p", t.pack_proof),
                     ("y", t.com_pack_y));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("sp.p", ("s", t.substr_proof), ("p", t.pack_proof),
                     ("y", t.com_pack_y));
}

struct ProveOutput {
  substr::Proof substr_proof;
  pack::Proof pack_proof;
  std::vector<Fr> y;
  Fr com_y_r;
  G1 com_pack_y;
  std::vector<Fr> pack_y;  
  Fr com_pack_y_r;
  Proof BuildProof() const {
    Proof ret;
    ret.substr_proof = substr_proof;
    ret.pack_proof = pack_proof;
    ret.com_pack_y = com_pack_y;
    return ret;
  }
};

inline void Prove(ProveOutput& output, h256_t seed, std::vector<Fr> const& x,
                  std::string const& k, G1 const& com_x, Fr const& com_x_r) {
  // prove y[i]=substr(x[i],k), i=[0,s)
  int64_t n = (int64_t)x.size();
  output.y.resize(x.size());
  std::vector<std::string> x_str(n);
  for (int64_t i = 0; i < n; ++i) {
    auto str = UnPackStrFromFr(x[i]);
    output.y[i] = str.find(k) != std::string::npos ? 1 : 0;
  }

  output.com_y_r = FrRand();
  auto com_y = PcComputeCommitment(output.y, output.com_y_r, true);

  auto& substr_proof = output.substr_proof;
  pc_utils::substr::Prove(substr_proof, seed, x, k, com_x, com_x_r, com_y,
                          output.com_y_r);
  assert(substr_proof.com_w.back() == com_y);

  // pack y to pack_y
  output.pack_y = FrBitsToFrs(output.y);
  output.com_pack_y_r = FrRand();
  output.com_pack_y = PcComputeCommitment(output.pack_y, output.com_pack_y_r);
  auto& pack_proof = output.pack_proof;
  pc_utils::pack::Prove(pack_proof, seed, output.y, com_y, output.com_y_r,
                        output.pack_y, output.com_pack_y, output.com_pack_y_r);
}

inline bool Verify(Proof const& proof, h256_t seed, int64_t n,
                   std::string const& k) {
  if (!substr::Verify(proof.substr_proof, seed, n, k)) return false;
  G1 const& com_y = proof.substr_proof.com_w.back();
  return pack::Verify(proof.pack_proof, seed, n, com_y, proof.com_pack_y);
}

inline bool Test() {
  auto seed = misc::RandH256();
  int64_t n = 1000;
  std::vector<Fr> x(n);
  FrRand(x);
  x.push_back(PackStrToFr("abcdefg"));
  auto com_x_r = FrRand();
  auto com_x = PcComputeCommitment(x, com_x_r);
  std::string k = "cde";
  ProveOutput output;
  Prove(output, seed, x, k, com_x, com_x_r);

  auto proof = output.BuildProof();
  if (proof.substr_proof.com_x() != com_x) {
    assert(false);
    return false;
  }
  if (!Verify(proof, seed, n, k)) {
    assert(false);
    return false;
  }

  std::cout << __FUNCTION__ << ": success\n";
  return true;
}

}  // namespace pc_utils::substrpack