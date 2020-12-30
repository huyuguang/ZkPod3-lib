#pragma once

#include "./details.h"

// proof of equality
// x: secret vector<Fr>
// y: secret vector<Fr>
// open com(g, x), com(g, y)
// prove x==y
// NOTE: same g

namespace clink {
struct Equality {
  struct Proof {
    G1 alpha;
    Fr z;
    bool operator==(Proof const& b) const {
      return alpha == b.alpha && z == b.z;
    }
    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("e.p", ("a", alpha), ("z", z));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("e.p", ("a", alpha), ("z", z));
    }
  };

  static void UpdateSeed(h256_t& seed, G1 const& com_x, G1 const& com_y,
                         G1 const& alpha, int64_t count) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_x);
    HashUpdate(hash, com_y);
    HashUpdate(hash, alpha);
    HashUpdate(hash, count);
    hash.Final(seed.data());
  }

  static void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                    G1 const& com_x, Fr const& rx, G1 const& com_y,
                    Fr const& ry, GetRefG1 const& get_gx) {
    (void)get_gx;
    assert(com_x == pc::ComputeCom(get_gx, x, rx));
    assert(com_y == pc::ComputeCom(get_gx, x, ry));
    Fr r = FrRand();
    proof.alpha = pc::PcH() * r;
    UpdateSeed(seed, com_x, com_y, proof.alpha, (int64_t)x.size());
    Fr c = H256ToFr(seed);
    proof.z = c * (rx - ry) + r;
  }

  static bool Verify(Proof const& proof, h256_t seed, G1 const& com_x,
                     G1 const& com_y, int64_t count) {
    UpdateSeed(seed, com_x, com_y, proof.alpha, count);
    Fr c = H256ToFr(seed);
    G1 left = pc::PcH() * proof.z;
    G1 right = (com_x - com_y) * c + proof.alpha;
    return left == right;
  }

  static bool Test();
};

bool Equality::Test() {
  auto seed = misc::RandH256();
  int64_t g_offset = 30;
  GetRefG1 get_gx = [g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[g_offset + i];
  };

  std::vector<Fr> x = {FrRand(), FrRand(), FrRand()};
  Fr rx = FrRand();
  Fr ry = FrRand();
  G1 com_x = pc::ComputeCom(get_gx, x, rx);
  G1 com_y = pc::ComputeCom(get_gx, x, ry);

  Proof proof;
  Prove(proof, seed, x, com_x, rx, com_y, ry, get_gx);

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

  bool success = Verify(proof, seed, com_x, com_y, (int64_t)x.size());
  std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink