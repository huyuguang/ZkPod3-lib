#pragma once

#include "./details.h"

// proof of opening
// x: secret vector<Fr>
// open com(x)
// prove know x

namespace clink {
struct Opening {
  struct Proof {
    G1 alpha;
    std::vector<Fr> z1;
    Fr z2;
    bool operator==(Proof const& b) const {
      return alpha == b.alpha && z1 == b.z1 && z2 == b.z2;
    }
    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("o.p", ("a", alpha), ("z1", z1), ("z2", z2));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("o.p", ("a", alpha), ("z1", z1), ("z2", z2));
    }
  };

  static void UpdateSeed(h256_t& seed, G1 const& com_x, G1 const& alpha,
                         int64_t count) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_x);
    HashUpdate(hash, alpha);
    HashUpdate(hash, count);
    hash.Final(seed.data());
  }

  static void Prove(Proof& proof, h256_t seed, std::vector<Fr> const& x,
                    G1 const& com_x,
                    Fr const& r, pc::GetRefG const& get_gx) {
    assert(com_x == pc::PcComputeCommitmentG(get_gx, x, r));
    std::vector<Fr> t1(x.size());
    FrRand(t1);
    Fr t2 = FrRand();
    proof.alpha = pc::PcComputeCommitmentG(get_gx, t1, t2);
    UpdateSeed(seed, com_x, proof.alpha, (int64_t)x.size());
    Fr c = H256ToFr(seed);
    proof.z1 = x * c + t1;
    proof.z2 = r * c + t2;
  }

  static bool Verify(Proof const& proof, h256_t seed, G1 const& com_x,
                     pc::GetRefG const& get_gx) {
    UpdateSeed(seed, com_x, proof.alpha, (int64_t)proof.z1.size());
    Fr c = H256ToFr(seed);
    G1 left = pc::PcComputeCommitmentG(get_gx, proof.z1, proof.z2);
    G1 right = com_x * c + proof.alpha;
    return left == right;
  }

  static bool Test();
};

bool Opening::Test() {
  auto seed = misc::RandH256();
  int64_t g_offset = 30;
  pc::GetRefG get_gx = [g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[g_offset + i];
  };

  std::vector<Fr> x(4);
  FrRand(x);
  Fr r = FrRand();
  G1 com_x = pc::PcComputeCommitmentG(get_gx, x, r);

  Proof proof;
  Prove(proof, seed, x, com_x, r, get_gx);

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  if (proof != proof2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  bool success = Verify(proof, seed, com_x, get_gx);
  std::cout << __FILE__ << " " << __FN__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink