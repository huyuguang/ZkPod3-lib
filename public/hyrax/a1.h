#pragma once

#include "./details.h"

// x, y, z: secret Fr
// open: com(gx,x), com(gy,y), com(gx,z), gx can equal gy
// prove: z = x*y
// proof size: 3 G1 and 5 Fr
// prove cost: 6 eccmul
// verify cost: 6 eccmul
namespace hyrax {
struct A1 {
  struct ProveInput {
    ProveInput(Fr const& x, Fr const& y, Fr const& z, G1 const& gx,
               G1 const& gy)
        : x(x), y(y), z(z), gx(gx), gy(gy) {
      assert(z == x * y);
    }
    Fr const x;
    Fr const y;
    Fr const z;
    G1 const gx;
    G1 const gy;
  };

  struct CommitmentPub {
    CommitmentPub() {}
    CommitmentPub(G1 const& x, G1 const& y, G1 const& z) : x(x), y(y), z(z) {}
    G1 x;  // com(x,r_x)
    G1 y;  // com(y,r_y)
    G1 z;  // com(z,r_z)
  };

  struct CommitmentSec {
    CommitmentSec() {}
    CommitmentSec(Fr const& r_x, Fr const& r_y, Fr const& r_z)
        : r_x(r_x), r_y(r_y), r_z(r_z) {}
    Fr r_x;
    Fr r_y;
    Fr r_z;
  };

  struct CommitmentExtPub {
    CommitmentExtPub() {}
    CommitmentExtPub(G1 const& alpha, G1 const& beta, G1 const& delta)
        : alpha(alpha), beta(beta), delta(delta) {}
    bool operator==(CommitmentExtPub const& right) const {
      return alpha == right.alpha && beta == right.beta && delta == right.delta;
    }
    bool operator!=(CommitmentExtPub const& right) const {
      return !(*this == right);
    }
    G1 alpha;  // com(b1, b2)
    G1 beta;   // com(b3, b4)
    G1 delta;  // com_pub.x * b3 + h * b5
  };

  struct CommitmentExtSec {
    Fr b1;
    Fr b2;
    Fr b3;
    Fr b4;
    Fr b5;
  };

  struct SubProof {
    Fr z1;
    Fr z2;
    Fr z3;
    Fr z4;
    Fr z5;
    bool operator==(SubProof const& right) const {
      return z1 == right.z1 && z2 == right.z2 && z3 == right.z3 &&
             z4 == right.z4 && z5 == right.z5;
    }
    bool operator!=(SubProof const& right) const { return !(*this == right); }
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;  // 3 G1
    SubProof sub_proof;            // 5 Fr
    bool operator==(Proof const& right) const {
      return com_ext_pub == right.com_ext_pub && sub_proof == right.sub_proof;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }
  };

  struct VerifyInput {
    VerifyInput(CommitmentPub const& com_pub, G1 const& gx, G1 const& gy)
        : com_pub(com_pub), gx(gx), gy(gy) {}
    CommitmentPub const& com_pub;
    G1 const gx;
    G1 const gy;
  };

  static bool VerifyInternal(VerifyInput const& input, Fr const& c,
                             CommitmentExtPub const& com_ext_pub,
                             SubProof const& sub_proof) {
    // Tick tick(__FN__);
    auto const& com_pub = input.com_pub;

    std::array<parallel::VoidTask, 3> tasks;
    bool ret0 = false;
    tasks[0] = [&ret0, &com_pub, &com_ext_pub, &c, &sub_proof, &input]() {
      G1 left = com_ext_pub.alpha + com_pub.x * c;
      G1 right = input.gx * sub_proof.z1 + pc::PcH() * sub_proof.z2;
      ret0 = left == right;
    };

    bool ret1 = false;
    tasks[1] = [&ret1, &com_pub, &com_ext_pub, &c, &sub_proof, &input]() {
      G1 left = com_ext_pub.beta + com_pub.y * c;
      G1 right = input.gy * sub_proof.z3 + pc::PcH() * sub_proof.z4;
      ret1 = left == right;
    };

    bool ret2 = false;
    tasks[2] = [&ret2, &com_pub, &com_ext_pub, &c, &sub_proof]() {
      G1 left = com_ext_pub.delta + com_pub.z * c;
      G1 right = com_pub.x * sub_proof.z3 + pc::PcH() * sub_proof.z5;
      ret2 = left == right;
    };

    parallel::Invoke(tasks, true);

    assert(ret0 && ret1 && ret2);
    return ret0 && ret1 && ret2;
  }

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                         ProveInput const& input) {
    // Tick tick(__FN__);
    auto const& gx = input.gx;
    auto const& gy = input.gy;
    auto const& gz = input.gx;
    auto const& h = pc::PcH();

    com_sec.r_x = FrRand();
    com_sec.r_y = FrRand();
    com_sec.r_z = FrRand();

    std::array<parallel::VoidTask, 3> tasks;
    tasks[0] = [&com_pub, &input, &com_sec, &gx, &h]() {
      com_pub.x = gx * input.x + h * com_sec.r_x;
    };
    tasks[1] = [&com_pub, &input, &com_sec, &gy, &h]() {
      com_pub.y = gy * input.y + h * com_sec.r_y;
    };
    tasks[2] = [&com_pub, &input, &com_sec, &gz, &h]() {
      com_pub.z = gz * input.z + h * com_sec.r_z;
    };

    parallel::Invoke(tasks, true);
  }

  static void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                   CommitmentExtSec& com_ext_sec,
                                   ProveInput const& input,
                                   CommitmentPub const& com_pub) {
    // Tick tick(__FN__);
    auto const& gx = input.gx;
    auto const& gy = input.gy;
    auto const& h = pc::PcH();

    com_ext_sec.b1 = FrRand();
    com_ext_sec.b2 = FrRand();
    com_ext_sec.b3 = FrRand();
    com_ext_sec.b4 = FrRand();
    com_ext_sec.b5 = FrRand();

    std::array<parallel::VoidTask, 3> tasks;
    tasks[0] = [&com_ext_pub, &com_ext_sec, &gx, &h]() {
      com_ext_pub.alpha = gx * com_ext_sec.b1 + h * com_ext_sec.b2;
    };
    tasks[1] = [&com_ext_pub, &com_ext_sec, &gy, &h]() {
      com_ext_pub.beta = gy * com_ext_sec.b3 + h * com_ext_sec.b4;
    };
    tasks[2] = [&com_ext_pub, &com_pub, &com_ext_sec, &h]() {
      com_ext_pub.delta = com_pub.x * com_ext_sec.b3 + h * com_ext_sec.b5;
    };
    parallel::Invoke(tasks, true);
  }

  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                         CommitmentExtPub const& com_ext_pub) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_pub.x);
    HashUpdate(hash, com_pub.y);
    HashUpdate(hash, com_pub.z);
    HashUpdate(hash, com_ext_pub.alpha);
    HashUpdate(hash, com_ext_pub.beta);
    HashUpdate(hash, com_ext_pub.delta);
    hash.Final(seed.data());
  }

  static void ComputeSubProof(SubProof& sub_proof, ProveInput const& input,
                              CommitmentSec const& com_sec,
                              CommitmentExtSec const& com_ext_sec,
                              Fr const& c) {
    sub_proof.z1 = com_ext_sec.b1 + c * input.x;
    sub_proof.z2 = com_ext_sec.b2 + c * com_sec.r_x;
    sub_proof.z3 = com_ext_sec.b3 + c * input.y;
    sub_proof.z4 = com_ext_sec.b4 + c * com_sec.r_y;
    sub_proof.z5 = com_ext_sec.b5 + c * (com_sec.r_z - com_sec.r_x * input.y);
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput input,
                    CommitmentPub com_pub, CommitmentSec com_sec) {
    // Tick tick(__FN__);

    CommitmentExtSec com_ext_sec;
    ComputeCommitmentExt(proof.com_ext_pub, com_ext_sec, input, com_pub);

    UpdateSeed(seed, com_pub, proof.com_ext_pub);
    Fr c = H256ToFr(seed);

    ComputeSubProof(proof.sub_proof, input, com_sec, com_ext_sec, c);
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    // Tick tick(__FN__);
    UpdateSeed(seed, input.com_pub, proof.com_ext_pub);
    Fr challenge = H256ToFr(seed);
    return VerifyInternal(input, challenge, proof.com_ext_pub, proof.sub_proof);
  }

  static bool Test();
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, A1::CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("a1.cep", ("a", t.alpha), ("b", t.beta), ("d", t.delta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A1::CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("a1.cep", ("a", t.alpha), ("b", t.beta), ("d", t.delta));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A1::SubProof const& t) {
  ar& YAS_OBJECT_NVP("a1.sp", ("z1", t.z1), ("z2", t.z2), ("z3", t.z3),
                     ("z4", t.z4), ("z5", t.z5));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A1::SubProof& t) {
  ar& YAS_OBJECT_NVP("a1.sp", ("z1", t.z1), ("z2", t.z2), ("z3", t.z3),
                     ("z4", t.z4), ("z5", t.z5));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A1::Proof const& t) {
  ar& YAS_OBJECT_NVP("a1.pf", ("cep", t.com_ext_pub), ("sp", t.sub_proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A1::Proof& t) {
  ar& YAS_OBJECT_NVP("a1.pf", ("cep", t.com_ext_pub), ("sp", t.sub_proof));
}

inline bool A1::Test() {
  Tick tick(__FN__);
  h256_t UpdateSeed = misc::RandH256();

  int64_t x_g_offset = 1;
  int64_t y_g_offset = 10;
  G1 const& gx = pc::PcG(x_g_offset);
  G1 const& gy = pc::PcG(y_g_offset);

  Fr x = FrRand();
  Fr y = FrRand();
  Fr z = x * y;
  ProveInput prove_input(x, y, z, gx, gy);

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prove_input);

  Proof proof;
  Prove(proof, UpdateSeed, prove_input, com_pub, com_sec);

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

  VerifyInput verify_input(com_pub, gx, gy);
  bool success = Verify(proof, UpdateSeed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace hyrax
