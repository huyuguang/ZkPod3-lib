#pragma once

#include "./details.h"

// a: public vector<Fr>, size = n
// x: secret vector<Fr>, size = n
// y: secret Fr, = <x,a>
// open: com(gx,x), com(gy,y)
// prove: y=<x,a>
// proof size: 2 G1 and n+2 Fr
// prove cost: mulexp(n)
// verify cost: mulexp(n)
namespace hyrax {
struct A2 {
  struct ProveInput {
    ProveInput(std::string const& tag, std::vector<Fr> const& x,
               std::vector<Fr> const& a, Fr const& y, GetRefG1 const& get_gx,
               G1 const& gy)
        : tag(tag), x(x), a(a), y(y), get_gx(get_gx), gy(gy) {
      assert(x.size() == a.size() && !a.empty());
      assert(y == InnerProduct(x, a));
    }
    int64_t n() const { return (int64_t)x.size(); }
    std::string to_string() const { return tag + ": " + std::to_string(n()); }

    std::string tag;
    std::vector<Fr> const& x;  // x.size = n // TODO: change to move
    std::vector<Fr> const& a;  // a.size = n
    Fr const y;                // y = <x, a>
    GetRefG1 const get_gx;
    G1 const gy;
  };

  struct CommitmentPub {
    CommitmentPub() {}
    CommitmentPub(G1 const& xi, G1 const& tau) : xi(xi), tau(tau) {}
    G1 xi;   // com(gx, x, r_xi)
    G1 tau;  // com(gy, y, r_tau)
    bool operator==(CommitmentPub const& right) const {
      return xi == right.xi && tau == right.tau;
    }

    bool operator!=(CommitmentPub const& right) const {
      return !(*this == right);
    }
  };

  struct CommitmentSec {
    CommitmentSec() {}
    CommitmentSec(Fr const& x, Fr const& t) : r_xi(x), r_tau(t) {}
    Fr r_xi;
    Fr r_tau;
  };

  struct CommitmentExtPub {
    CommitmentExtPub() {}
    CommitmentExtPub(G1 const& delta, G1 const& beta)
        : delta(delta), beta(beta) {}
    G1 delta;  // com(d, r_delta)
    G1 beta;   // com(<a,d>, r_beta)
    bool operator==(CommitmentExtPub const& right) const {
      return delta == right.delta && beta == right.beta;
    }

    bool operator!=(CommitmentExtPub const& right) const {
      return !(*this == right);
    }
  };

  struct CommitmentExtSec {
    std::vector<Fr> d;  // size = n
    Fr r_beta;
    Fr r_delta;
  };

  struct SubProof {
    std::vector<Fr> z;  // z.size = n
    Fr z_delta;
    Fr z_beta;
    int64_t n() const { return (int64_t)z.size(); }
    bool operator==(SubProof const& right) const {
      return z == right.z && z_delta == right.z_delta && z_beta == right.z_beta;
    }

    bool operator!=(SubProof const& right) const { return !(*this == right); }
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;  // 2 G1
    SubProof sub_proof;            // n+2 Fr
    int64_t n() const { return sub_proof.n(); }
    bool operator==(Proof const& right) const {
      return com_ext_pub == right.com_ext_pub && sub_proof == right.sub_proof;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }
  };

  struct VerifyInput {
    VerifyInput(std::string const& tag, std::vector<Fr> const& a,
                CommitmentPub const& com_pub, GetRefG1 const& get_gx,
                G1 const& gy)
        : tag(tag), a(a), com_pub(com_pub), get_gx(get_gx), gy(gy) {}
    std::string tag;
    std::vector<Fr> const& a;  // a.size = n
    CommitmentPub const& com_pub;
    GetRefG1 const get_gx;
    G1 const gy;
    int64_t n() const { return (int64_t)a.size(); }
    std::string to_string() const { return tag + ": " + std::to_string(n()); }
  };

  // com(n) + com(1) + ip(n)
  static bool VerifyInternal(VerifyInput const& input, Fr const& challenge,
                             CommitmentExtPub const& com_ext_pub,
                             SubProof const& sub_proof) {
    // Tick tick(__FN__);
    auto const& com_pub = input.com_pub;

    std::array<parallel::VoidTask, 2> tasks;
    bool ret1 = false;
    tasks[0] = [&ret1, &com_pub, &com_ext_pub, &challenge, &sub_proof,
                &input]() {
      auto const& xi = com_pub.xi;
      auto const& delta = com_ext_pub.delta;
      G1 left = xi * challenge + delta;
      G1 right = pc::ComputeCom(input.get_gx, sub_proof.z, sub_proof.z_delta);
      ret1 = left == right;
    };

    bool ret2 = false;
    tasks[1] = [&ret2, &com_pub, &com_ext_pub, &challenge, &sub_proof,
                &input]() {
      auto const& tau = com_pub.tau;
      auto const& beta = com_ext_pub.beta;
      G1 left = tau * challenge + beta;
      auto ip_za = InnerProduct(sub_proof.z, input.a);
      G1 right = pc::ComputeCom(input.gy, ip_za, sub_proof.z_beta);
      ret2 = left == right;
    };
    parallel::Invoke(tasks, true);

    assert(ret1 && ret2);
    return ret1 && ret2;
  }

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec const& com_sec,
                         ProveInput const& input) {
    // Tick tick(__FN__);
    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&com_pub, &input, &com_sec]() {
      com_pub.xi = pc::ComputeCom(input.get_gx, input.x, com_sec.r_xi);
    };
    tasks[1] = [&com_pub, &input, &com_sec]() {
      com_pub.tau = pc::ComputeCom(input.gy, input.y, com_sec.r_tau);
    };
    parallel::Invoke(tasks, true);
  }

  // com(n) + com(1) + ip(n)
  static void ComputeCommitmentExt(CommitmentExtPub& com_ext_pub,
                                   CommitmentExtSec& com_ext_sec,
                                   ProveInput const& input) {
    // Tick tick(__FN__);
    auto n = input.n();
    com_ext_sec.d.resize(n);
    FrRand(com_ext_sec.d.data(), n);
    com_ext_sec.r_beta = FrRand();
    com_ext_sec.r_delta = FrRand();

    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&com_ext_pub, &com_ext_sec, &input]() {
      com_ext_pub.delta =
          pc::ComputeCom(input.get_gx, com_ext_sec.d, com_ext_sec.r_delta);
    };
    tasks[1] = [&com_ext_pub, &input, &com_ext_sec]() {
      com_ext_pub.beta = pc::ComputeCom(
          input.gy, InnerProduct(input.a, com_ext_sec.d), com_ext_sec.r_beta);
    };
    parallel::Invoke(tasks, true);

    // std::cout << Tick::GetIndentString() << "multiexp(" << input.n() <<
    // ")\n";
  }

  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                         CommitmentExtPub const& com_ext_pub) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_pub.xi);
    HashUpdate(hash, com_pub.tau);
    HashUpdate(hash, com_ext_pub.beta);
    HashUpdate(hash, com_ext_pub.delta);
    hash.Final(seed.data());
  }

  static void ComputeSubProof(SubProof& sub_proof, ProveInput const& input,
                              CommitmentSec const& com_sec,
                              CommitmentExtSec const& com_ext_sec,
                              Fr const& challenge) {
    // z = c * x + d
    sub_proof.z = input.x * challenge + com_ext_sec.d;
    sub_proof.z_delta = challenge * com_sec.r_xi + com_ext_sec.r_delta;
    sub_proof.z_beta = challenge * com_sec.r_tau + com_ext_sec.r_beta;
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input,
                    CommitmentPub const& com_pub,
                    CommitmentSec const& com_sec) {
    Tick tick(__FN__, input.to_string());

    CommitmentExtSec com_ext_sec;
    ComputeCommitmentExt(proof.com_ext_pub, com_ext_sec, input);

    UpdateSeed(seed, com_pub, proof.com_ext_pub);
    Fr challenge = H256ToFr(seed);

    ComputeSubProof(proof.sub_proof, input, com_sec, com_ext_sec, challenge);
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FN__, input.to_string());
    if (input.a.size() != proof.sub_proof.z.size() || input.a.empty())
      return false;

    UpdateSeed(seed, input.com_pub, proof.com_ext_pub);
    Fr challenge = H256ToFr(seed);

    return VerifyInternal(input, challenge, proof.com_ext_pub, proof.sub_proof);
  }

  static bool Test(int64_t n);
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, A2::CommitmentPub const& t) {
  ar& YAS_OBJECT_NVP("a2.cp", ("xi", t.xi), ("tau", t.tau));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A2::CommitmentPub& t) {
  ar& YAS_OBJECT_NVP("a2.cp", ("xi", t.xi), ("tau", t.tau));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A2::CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A2::CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A2::SubProof const& t) {
  ar& YAS_OBJECT_NVP("a2.sp", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A2::SubProof& t) {
  ar& YAS_OBJECT_NVP("a2.sp", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A2::Proof const& t) {
  ar& YAS_OBJECT_NVP("a2.pf", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A2::Proof& t) {
  ar& YAS_OBJECT_NVP("a2.pf", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

bool A2::Test(int64_t n) {
  Tick tick(__FN__);
  std::cout << "n = " << n << "\n";
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  h256_t UpdateSeed = misc::RandH256();

  int64_t x_g_offset = 30;
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  auto z = InnerProduct(x, a);
  ProveInput prove_input("test", x, a, z, get_gx, pc::PcU());

  CommitmentPub com_pub;
  CommitmentSec com_sec(FrRand(), FrRand());
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

  VerifyInput verify_input("test", a, com_pub, get_gx, pc::PcU());
  bool success = Verify(proof, UpdateSeed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace hyrax
