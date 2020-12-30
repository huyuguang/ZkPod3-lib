#pragma once

#include "bp/bp.h"
#include "clink/equal_ip.h"
#include "groth09/details.h"
#include "hyrax/a3.h"
#include "utils/fst.h"

// NOTE: it does not exist in gro09 paper, just use the name.
// t: public vector<Fr>, size = n
// X, Y: secret vector<Fr>, size = n
// z: secret Fr
// open: com(gx, X), com(gy, Y), com(gz, z). gx, gy maybe overlapped, gz can not
// gyt MUST NOT overlapped with gy
// prove: z = <X, Y o t>
// proof size: 2log(n) Fr and 4 G1
// prove cost:
// verify cost:

namespace groth09 {

struct Sec51c {
  struct ProveInput {
    std::vector<Fr> const& x;   // size = n
    std::vector<Fr> const& y;   // size = n
    std::vector<Fr> const& t;   // size = n
    std::vector<Fr> const& yt;  // yt = y o t
    Fr const z;                 // z = <x, y o t>
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 get_gyt;
    G1 const& gz;

    int64_t n() const { return (int64_t)x.size(); }
    ProveInput(std::vector<Fr> const& x, std::vector<Fr> const& y,
               std::vector<Fr> const& t, std::vector<Fr> const& yt, Fr const& z,
               GetRefG1 const& get_gx, GetRefG1 const& get_gy, G1 const& gz)
        : x(x),
          y(y),
          t(t),
          yt(yt),
          z(z),
          get_gx(get_gx),
          get_gy(get_gy),
          gz(gz) {
      assert(!x.empty());
      assert(x.size() == y.size());
      assert(x.size() == t.size());
      assert(x.size() == yt.size());
      assert(yt == HadamardProduct(y, t));
      assert(z == InnerProduct(x, yt));
      get_gyt = [this](int64_t i) -> G1 const& {
        size_t offset = pc::Base::GSize() - this->n();
        return pc::PcG()[offset + i];
      };
    }
  };

  struct CommitmentPub {
    CommitmentPub() {}
    CommitmentPub(G1 const& a, G1 const& b, G1 const& c) : a(a), b(b), c(c) {}
    G1 a;
    G1 b;
    G1 c;
  };

  struct CommitmentSec {
    CommitmentSec() {}
    CommitmentSec(Fr const& r, Fr const& s, Fr const& t) : r(r), s(s), t(t) {}
    Fr r;
    Fr s;
    Fr t;
  };

  struct Proof {
    G1 com_yt;
    clink::EqualIp<hyrax::A3>::Proof proof_eip;
    bp::p31::Proof proof_bp;
    bool CheckFormat() const { return true; }
    bool operator==(Proof const& right) const {
      return com_yt == right.com_yt && proof_eip == right.proof_eip &&
             proof_bp == right.proof_bp;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }
  };

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input,
                    CommitmentPub const& com_pub,
                    CommitmentSec const& com_sec) {
    Tick tick(__FN__, std::to_string(input.n()));
    UpdateSeed(seed, com_pub, input.t);
    Fr com_yt_r = FrRand();
    proof.com_yt = pc::ComputeCom(input.get_gyt, input.yt, com_yt_r);
    UpdateSeed(seed, proof.com_yt);

    std::array<parallel::VoidTask, 2> tasks;

    tasks[0] = [&proof, &seed, &input, &com_pub, &com_sec, &com_yt_r]() {
      // prove yt = <y o t>
      ProveYt(proof, seed, input, com_pub, com_sec, com_yt_r);
    };

    tasks[1] = [&input, &com_pub, &com_sec, &proof, &com_yt_r, &seed]() {
      // now we have com(yt_g_offset, yt, com_yt_r), want to prove <x, yt> == z
      std::vector<G1> g1 = pc::CopyG(input.get_gx, input.n());
      std::vector<G1> g2 = pc::CopyG(input.get_gyt, input.n());
      std::vector<Fr> x_copy = input.x;
      std::vector<Fr> yt_copy = input.yt;
      bp::p31::ProveInput input_bp(std::move(g1), std::move(g2), pc::PcH(),
                                   input.gz, std::move(x_copy),
                                   std::move(yt_copy), input.z);
      bp::p31::CommitmentPub com_pub_bp;
      bp::p31::CommitmentSec com_sec_bp;
      com_sec_bp.alpha = com_sec.r + com_yt_r;
      com_sec_bp.beta = com_sec.t;
      com_pub_bp.p = com_pub.a + proof.com_yt;
      com_pub_bp.q = com_pub.c;
      bp::p31::Prove(proof.proof_bp, seed, std::move(input_bp), com_pub_bp,
                     com_sec_bp);
    };
    parallel::Invoke(tasks);
  }

  struct VerifyInput {
    VerifyInput(std::vector<Fr> const& t, CommitmentPub const& com_pub,
                GetRefG1 const& get_gx, GetRefG1 const& get_gy, G1 const& gz)
        : t(t), com_pub(com_pub), get_gx(get_gx), get_gy(get_gy), gz(gz) {
      get_gyt = [this](int64_t i) -> G1 const& {
        size_t offset = pc::Base::GSize() - this->n();
        return pc::PcG()[offset + i];
      };
    }
    int64_t n() const { return (int64_t)t.size(); }
    std::vector<Fr> const& t;
    CommitmentPub const& com_pub;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gy;
    GetRefG1 get_gyt;
    G1 const& gz;
  };

  static bool VerifyYt(Proof const& proof, h256_t seed,
                       VerifyInput const& input) {
    std::vector<Fr> e(input.n());
    ComputeFst1(seed, "sec51c", e);
    auto et = HadamardProduct(e, input.t);

    using eip = clink::EqualIp<hyrax::A3>;
    eip::VerifyInput input_eip(e, proof.com_yt, input.get_gyt, et,
                               input.com_pub.b, input.get_gy);
    return eip::Verify(seed, proof.proof_eip, input_eip);
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    UpdateSeed(seed, input.com_pub, input.t);
    UpdateSeed(seed, proof.com_yt);

    std::array<std::atomic<bool>, 2> rets;
    std::array<parallel::VoidTask, 2> tasks;

    tasks[0] = [&rets, &proof, &seed, &input]() {
      rets[0] = VerifyYt(proof, seed, input);
      assert(rets[0]);
    };

    tasks[1] = [&rets, &input, &proof, &seed]() {
      std::vector<G1> g1 = pc::CopyG(input.get_gx, input.n());
      std::vector<G1> g2 = pc::CopyG(input.get_gyt, input.n());
      bp::p31::CommitmentPub com_pub_bp;
      com_pub_bp.p = input.com_pub.a + proof.com_yt;
      com_pub_bp.q = input.com_pub.c;
      bp::p31::VerifyInput input_p3(std::move(g1), std::move(g2), pc::PcH(),
                                    input.gz, com_pub_bp);
      rets[1] = bp::p31::Verify(proof.proof_bp, seed, std::move(input_p3));
    };

    parallel::Invoke(tasks);

    if (!rets[0] || !rets[1]) {
      assert(false);
      return false;
    }
    return true;
  }

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                         ProveInput const& input) {
    // Tick tick(__FN__);
    com_sec.r = FrRand();
    com_sec.s = FrRand();
    com_sec.t = FrRand();

    std::array<parallel::VoidTask, 3> tasks;
    tasks[0] = [&com_pub, &input, &com_sec]() {
      com_pub.a = pc::ComputeCom(input.get_gx, input.x, com_sec.r);
    };
    tasks[1] = [&com_pub, &input, &com_sec]() {
      com_pub.b = pc::ComputeCom(input.get_gy, input.y, com_sec.s);
    };
    tasks[2] = [&com_pub, &input, &com_sec]() {
      com_pub.c = pc::ComputeCom(input.gz, input.z, com_sec.t);
    };
    parallel::Invoke(tasks);
  }

  static bool Test(int64_t n);

 private:
  static void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                         std::vector<Fr> const& t) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, com_pub.a);
    HashUpdate(hash, com_pub.b);
    HashUpdate(hash, com_pub.c);
    HashUpdate(hash, t);
    hash.Final(seed.data());
  }

  static void UpdateSeed(h256_t& seed, G1 const& g) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, g);
    hash.Final(seed.data());
  }

  static void ProveYt(Proof& proof, h256_t seed, ProveInput const& input,
                      CommitmentPub const& com_pub,
                      CommitmentSec const& com_sec, Fr const& com_yt_r) {
    Tick tick(__FN__, std::to_string(input.n()));
    std::vector<Fr> e(input.n());
    ComputeFst1(seed, "sec51c", e);

    Fr yte = InnerProduct(input.yt, e);
    auto et = HadamardProduct(e, input.t);
    DCHECK(yte == InnerProduct(input.y, et), "");

    // prove <yt, e> == <y, et>
    using eip = clink::EqualIp<hyrax::A3>;
    eip::ProveInput input_eip(input.yt, e, proof.com_yt, com_yt_r,
                              input.get_gyt, input.y, et, com_pub.b, com_sec.s,
                              input.get_gy, yte);
    eip::Prove(proof.proof_eip, seed, input_eip);
  }
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, Sec51c::Proof const& t) {
  ar& YAS_OBJECT_NVP("51c.pf", ("c", t.com_yt), ("eip", t.proof_eip),
                     ("bp", t.proof_bp));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Sec51c::Proof& t) {
  ar& YAS_OBJECT_NVP("51c.pf", ("c", t.com_yt), ("eip", t.proof_eip),
                     ("bp", t.proof_bp));
}

bool Sec51c::Test(int64_t n) {
  Tick tick(__FN__);
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> y(n);
  FrRand(y.data(), n);
  std::vector<Fr> t(n);
  FrRand(t.data(), t.size());

  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 10;
  int64_t y_g_offset = 30;
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  GetRefG1 get_gy = [y_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[y_g_offset + i];
  };

  std::vector<Fr> yt(n);
  yt = HadamardProduct(y, t);
  Fr z = InnerProduct(x, yt);
  ProveInput prove_input(x, y, t, yt, z, get_gx, get_gy, pc::PcU());

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prove_input);

  Proof proof;
  Prove(proof, seed, prove_input, com_pub, com_sec);

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

  VerifyInput verify_input(t, com_pub, get_gx, get_gy, pc::PcU());
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace groth09
