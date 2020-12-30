#pragma once

#include "./protocol1.h"

// protocol 3.1
// a, b: secret vector<Fr>
// c = <a,b>
// open com(a,b): P=h*alpha + g1*a + g2*b
// open com(c): Q=h*beta + u * c
// prove: c = <a,b>

namespace bp::p31 {

struct ProveInput {
  ProveInput(std::vector<G1>&& ig1, std::vector<G1>&& ig2, G1 const& h,
             G1 const& u, std::vector<Fr>&& ia, std::vector<Fr>&& ib,
             Fr const& c)
      : n((int64_t)ig1.size()),
        g1(std::move(ig1)),
        g2(std::move(ig2)),
        h(h),
        u(u),
        a(std::move(ia)),
        b(std::move(ib)),
        c(c) {
    assert(g1.size() == g2.size());
    assert(g1.size() == a.size());
    assert(g1.size() == b.size());
    assert(c == InnerProduct(a, b));
  }
  int64_t const n;
  std::vector<G1> g1;
  std::vector<G1> g2;
  G1 h;
  G1 u;
  std::vector<Fr> a;
  std::vector<Fr> b;
  Fr c;
};

struct CommitmentPub {
  G1 p;
  G1 q;
};

struct CommitmentSec {
  Fr alpha;
  Fr beta;
};

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProveInput const& input) {
  com_sec.alpha = FrRand();
  com_sec.beta = FrRand();

  // q = h*beta + u * c
  com_pub.q = input.h * com_sec.beta + input.u * input.c;

  // p = h*alpha + g1*a + g2*b
  auto get_g = [&input](int64_t i) -> G1 const& {
    if (i == 0) return input.h;
    if (i > 0 && i <= input.n) return input.g1[i - 1];
    return input.g2[i - 1 - input.n];
  };
  auto get_f = [&input, &com_sec](int64_t i) -> Fr const& {
    if (i == 0) return com_sec.alpha;
    if (i > 0 && i <= input.n) return input.a[i - 1];
    return input.b[i - 1 - input.n];
  };
  com_pub.p = MultiExpBdlo12<G1>(get_g, get_f, input.n * 2 + 1);
}

struct CommitmentExtPub {
  G1 r;
  G1 t1;
  G1 t2;
};

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.r == right.r && left.t1 == right.t1 && left.t2 == right.t2;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("bp3.cep", ("r", t.r), ("t1", t.t1), ("t2", t.t2));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("bp3.cep", ("r", t.r), ("t1", t.t1), ("t2", t.t2));
}

struct CommitmentExtSec {
  CommitmentExtSec(int64_t n) : da(n), db(n) {
    FrRand(da);
    FrRand(db);
    rou = FrRand();
    tau1 = FrRand();
    tau2 = FrRand();
  }
  std::vector<Fr> da;
  std::vector<Fr> db;
  Fr rou;
  Fr tau1;
  Fr tau2;
};

inline void ComputeComExt(CommitmentExtPub& com_ext_pub,
                          CommitmentExtSec const& com_ext_sec,
                          ProveInput const& input) {
  auto get_g = [&input](int64_t i) -> G1 const& {
    if (i == 0) return input.h;
    if (i > 0 && i <= input.n) return input.g1[i - 1];
    return input.g2[i - 1 - input.n];
  };
  auto get_f = [&input, &com_ext_sec](int64_t i) -> Fr const& {
    if (i == 0) return com_ext_sec.rou;
    if (i > 0 && i <= input.n) return com_ext_sec.da[i - 1];
    return com_ext_sec.db[i - 1 - input.n];
  };
  com_ext_pub.r = MultiExpBdlo12<G1>(get_g, get_f, input.n * 2 + 1);
  Fr t1 = InnerProduct(input.a, com_ext_sec.db) +
          InnerProduct(com_ext_sec.da, input.b);
  Fr t2 = InnerProduct(com_ext_sec.da, com_ext_sec.db);
  com_ext_pub.t1 = input.h * com_ext_sec.tau1 + input.u * t1;
  com_ext_pub.t2 = input.h * com_ext_sec.tau2 + input.u * t2;
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub, int64_t n) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.p);
  HashUpdate(hash, com_pub.q);
  HashUpdate(hash, com_ext_pub.r);
  HashUpdate(hash, com_ext_pub.t1);
  HashUpdate(hash, com_ext_pub.t2);
  HashUpdate(hash, n);
  hash.Final(seed.data());
}

struct Proof {
  Fr c;
  Fr mu;
  Fr tau;
  CommitmentExtPub com_ext_pub;
  p1::Proof p1;
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.c == right.c && left.mu == right.mu && left.tau == right.tau &&
         left.com_ext_pub == right.com_ext_pub && left.p1 == right.p1;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("bp3.pf", ("c", t.c), ("m", t.mu), ("t", t.tau),
                     ("cep", t.com_ext_pub), ("p1", t.p1));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("bp3.pf", ("c", t.c), ("m", t.mu), ("t", t.tau),
                     ("cep", t.com_ext_pub), ("p1", t.p1));
}

inline void Prove(Proof& proof, h256_t seed, ProveInput&& input,
                  CommitmentPub const& com_pub, CommitmentSec const& com_sec) {
  Tick tick(__FN__, std::to_string(input.n));
  CommitmentExtSec com_ext_sec(input.n);
  ComputeComExt(proof.com_ext_pub, com_ext_sec, input);
  CommitmentExtPub& com_ext_pub = proof.com_ext_pub;

  UpdateSeed(seed, com_pub, com_ext_pub, input.n);
  Fr x = H256ToFr(seed);

  auto a2 = input.a + com_ext_sec.da * x;
  auto b2 = input.b + com_ext_sec.db * x;

  proof.c = InnerProduct(a2, b2);
  proof.mu = com_sec.alpha + x * com_ext_sec.rou;
  proof.tau = com_sec.beta + x * com_ext_sec.tau1 + x * x * com_ext_sec.tau2;

  G1 p = com_pub.p + com_ext_pub.r * x + input.h * (-proof.mu);
  p1::Prove(proof.p1, seed, std::move(input.g1), std::move(input.g2),
            std::move(a2), std::move(b2), p, proof.c);
}

struct VerifyInput {
  VerifyInput(std::vector<G1>&& ig1, std::vector<G1>&& ig2, G1 const& h,
              G1 const& u, CommitmentPub const& com_pub)
      : n((int64_t)ig1.size()),
        g1(std::move(ig1)),
        g2(std::move(ig2)),
        h(h),
        u(u),
        com_pub(com_pub) {
    assert(g1.size() == g2.size());
  }

  int64_t const n;
  std::vector<G1> g1;
  std::vector<G1> g2;
  G1 h;
  G1 u;
  CommitmentPub const& com_pub;
};

inline bool Verify(Proof const& proof, h256_t seed, VerifyInput&& input) {
  UpdateSeed(seed, input.com_pub, proof.com_ext_pub, input.n);
  Fr x = H256ToFr(seed);
  G1 left = input.com_pub.q + proof.com_ext_pub.t1 * x +
            proof.com_ext_pub.t2 * (x * x);
  G1 right = input.h * proof.tau + input.u * proof.c;
  if (left != right) {
    assert(false);
    return false;
  }
  G1 p = input.com_pub.p + proof.com_ext_pub.r * x + input.h * (-proof.mu);
  if (p != proof.p1.p) {
    assert(false);
    return false;
  }
  return p1::Verify(seed, std::move(input.g1), std::move(input.g2), proof.p1);
}

inline bool Test(int64_t n) {
  Tick tick(__FN__);
  h256_t seed = misc::RandH256();
  std::vector<G1> g1(n);
  G1Rand(g1.data(), n);

  std::vector<G1> g2(n);
  G1Rand(g2.data(), n);

  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  std::vector<Fr> b(n);
  FrRand(b.data(), n);

  Fr c = InnerProduct(a, b);

  G1 h = G1Rand();
  G1 u = G1Rand();

  auto g1_copy = g1;
  auto g2_copy = g2;
  ProveInput prove_input(std::move(g1), std::move(g2), h, u, std::move(a),
                         std::move(b), c);

  CommitmentSec com_sec;
  CommitmentPub com_pub;
  ComputeCom(com_pub, com_sec, prove_input);

  Proof proof;
  Prove(proof, seed, std::move(prove_input), com_pub, com_sec);

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

  VerifyInput verify_input(std::move(g1_copy), std::move(g2_copy), h, u,
                           com_pub);
  bool success = Verify(proof, seed, std::move(verify_input));
  std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace bp::p31