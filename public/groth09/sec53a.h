#pragma once

#include <map>
#include <memory>

#include "groth09/details.h"
#include "groth09/sec51a.h"
#include "utils/fst.h"

// x, y: secret matric<Fr>, size = m*n
// z: secret Fr
// open: com(gx, x1), com(gy, y1), com(gx, x2), com(gy, y2) ... com(gz, z)
// prove: z = <x1,y1> + <x2,y2>...
// proof size: 2*log(m)+4 G1, 2n+3 Fr
// prove cost: 2*log(m)*mulexp(n)
// verify cost: mulexp(n)
// base on sec51a
namespace groth09::sec53a {

struct CommitmentPub {
  std::vector<G1> a;  // a.size = m
  std::vector<G1> b;  // b.size = m
  G1 c;
  int64_t m() const { return a.size(); }
  bool CheckFormat() const {
    if (a.empty() || a.size() != b.size()) return false;
    return misc::Pow2UB(m()) == (uint64_t)m();
  }
  void Align() {
    int64_t old_m = a.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      auto const& g0 = G1Zero();
      a.resize(new_m);
      std::fill(a.begin() + old_m, a.end(), g0);
      b.resize(new_m);
      std::fill(b.begin() + old_m, b.end(), g0);
    }
  }
};

struct CommitmentSec {
  std::vector<Fr> r;  // r.size = m
  std::vector<Fr> s;  // s.size = m
  Fr t;
  void Align() {
    int64_t old_m = r.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      static const Fr f0 = FrZero();
      r.resize(new_m);
      std::fill(r.begin() + old_m, r.end(), f0);
      s.resize(new_m);
      std::fill(s.begin() + old_m, s.end(), f0);
    }
  }
};

struct VerifierInput {
  VerifierInput(CommitmentPub const& com_pub, int64_t x_g_offset,
                int64_t y_g_offset, int64_t z_g_offset)
      : com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {}
  CommitmentPub const& com_pub;
  int64_t g_offset = 0;
  int64_t m() const { return com_pub.m(); }
  bool CheckFormat() const { return com_pub.CheckFormat(); }
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
};

struct ProverInput {
  std::vector<std::vector<Fr>> x;
  std::vector<std::vector<Fr>> y;
  Fr z;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;

  int64_t m() const { return x.size(); }
  int64_t n() const { return x[0].size(); }

  ProverInput(std::vector<std::vector<Fr>> ix, std::vector<std::vector<Fr>> iy,
              Fr const& z, int64_t x_g_offset, int64_t y_g_offset,
              int64_t z_g_offset)
      : x(std::move(ix)),
        y(std::move(iy)),
        z(z),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {
#ifdef _DEBUG
    assert(!x.empty() && x.size() == y.size());
    Fr check_z = FrZero();
    for (int64_t i = 0; i < m(); ++i) {
      assert(x[i].size() == (size_t)n());
      assert(y[i].size() == (size_t)n());
      check_z += InnerProduct(x[i], y[i]);
    }
    assert(z == check_z);
#endif
  }

  // pad some trivial values
  void Align() {
    // Tick tick(__FUNCTION__);
    int64_t old_m = m();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (old_m == new_m) return;

    static const Fr f0 = FrZero();

    x.resize(new_m);
    y.resize(new_m);

    for (int64_t i = old_m; i < new_m; ++i) {
      auto& x_i = x[i];
      x_i.resize(n());
      std::fill(x_i.begin(), x_i.end(), f0);
      auto& y_i = y[i];
      y_i.resize(n());
      std::fill(y_i.begin(), y_i.end(), f0);
    }
  }

  void Update(Fr const& sigma_xy1, Fr const& sigma_xy2, Fr const& e,
              Fr const& ee) {
    // Tick tick(__FUNCTION__);
    auto m2 = m() / 2;
    for (int64_t i = 0; i < m2; ++i) {
      x[i] = x[2 * i + 1] * e + x[2 * i];
    }
    x.resize(m2);

    for (int64_t i = 0; i < m2; ++i) {
      y[i] = y[2 * i] * e + y[2 * i + 1];
    }
    y.resize(m2);

    z = sigma_xy1 * ee + z * e + sigma_xy2;
  }
};

struct CommitmentExtPub {
  // 5.3 Recursive
  std::vector<G1> cl;  // size = log(m)
  std::vector<G1> cu;  // size = log(m)
  bool CheckFormat(int64_t check_m) const {
    if (cl.size() != cu.size()) return false;
    return m() == check_m;
  }
  int64_t m() const { return 1LL << cl.size(); }
};

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.cl == right.cl && left.cu == right.cu;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("53.cep", ("cl", t.cl), ("cu", t.cu));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("53.cep", ("cl", t.cl), ("cu", t.cu));
}

struct RomProof {
  CommitmentExtPub com_ext_pub;   // 2*log(m) G1
  sec51a::RomProof rom_proof_51;  // 4 G1, 2n+3 Fr

  int64_t n() const { return rom_proof_51.n(); }
  int64_t m() const { return com_ext_pub.m(); }

  bool CheckFormat(int64_t check_m) const {
    return com_ext_pub.CheckFormat(check_m) && rom_proof_51.CheckFormat();
  }
};

inline bool operator==(RomProof const& left, RomProof const& right) {
  return left.com_ext_pub == right.com_ext_pub &&
         left.rom_proof_51 == right.rom_proof_51;
}

inline bool operator!=(RomProof const& left, RomProof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, RomProof const& t) {
  ar& YAS_OBJECT_NVP("53.rp", ("c", t.com_ext_pub), ("r", t.rom_proof_51));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, RomProof& t) {
  ar& YAS_OBJECT_NVP("53.rp", ("c", t.com_ext_pub), ("r", t.rom_proof_51));
}

inline void ComputeCom(ProverInput const& input, CommitmentPub* com_pub,
                       CommitmentSec const& com_sec) {
  // Tick tick(__FUNCTION__);
  auto const m = input.m();
  auto const n = input.n();

  com_pub->a.resize(m);
  com_pub->b.resize(m);

  auto parallel_f = [&input, &com_pub, &com_sec](int64_t i) mutable {
    com_pub->a[i] =
        PcComputeCommitmentG(input.x_g_offset, input.x[i], com_sec.r[i]);
    com_pub->b[i] =
        PcComputeCommitmentG(input.y_g_offset, input.y[i], com_sec.s[i]);
  };
  parallel::For(m, parallel_f, n < 16 * 1024);

  com_pub->c = PcComputeCommitmentG(input.z_g_offset,input.z, com_sec.t);
}

inline void ComputeCom(ProverInput const& input, CommitmentPub* com_pub,
                       CommitmentSec* com_sec) {
  // Tick tick(__FUNCTION__);
  auto const m = input.m();
  com_sec->r.resize(m);
  FrRand(com_sec->r.data(), m);

  com_sec->s.resize(m);
  FrRand(com_sec->s.data(), m);

  com_sec->t = FrRand();

  ComputeCom(input, com_pub, *com_sec);
}

inline void RomProveFinal(RomProof& rom_proof, h256_t const& seed,
                          ProverInput const& input,
                          CommitmentPub const& com_pub,
                          CommitmentSec const& com_sec) {
  // Tick tick(__FUNCTION__);
  assert(input.m() == 1);

  sec51a::ProverInput input_51(input.x[0], input.y[0], input.z,input.x_g_offset,
                               input.y_g_offset,input.z_g_offset);
  sec51a::CommitmentPub com_pub_51(com_pub.a[0], com_pub.b[0], com_pub.c);
  sec51a::CommitmentSec com_sec_51(com_sec.r[0], com_sec.s[0], com_sec.t);
  sec51a::RomProve(rom_proof.rom_proof_51, seed, input_51, com_pub_51,
                   com_sec_51);
}

inline void ComputeSigmaXY(ProverInput const& input, Fr* sigma_xy1,
                           Fr* sigma_xy2) {
  // Tick tick(__FUNCTION__);
  int64_t m = input.m();
  auto m2 = m / 2;
  std::vector<Fr> xy1(m2, FrZero());
  std::vector<Fr> xy2(m2, FrZero());
  auto parallel_f = [&input, &xy1, &xy2](int64_t i) {
    auto const& x1 = input.x[2 * i + 1];
    auto const& y1 = input.y[2 * i];
    xy1[i] = InnerProduct(x1, y1);

    auto const& x2 = input.x[2 * i];
    auto const& y2 = input.y[2 * i + 1];
    xy2[i] = InnerProduct(x2, y2);
  };
  parallel::For(m2, parallel_f, m2 < 1024);

  *sigma_xy1 = parallel::Accumulate(xy1.begin(), xy1.end(), FrZero());
  *sigma_xy2 = parallel::Accumulate(xy2.begin(), xy2.end(), FrZero());
}

inline void UpdateCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                      Fr const& tl, Fr const& tu, G1 const& cl, G1 const& cu,
                      Fr const& e, Fr const& ee) {
  // Tick tick(__FUNCTION__);
  CommitmentPub com_pub2;
  CommitmentSec com_sec2;
  auto m2 = com_pub.a.size() / 2;
  com_pub2.a.resize(m2);
  com_pub2.b.resize(m2);
  com_sec2.r.resize(m2);
  com_sec2.s.resize(m2);

  auto parallel_f = [&com_pub, &com_sec, &com_pub2, &com_sec2, &e](int64_t i) {
    auto& a2 = com_pub2.a;
    auto const& a = com_pub.a;
    a2[i] = a[2 * i] + a[2 * i + 1] * e;

    auto& b2 = com_pub2.b;
    auto const& b = com_pub.b;
    b2[i] = b[2 * i] * e + b[2 * i + 1];

    auto& r2 = com_sec2.r;
    auto const& r = com_sec.r;
    r2[i] = r[2 * i] + r[2 * i + 1] * e;

    auto& s2 = com_sec2.s;
    auto const& s = com_sec.s;
    s2[i] = s[2 * i] * e + s[2 * i + 1];
  };
  parallel::For((int64_t)m2, parallel_f, m2 < 1024);

  com_pub2.c = cl * ee + com_pub.c * e + cu;
  com_sec2.t = tl * ee + com_sec.t * e + tu;
  com_pub = std::move(com_pub2);
  com_sec = std::move(com_sec2);
}

inline Fr ComputeChallenge(h256_t const& seed, CommitmentPub const& com_pub,
                           G1 const& cl, G1 const& cu) {
  // Tick tick(__FUNCTION__);
  CryptoPP::Keccak_256 hash;
  h256_t digest;
  HashUpdate(hash, seed);
  HashUpdate(hash, cl);
  HashUpdate(hash, cu);
  HashUpdate(hash, com_pub.a);
  HashUpdate(hash, com_pub.b);
  HashUpdate(hash, com_pub.c);
  hash.Final(digest.data());
  return H256ToFr(digest);
}

// pad some trivial value
inline void AlignData(ProverInput& input, CommitmentPub& com_pub,
                      CommitmentSec& com_sec) {
  // Tick tick(__FUNCTION__);
  input.Align();
  com_sec.Align();
  com_pub.Align();
}

inline void RomProveRecursive(RomProof& rom_proof, h256_t& seed,
                              ProverInput& input, CommitmentPub& com_pub,
                              CommitmentSec& com_sec) {
  // Tick tick(__FUNCTION__);
  assert(input.m() > 1);

  Fr sigma_xy1, sigma_xy2;
  ComputeSigmaXY(input, &sigma_xy1, &sigma_xy2);

  // compute cl, cu
  Fr tl = FrRand();
  Fr tu = FrRand();
  G1 cl = PcComputeCommitmentG(input.z_g_offset,sigma_xy1, tl);
  G1 cu = PcComputeCommitmentG(input.z_g_offset,sigma_xy2, tu);
  rom_proof.com_ext_pub.cl.push_back(cl);
  rom_proof.com_ext_pub.cu.push_back(cu);

  // rom challenge
  Fr e = ComputeChallenge(seed, com_pub, cl, cu);
  Fr ee = e * e;
  seed = FrToBin(e);

  input.Update(sigma_xy1, sigma_xy2, e, ee);

  UpdateCom(com_pub, com_sec, tl, tu, cl, cu, e, ee);
  // debug check com_pub2 and com_sec2
#ifdef _DEBUG
  CommitmentPub check_com_pub;
  ComputeCom(input, &check_com_pub, com_sec);
  assert(check_com_pub.a == com_pub.a);
  assert(check_com_pub.b == com_pub.b);
  assert(check_com_pub.c == com_pub.c);
#endif
}

inline void RomProve(RomProof& rom_proof, h256_t seed, ProverInput input,
                     CommitmentPub com_pub, CommitmentSec com_sec) {
  // Tick tick(__FUNCTION__);
  assert(PcBase::kGSize >= input.n());

  while (input.m() > 1) {
    RomProveRecursive(rom_proof, seed, input, com_pub, com_sec);
  }
  return RomProveFinal(rom_proof, seed, input, com_pub, com_sec);
}

inline bool RomVerify(RomProof const& rom_proof, h256_t seed,
                      VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  if (!rom_proof.CheckFormat(input.m())) {
    assert(false);
    return false;
  }

  CommitmentPub com_pub = input.com_pub;

  for (size_t loop = 0; loop < rom_proof.com_ext_pub.cl.size(); ++loop) {
    // rom challenge
    auto const& cl = rom_proof.com_ext_pub.cl[loop];
    auto const& cu = rom_proof.com_ext_pub.cu[loop];
    Fr e = ComputeChallenge(seed, com_pub, cl, cu);
    Fr ee = e * e;
    seed = FrToBin(e);

    std::vector<G1> a2(com_pub.m() / 2);
    std::vector<G1> b2(com_pub.m() / 2);
    G1 c2;

    auto m2 = com_pub.m() / 2;
    auto parallel_f = [&com_pub, &a2, &b2, &e](int64_t i) {
      auto const& a = com_pub.a;
      a2[i] = a[2 * i] + a[2 * i + 1] * e;

      auto const& b = com_pub.b;
      b2[i] = b[2 * i] * e + b[2 * i + 1];
    };
    parallel::For(m2, parallel_f, m2 < 1024);

    c2 = cl * ee + com_pub.c * e + cu;

    com_pub.a = std::move(a2);
    com_pub.b = std::move(b2);
    com_pub.c = std::move(c2);
  }

  assert(com_pub.m() == 1);

  sec51a::CommitmentPub com_pub_51(com_pub.a[0], com_pub.b[0], com_pub.c);
  sec51a::VerifierInput verifier_input_51(com_pub_51, input.x_g_offset,
                                          input.y_g_offset, input.z_g_offset);
  return sec51a::RomVerify(rom_proof.rom_proof_51, seed, verifier_input_51);
}

inline bool TestRom(int64_t m, int64_t n) {
  Tick tick(__FUNCTION__);
  std::cout << "m=" << m << ", n=" << n << "\n";

  std::vector<std::vector<Fr>> x(m);
  for (auto& i : x) {
    i.resize(n);
    FrRand(i.data(), n);
  }

  std::vector<std::vector<Fr>> y(m);
  for (auto& i : y) {
    i.resize(n);
    FrRand(i.data(), n);
  }

  Fr z = FrZero();
  for (int64_t i = 0; i < m; ++i) {
    z += InnerProduct(x[i], y[i]);
  }

  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 20;
  int64_t y_g_offset = 240;
  int64_t z_g_offset = -1;
  ProverInput prover_input(std::move(x), std::move(y), z, x_g_offset,
                           y_g_offset, z_g_offset);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(prover_input, &com_pub, &com_sec);

  AlignData(prover_input, com_pub, com_sec);

  RomProof rom_proof;
  RomProve(rom_proof, seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(com_pub, x_g_offset, y_g_offset, z_g_offset);
  bool success = RomVerify(rom_proof, seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}

}  // namespace groth09::sec53a