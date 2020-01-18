#pragma once

#include "./details.h"

// recursive version of a2
// a: public vector<Fr>, size = n
// x: secret vector<Fr>, size = n
// y: secret Fr, = <x,a>
// open: com(gx,x), com(gy,y)
// prove: y=<x,a>
// proof size: 2 G1 and 2log(n) Fr

namespace hyrax::a3 {

struct ProverInput {
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& a, Fr const& y,
              int64_t x_g_offset, int64_t y_g_offset)
      : x(x), a(a), y(y), x_g_offset(x_g_offset), y_g_offset(y_g_offset) {
    assert(x.size() == a.size() && !a.empty());
    assert(y == InnerProduct(x, a));
  }
  int64_t n() const { return (int64_t)x.size(); }
  std::vector<Fr> const& x;  // x.size = n
  std::vector<Fr> const& a;  // a.size = n
  Fr const y;                // y = <x, a>
  int64_t const x_g_offset;
  int64_t const y_g_offset;
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& xi, G1 const& tau) : xi(xi), tau(tau) {}
  G1 xi;   // com(x,r_xi)
  G1 tau;  // com(y,r_tau)
};
inline bool operator==(CommitmentPub const& left, CommitmentPub const& right) {
  return left.xi == right.xi && left.tau == right.tau;
}

inline bool operator!=(CommitmentPub const& left, CommitmentPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentPub const& t) {
  ar& YAS_OBJECT_NVP("a3.cp", ("xi", t.xi), ("tau", t.tau));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentPub& t) {
  ar& YAS_OBJECT_NVP("a3.cp", ("xi", t.xi), ("tau", t.tau));
}

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
  // recursive rounds
  std::vector<G1> gamma_neg_1;  // size=log(n)
  std::vector<G1> gamma_pos_1;
  // final round
  G1 delta;  // com(d, r_delta)
  G1 beta;   // com(<a,d>, r_beta)
};

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.delta == right.delta && left.beta == right.beta &&
         left.gamma_neg_1 == right.gamma_neg_1 &&
         left.gamma_pos_1 == right.gamma_pos_1;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("a3.cep", ("gn1", t.gamma_neg_1), ("gp1", t.gamma_pos_1),
                     ("delta", t.delta), ("beta", t.beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("a3.cep", ("gn1", t.gamma_neg_1), ("gp1", t.gamma_pos_1),
                     ("delta", t.delta), ("beta", t.beta));
}

struct CommitmentExtSec {
  // recursive rounds
  std::vector<Fr> r_gamma_neg_1;  // size=log(n)
  std::vector<Fr> r_gamma_pos_1;
  // finnal round
  Fr d;
  Fr r_beta;
  Fr r_delta;
};

struct Proof {
  Fr z1;
  Fr z2;
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.z1 == right.z1 && left.z2 == right.z2;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("a3.pf", ("z1", t.z1), ("z2", t.z2));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("a3.pf", ("z1", t.z1), ("z2", t.z2));
}

struct RomProof {
  CommitmentExtPub com_ext_pub;
  Proof proof;
  int64_t aligned_n() const { return 1LL << com_ext_pub.gamma_neg_1.size(); }
};

inline bool operator==(RomProof const& left, RomProof const& right) {
  return left.com_ext_pub == right.com_ext_pub && left.proof == right.proof;
}

inline bool operator!=(RomProof const& left, RomProof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, RomProof const& t) {
  ar& YAS_OBJECT_NVP("a3.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, RomProof& t) {
  ar& YAS_OBJECT_NVP("a3.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

struct VerifierInput {
  VerifierInput(std::vector<Fr> const& a, CommitmentPub const& com_pub,
                int64_t x_g_offset, int64_t y_g_offset)
      : a(a),
        com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset) {}
  std::vector<Fr> const& a;  // a.size = n
  CommitmentPub const& com_pub;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
};

inline void UpdateSeed(h256_t& seed, std::vector<Fr> const& a,
                       CommitmentPub const& com_pub) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, a);
  HashUpdate(hash, com_pub.xi);
  HashUpdate(hash, com_pub.tau);
  hash.Final(seed.data());
}

inline void UpdateSeed(h256_t& seed, G1 const& a, G1 const& b) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, a);
  HashUpdate(hash, b);
  hash.Final(seed.data());
}

// TODO: maybe can optimize
inline std::vector<G1> FuncO(std::vector<G1> const& g1, Fr const& k,
                             std::vector<G1> const& g2, Fr const& l) {
  assert(g1.size() == g2.size());
  std::vector<G1> g(g1.size());
  for (size_t i = 0; i < g1.size(); ++i) {
    g[i] = g1[i] * k + g2[i] * l;
  }
  return g;
}

template <typename T>
void Divide(std::vector<T> const& t, std::vector<T>& t1, std::vector<T>& t2,
            T const& t0) {
  auto n = t.size();
  auto half = misc::Pow2UB(n) / 2;
  t1.resize(half);
  t2.resize(half);
  std::copy(t.begin(), t.begin() + half, t1.begin());
  std::copy(t.begin() + half, t.end(), t2.begin());
  std::fill(t2.begin() + (n - half), t2.end(), t0);
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec const& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  std::array<parallel::Task, 2> tasks;
  tasks[0] = [&com_pub, &input, &com_sec]() {
    com_pub.xi = PcComputeCommitmentG(input.x_g_offset, input.x, com_sec.r_xi);
  };
  tasks[1] = [&com_pub, &input, &com_sec]() {
    com_pub.tau = PcComputeCommitmentG(input.y_g_offset,input.y, com_sec.r_tau);
  };
  parallel::Invoke(tasks, true);
}

inline void RomProve(RomProof& rom_proof, h256_t seed, ProverInput input,
                     CommitmentPub com_pub, CommitmentSec com_sec) {
  UpdateSeed(seed, input.a, com_pub);

  auto x = input.x;
  auto a = input.a;
  auto y = input.y;
  int64_t n = x.size();
  int64_t round = (int64_t)misc::Log2UB(n);
  auto gx = GetPcBase().CopyG(input.x_g_offset, n);
  gx.resize(misc::Pow2UB(n));
  std::fill(gx.begin() + n, gx.end(), G1Zero());
  auto const& h = PcH();
  rom_proof.com_ext_pub.gamma_neg_1.resize(round);
  rom_proof.com_ext_pub.gamma_pos_1.resize(round);
  CommitmentExtSec com_ext_sec;
  com_ext_sec.r_gamma_neg_1.resize(round);
  com_ext_sec.r_gamma_pos_1.resize(round);
  auto r_gamma = com_sec.r_tau + com_sec.r_xi;
  CommitmentExtPub& com_ext_pub = rom_proof.com_ext_pub;
  auto const& gy = PcG(input.y_g_offset);
  G1 gamma = com_pub.xi + com_pub.tau;

  // recursive round
  for (int64_t loop = 0; loop < round; ++loop) {
    std::vector<Fr> x1, x2;
    Divide(x, x1, x2, FrZero());
    std::vector<Fr> a1, a2;
    Divide(a, a1, a2, FrZero());
    std::vector<G1> g1, g2;
    Divide(gx, g1, g2, G1Zero());

    auto& gamma_neg_1 = com_ext_pub.gamma_neg_1[loop];
    auto& gamma_pos_1 = com_ext_pub.gamma_pos_1[loop];
    auto& r_gamma_neg_1 = com_ext_sec.r_gamma_neg_1[loop];
    r_gamma_neg_1 = FrRand();
    auto& r_gamma_pos_1 = com_ext_sec.r_gamma_pos_1[loop];
    r_gamma_pos_1 = FrRand();
    auto x1_a2 = InnerProduct(x1, a2);
    gamma_neg_1 = h * r_gamma_neg_1 + gy * x1_a2;
    gamma_neg_1 += MultiExpBdlo12(g2, x1);
    auto x2_a1 = InnerProduct(x2, a1);
    gamma_pos_1 = h * r_gamma_pos_1 + gy * x2_a1;
    gamma_pos_1 += MultiExpBdlo12(g1, x2);

    UpdateSeed(seed, gamma_neg_1, gamma_pos_1);
    Fr c = H256ToFr(seed);
    Fr cc = c * c;
    Fr c_inv = FrInv(c);
    Fr cc_inv = FrInv(cc);

    gamma += gamma_neg_1 * cc + gamma_pos_1 * cc_inv;
    a = a1 * c_inv + a2 * c;
    gx = FuncO(g1, c_inv, g2, c);
    x = x1 * c + x2 * c_inv;
    y += cc * x1_a2 + cc_inv * x2_a1;
    r_gamma += r_gamma_neg_1 * cc + r_gamma_pos_1 * cc_inv;
  }

  assert(gx.size() == 1);
  assert(x.size() == 1);
  assert(a.size() == 1);
  assert(y == x[0] * a[0]);
  assert(gamma == gx[0] * x[0] + gy * y + h * r_gamma);
  std::cout << "gx[0]: " << gx[0] << "\n";
  std::cout << "a[0]: " << a[0] << "\n";
  std::cout << "gamma: " << gamma << "\n";

  // final round
  com_ext_sec.d = FrRand();
  com_ext_sec.r_beta = FrRand();
  com_ext_sec.r_delta = FrRand();
  com_ext_pub.delta = gx[0] * com_ext_sec.d + h * com_ext_sec.r_delta;
  com_ext_pub.beta = gy * com_ext_sec.d + h * com_ext_sec.r_beta;
  UpdateSeed(seed, com_ext_pub.delta, com_ext_pub.beta);
  Fr c = H256ToFr(seed);
  std::cout << c << "\n";
  rom_proof.proof.z1 = com_ext_sec.d + c * y;
  rom_proof.proof.z2 =
      a[0] * (c * r_gamma + com_ext_sec.r_beta) + com_ext_sec.r_delta;
}

inline bool RomVerify(RomProof const& rom_proof, h256_t seed,
                      VerifierInput const& input) {
  assert(PcBase::kGSize >= rom_proof.aligned_n());
  auto n = input.a.size();
  if (!n || (int64_t)misc::Pow2UB(n) != rom_proof.aligned_n()) {
    return false;
  }
  CommitmentPub const& com_pub = input.com_pub;
  auto a = input.a;
  CommitmentExtPub const& com_ext_pub = rom_proof.com_ext_pub;
  int64_t round = (int64_t)misc::Log2UB(n);
  auto gx = GetPcBase().CopyG(input.x_g_offset, n);
  gx.resize(misc::Pow2UB(n));
  std::fill(gx.begin() + n, gx.end(), G1Zero());
  auto const& h = PcH();
  auto const& gy = PcG(input.y_g_offset);
  G1 gamma = com_pub.xi + com_pub.tau;

  UpdateSeed(seed, a, com_pub);

  // recursive round
  for (int64_t loop = 0; loop < round; ++loop) {
    std::vector<Fr> a1, a2;
    Divide(a, a1, a2, FrZero());
    std::vector<G1> g1, g2;
    Divide(gx, g1, g2, G1Zero());

    auto const& gamma_neg_1 = com_ext_pub.gamma_neg_1[loop];
    auto const& gamma_pos_1 = com_ext_pub.gamma_pos_1[loop];

    UpdateSeed(seed, gamma_neg_1, gamma_pos_1);
    Fr c = H256ToFr(seed);
    Fr cc = c * c;
    Fr c_inv = FrInv(c);
    Fr cc_inv = FrInv(cc);

    gamma += gamma_neg_1 * cc + gamma_pos_1 * cc_inv;
    a = a1 * c_inv + a2 * c;
    gx = FuncO(g1, c_inv, g2, c);
  }

  assert(gx.size() == 1);
  assert(a.size() == 1);
  std::cout << "gx[0]: " << gx[0] << "\n";
  std::cout << "a[0]: " << a[0] << "\n";
  std::cout << "gamma: " << gamma << "\n";

  // final round
  UpdateSeed(seed, com_ext_pub.delta, com_ext_pub.beta);
  Fr c = H256ToFr(seed);
  std::cout << c << "\n";
  auto const& proof = rom_proof.proof;
  G1 left = (gamma * c + com_ext_pub.beta) * a[0] + com_ext_pub.delta;
  G1 right = (gx[0] + gy * a[0]) * proof.z1 + h * proof.z2;
  if (left != right) {
    assert(false);
    return false;
  }
  return true;
}

inline bool TestRom(int64_t n) {
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  h256_t UpdateSeed = misc::RandH256();

  int64_t x_g_offset = 20;
  int64_t y_g_offset = -1;
  auto y = InnerProduct(x, a);
  ProverInput prover_input(x, a, y, x_g_offset, y_g_offset);

  CommitmentPub com_pub;
  CommitmentSec com_sec(FrRand(), FrRand());
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, UpdateSeed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(a, com_pub, x_g_offset, y_g_offset);
  bool success = RomVerify(rom_proof, UpdateSeed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace hyrax::a3

// TODO: func(A const& a); func(move(a)), what will happen?