#pragma once

#include <numeric>

#include "groth09/details.h"
#include "utils/fst.h"

// X, Y: secret vector<Fr>, size = n
// z: secret Fr
// open: com(gx, X), com(gy, Y), com(gz, z)
// prove: z = <X, Y>
// proof size: (2n+3) Fr and 4 G1
// prove cost: 2*mulexp(n)
// verify cost: mulexp(n)
namespace groth09::sec51a {

struct ProverInput {
  std::vector<Fr> const& x;  // size = n
  std::vector<Fr> const& y;  // size = n
  Fr const z;                // z = <x,y>
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& y, Fr const& z,
              int64_t x_g_offset, int64_t y_g_offset, int64_t z_g_offset)
      : x(x),
        y(y),
        z(z),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {
    assert(x.size() == y.size() && !x.empty());
    assert(z == InnerProduct(x, y));
  }

  int64_t n() const { return x.size(); }
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

struct CommitmentExtPub {
  CommitmentExtPub() {}
  CommitmentExtPub(G1 const& ad, G1 const& bd, G1 const& c1, G1 const& c0)
      : ad(ad), bd(bd), c1(c1), c0(c0) {}
  G1 ad;
  G1 bd;
  G1 c1;
  G1 c0;
};

inline bool operator==(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return left.ad == right.ad && left.bd == right.bd && left.c1 == right.c1 &&
         left.c0 == right.c0;
}

inline bool operator!=(CommitmentExtPub const& left,
                       CommitmentExtPub const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

struct CommitmentExtSec {
  std::vector<Fr> dx;
  std::vector<Fr> dy;
  Fr dz;
  Fr rd;
  Fr sd;
  Fr t1;
  Fr t0;
};

struct Proof {
  std::vector<Fr> fx;  // fx.size = n
  std::vector<Fr> fy;  // fy.size = n
  Fr rx;
  Fr sy;
  Fr tz;
  int64_t n() const { return (int64_t)fx.size(); }
};

inline bool operator==(Proof const& left, Proof const& right) {
  return left.fx == right.fx && left.fy == right.fy && left.rx == right.rx &&
         left.sy == right.sy && left.tz == right.tz;
}

inline bool operator!=(Proof const& left, Proof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("51.pf", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("51.pf", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

struct RomProof {
  CommitmentExtPub com_ext_pub;  // 4 G1
  Proof proof;                   // (2n+3) Fr
  bool CheckFormat() const {
    return true;  // TODO:
  }
  int64_t n() const { return proof.n(); }
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
  ar& YAS_OBJECT_NVP("51.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, RomProof& t) {
  ar& YAS_OBJECT_NVP("51.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  com_sec.r = FrRand();
  com_sec.s = FrRand();
  com_sec.t = FrRand();

  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_pub, &input, &com_sec]() {
    com_pub.a = PcComputeCommitmentG(input.x_g_offset, input.x, com_sec.r);
  };
  tasks[1] = [&com_pub, &input, &com_sec]() {
    com_pub.b = PcComputeCommitmentG(input.y_g_offset, input.y, com_sec.s);
  };
  tasks[2] = [&com_pub, &input, &com_sec]() {
    com_pub.c = PcComputeCommitmentG(input.z_g_offset, input.z, com_sec.t);
  };
  parallel::Invoke(tasks);
}

inline void ComputeComExt(CommitmentExtPub& com_ext_pub,
                          CommitmentExtSec& com_ext_sec,
                          ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  auto const n = input.n();
  com_ext_sec.dx.resize(n);
  FrRand(com_ext_sec.dx.data(), n);
  com_ext_sec.dy.resize(n);
  FrRand(com_ext_sec.dy.data(), n);
  com_ext_sec.dz = InnerProduct(com_ext_sec.dx, com_ext_sec.dy);

  com_ext_sec.rd = FrRand();
  com_ext_sec.sd = FrRand();
  com_ext_sec.t1 = FrRand();
  com_ext_sec.t0 = FrRand();

  Fr xdy_dxy = InnerProduct(input.x, com_ext_sec.dy) +
               InnerProduct(com_ext_sec.dx, input.y);

  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_ext_pub, &com_ext_sec, &input]() {
    com_ext_pub.ad =
        PcComputeCommitmentG(input.x_g_offset, com_ext_sec.dx, com_ext_sec.rd);
  };
  tasks[1] = [&com_ext_pub, &com_ext_sec, &input]() {
    com_ext_pub.bd =
        PcComputeCommitmentG(input.y_g_offset, com_ext_sec.dy, com_ext_sec.sd);
  };
  tasks[2] = [&com_ext_pub, &com_ext_sec, &xdy_dxy, &input]() {
    com_ext_pub.c1 =
        PcComputeCommitmentG(input.z_g_offset, xdy_dxy, com_ext_sec.t1);
    com_ext_pub.c0 =
        PcComputeCommitmentG(input.z_g_offset, com_ext_sec.dz, com_ext_sec.t0);
  };
  parallel::Invoke(tasks);
}

inline void ComputeProof(Proof& proof, ProverInput const& input,
                         CommitmentSec const& com_sec,
                         CommitmentExtSec const& com_ext_sec,
                         Fr const& challenge) {
  // Tick tick(__FUNCTION__);
  auto n = input.n();
  proof.fx.resize(n);
  proof.fy.resize(n);
  auto parallel_f = [&proof, &challenge, &input, &com_ext_sec](int64_t i) {
    // fx = e * x + dx
    proof.fx[i] = challenge * input.x[i] + com_ext_sec.dx[i];
    // fy = e * y + dy
    proof.fy[i] = challenge * input.y[i] + com_ext_sec.dy[i];
  };
  parallel::For(n, parallel_f, n < 16 * 1024);

  proof.rx = challenge * com_sec.r + com_ext_sec.rd;
  proof.sy = challenge * com_sec.s + com_ext_sec.sd;
  proof.tz = challenge * challenge * com_sec.t + challenge * com_ext_sec.t1 +
             com_ext_sec.t0;
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const com_pub,
                       CommitmentExtPub const& com_ext_pub) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.a);
  HashUpdate(hash, com_pub.b);
  HashUpdate(hash, com_pub.c);
  HashUpdate(hash, com_ext_pub.ad);
  HashUpdate(hash, com_ext_pub.bd);
  HashUpdate(hash, com_ext_pub.c1);
  HashUpdate(hash, com_ext_pub.c0);
  hash.Final(seed.data());
}

struct VerifierInput {
  VerifierInput(CommitmentPub const& com_pub, int64_t x_g_offset,
                int64_t y_g_offset, int64_t z_g_offset)
      : com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {}
  CommitmentPub const& com_pub;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
};

inline bool VerifyInternal(VerifierInput const& input, Fr const& challenge,
                           CommitmentExtPub const& com_ext_pub,
                           Proof const& proof) {
  auto const n = proof.fx.size();
  assert(n == proof.fy.size());
  (void)n;

  // std::cout << Tick::GetIndentString() << "multiexp(" << n << ")\n";

  auto const& com_pub = input.com_pub;
  auto const& e = challenge;
  std::vector<int64_t> rets;
  std::vector<parallel::Task> tasks;
  if (input.x_g_offset == input.y_g_offset) {
    // if x_g_offset == y_g_offset, we can combine the verification
    tasks.resize(2);
    rets.resize(2);
    tasks[0] = [&rets, &com_pub, &e, &com_ext_pub, &proof, &input]() {
      auto g_offset = input.x_g_offset;
      Fr alpha = FrRand();
      // (a^e * a_d)^alpha * b^e * b_d
      G1 left = (com_pub.a * e + com_ext_pub.ad) * alpha +
                (com_pub.b * e + com_ext_pub.bd);
      // com(alpha * fx + fy,alpha * rx + sy)
      std::vector<Fr> alpha_fx_fy = proof.fx * alpha + proof.fy;
      Fr alpha_rx_sy = alpha * proof.rx + proof.sy;
      G1 right = PcComputeCommitmentG(g_offset, alpha_fx_fy, alpha_rx_sy);
      rets[0] = left == right;
      assert(rets[0]);
    };
  } else {
    tasks.resize(3);
    rets.resize(3);
    tasks[0] = [&rets, &com_pub, &e, &com_ext_pub, &proof, &input]() {
      // a^e * a_d
      G1 left = (com_pub.a * e + com_ext_pub.ad);
      // com(fx,rx)
      G1 right = PcComputeCommitmentG(input.x_g_offset, proof.fx, proof.rx);
      rets[0] = left == right;
      assert(rets[0]);
    };
    tasks[1] = [&rets, &com_pub, &e, &com_ext_pub, &proof, &input]() {
      // b^e * b_d
      G1 left = com_pub.b * e + com_ext_pub.bd;
      // com(fy,sy)
      G1 right = PcComputeCommitmentG(input.y_g_offset, proof.fy, proof.sy);
      rets[1] = left == right;
      assert(rets[1]);
    };
  }

  tasks.back() = [&rets, &com_pub, &e, &com_ext_pub, &proof, &input]() {
    // c^(e^2) * c_1^e * c_0 == com(f_x * f_y , t_z)
    Fr e2_square = e * e;
    G1 left = com_pub.c * e2_square + com_ext_pub.c1 * e + com_ext_pub.c0;
    Fr fz = InnerProduct(proof.fx, proof.fy);

    G1 right = PcComputeCommitmentG(input.z_g_offset, fz, proof.tz);
    rets.back() = left == right;
    assert(rets.back());
  };

  parallel::Invoke(tasks);

  auto is_true = [](int64_t const& r) { return !!r; };
  auto all_success = std::all_of(rets.begin(), rets.end(), is_true);
  assert(all_success);
  return all_success;
}

inline void RomProve(RomProof& rom_proof, h256_t seed, ProverInput const& input,
                     CommitmentPub const& com_pub,
                     CommitmentSec const& com_sec) {
  // Tick tick(__FUNCTION__);

  assert(PcBase::kGSize >= input.n());

  CommitmentExtSec com_ext_sec;
  ComputeComExt(rom_proof.com_ext_pub, com_ext_sec, input);

  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  ComputeProof(rom_proof.proof, input, com_sec, com_ext_sec, challenge);
}

inline bool RomVerify(RomProof const& rom_proof, h256_t seed,
                      VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  assert(PcBase::kGSize >= rom_proof.n());

  UpdateSeed(seed, input.com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  return VerifyInternal(input, challenge, rom_proof.com_ext_pub,
                        rom_proof.proof);
}

inline bool TestRom(int64_t n) {
  Tick tick(__FUNCTION__);
  std::cout << "n=" << n << "\n";
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> y(n);
  FrRand(y.data(), n);
  Fr z = InnerProduct(x, y);
  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 5;
  int64_t y_g_offset = 10;
  int64_t z_g_offset = -1;
  ProverInput prover_input(x, y, z, x_g_offset, y_g_offset, z_g_offset);

  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(com_pub, x_g_offset, y_g_offset, z_g_offset);
  bool success = RomVerify(rom_proof, seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}

}  // namespace groth09::sec51a
