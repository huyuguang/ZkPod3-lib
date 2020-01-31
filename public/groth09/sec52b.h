#pragma once

#include <map>
#include <memory>

#include "./details.h"
#include "./sec51b.h"
#include "utils/fst.h"

// t is a public vector which size = n
// x, y are secret matric which size = m*n
// commit {x1}, {y1}...., z
// open: com(gx, x1), com(gx, x2)...com(gy, y1)... com(gz, z)
// prove z = <x1,y1 o t> + <x2,y2 o t>...
namespace groth09::sec52b {

struct CommitmentPub {
  std::vector<G1> a;  // a.size = m
  std::vector<G1> b;  // b.size = m
  G1 c;
  int64_t m() const { return a.size(); }
  bool CheckFormat() const { return !a.empty() && (a.size() == b.size()); }
};

struct CommitmentSec {
  std::vector<Fr> r;  // r.size = m
  std::vector<Fr> s;  // s.size = m
  Fr t;
};

struct VerifierInput {
  VerifierInput(std::vector<Fr> const& t, CommitmentPub const& com_pub,
                int64_t x_g_offset, int64_t y_g_offset, int64_t z_g_offset)
      : t(t),
        com_pub(com_pub),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset) {}
  std::vector<Fr> const& t;
  CommitmentPub const& com_pub;
  int64_t m() const { return com_pub.m(); }
  bool CheckFormat(int64_t check_n) const {
    if (t.size() != (uint64_t)check_n) return false;
    return com_pub.CheckFormat();
  }
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
};

struct ProverInput {
  std::vector<std::vector<Fr>> const& x;   // m*n
  std::vector<std::vector<Fr>> const& y;   // m*n
  std::vector<Fr> const& t;                // n
  std::vector<std::vector<Fr>> const& yt;  // m*n
  Fr const z;
  int64_t const x_g_offset;
  int64_t const y_g_offset;
  int64_t const z_g_offset;
  std::vector<Fr> sum_xy;  // m * 2 - 1

  int64_t m() const { return x.size(); }
  int64_t n() const { return x[0].size(); }

  ProverInput(std::vector<std::vector<Fr>> const& x,
              std::vector<std::vector<Fr>> const& y, std::vector<Fr> const& t,
              std::vector<std::vector<Fr>> const& yt, Fr const& z,
              int64_t x_g_offset, int64_t y_g_offset, int64_t z_g_offset)
      : x(x),
        y(y),
        t(t),
        yt(yt),
        z(z),
        x_g_offset(x_g_offset),
        y_g_offset(y_g_offset),
        z_g_offset(z_g_offset),
        sum_xy(BuildSumXY()) {
#ifdef _DEBUG
    assert(x.size() == y.size() && !x.empty());
    assert(x.size() == yt.size());
    assert(t.size() == (size_t)n());
    Fr check_z = FrZero();
    for (int64_t i = 0; i < m(); ++i) {
      assert(x[i].size() == (size_t)n());
      assert(y[i].size() == (size_t)n());
      assert(yt[i].size() == (size_t)n());
      assert(yt[i] == details::HadamardProduct(y[i], t));
      check_z += InnerProduct(x[i], yt[i]);
    }
    assert(z == check_z);
#endif
  }

 private:
  Fr ComputeSumOfXY(int64_t l) {
    Fr ret = FrZero();
    // l = m + i - j - 1
    int64_t min_j = m() - l - 1;
    int64_t max_j = m() * 2 - l - 1;
    if (min_j < 0) min_j = 0;
    if (max_j > m()) max_j = m();
    for (int64_t j = min_j; j < max_j; ++j) {
      int64_t i = l + j + 1LL - m();
      assert(i >= 0 && i < m());
      ret += InnerProduct(x[i], yt[j]);
    }
    return ret;
  }

  std::vector<Fr> BuildSumXY() {
    Tick tick(__FUNCTION__);
    std::vector<Fr> ret(m() * 2 - 1);

    auto parallel_f = [this, &ret](int64_t l) { ret[l] = ComputeSumOfXY(l); };
    parallel::For(m() * 2 - 1, parallel_f);
    return ret;
  }
};

struct CommitmentExtPub {
  std::vector<G1> cl;  // cl.size = 2m-1
  bool CheckFormat(int64_t check_m) const {
    return (int64_t)cl.size() == check_m * 2 - 1;
  }
  int64_t m() const { return (cl.size() + 1) / 2; }
};

struct CommitmentExtSec {
  std::vector<Fr> tl;  // tl.size = 2m
};

struct RomProof {
  CommitmentExtPub com_ext_pub;
  sec51b::RomProof rom_proof_51;
  int64_t n() const { return rom_proof_51.n(); }
  int64_t m() const { return com_ext_pub.m(); }

  bool CheckFormat(int64_t check_m) const {
    return com_ext_pub.CheckFormat(check_m) && rom_proof_51.CheckFormat();
  }
};

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  Tick tick(__FUNCTION__);
  auto const m = input.m();
  // auto const n = input.n();
  com_sec.r.resize(m);
  FrRand(com_sec.r.data(), m);

  com_sec.s.resize(m);
  FrRand(com_sec.s.data(), m);

  com_sec.t = FrRand();

  com_pub.a.resize(m);
  com_pub.b.resize(m);

  auto parallel_f = [&input, &com_pub, &com_sec](int64_t i) {
    com_pub.a[i] =
        PcComputeCommitmentG(input.x_g_offset, input.x[i], com_sec.r[i]);
    com_pub.b[i] =
        PcComputeCommitmentG(input.y_g_offset, input.y[i], com_sec.s[i]);
  };
  parallel::For(m, parallel_f);

  com_pub.c = PcComputeCommitmentG(input.z_g_offset, input.z, com_sec.t);
}

inline void ComputeComExt(CommitmentExtPub& com_ext_pub,
                          CommitmentExtSec& com_ext_sec,
                          ProverInput const& input,
                          CommitmentPub const& com_pub,
                          CommitmentSec const& com_sec) {
  Tick tick(__FUNCTION__);

  auto const m = input.m();
  com_ext_sec.tl.resize(m * 2);
  FrRand(com_ext_sec.tl.data(), m * 2);
  com_ext_sec.tl[m - 1] = com_sec.t;

  com_ext_pub.cl.resize(m * 2 - 1);

  auto parallel_f = [&input, &com_ext_pub, &com_ext_sec](int64_t l) {
    com_ext_pub.cl[l] = PcComputeCommitmentG(input.z_g_offset, input.sum_xy[l],
                                             com_ext_sec.tl[l]);
  };
  parallel::For(m * 2 - 1, parallel_f);

  if (com_ext_pub.cl[m - 1] != com_pub.c) {
    assert(false);
    throw std::runtime_error(__FUNCTION__);
  }
}

// compute challenge1 by commitment1 and commitment2
inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub,
                       CommitmentExtPub const& com_ext_pub) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.a);
  HashUpdate(hash, com_pub.b);
  HashUpdate(hash, com_pub.c);
  HashUpdate(hash, com_ext_pub.cl);
  hash.Final(seed.data());
}

inline void RomProve(RomProof& rom_proof, h256_t const& common_seed,
                     ProverInput const& input, CommitmentPub const& com_pub,
                     CommitmentSec const& com_sec) {
  Tick tick(__FUNCTION__);
  auto m = input.m();
  auto n = input.n();
  assert(PcBase::kGSize >= input.n());

  CommitmentExtSec com_ext_sec;
  ComputeComExt(rom_proof.com_ext_pub, com_ext_sec, input, com_pub, com_sec);

  auto seed = common_seed;
  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr e = H256ToFr(seed);

  std::vector<Fr> e_pow;
  std::vector<Fr> e_pow_reverse;
  details::ComputePowOfE(e, m, e_pow, e_pow_reverse);

  std::vector<Fr> x(n);
  std::fill(x.begin(), x.end(), FrZero());

  auto parallel_fx = [m, &x, &e_pow, &input](int64_t j) {
    for (int64_t i = 0; i < m; ++i) {
      x[j] += e_pow[i] * input.x[i][j];
    }
  };
  parallel::For(n, parallel_fx);

  std::vector<Fr> y(n);
  std::fill(y.begin(), y.end(), FrZero());

  auto parallel_fy = [m, &y, &e_pow_reverse, &input](int64_t j) {
    for (int64_t i = 0; i < m; ++i) {
      y[j] += e_pow_reverse[i] * input.y[i][j];
    }
  };
  parallel::For(n, parallel_fy);

  auto const& t = input.t;
  std::vector<Fr> yt = details::HadamardProduct(y, t);
  auto z = InnerProduct(x, yt);
  sec51b::ProverInput input_51(x, y, t, yt, z, input.x_g_offset,
                               input.y_g_offset, input.z_g_offset);

  assert(z == InnerProduct(input.sum_xy, e_pow));

  sec51b::CommitmentPub com_pub_51;
  sec51b::CommitmentSec com_sec_51;
  com_pub_51.a = MultiExpBdlo12(com_pub.a, e_pow);
  com_pub_51.b = MultiExpBdlo12(com_pub.b, e_pow_reverse);
  com_pub_51.c = MultiExpBdlo12(rom_proof.com_ext_pub.cl, e_pow);

  com_sec_51.r = FrZero();
  for (int64_t i = 0; i < m; ++i) {
    com_sec_51.r += e_pow[i] * com_sec.r[i];
  }

  com_sec_51.s = FrZero();
  for (int64_t i = 0; i < m; ++i) {
    com_sec_51.s += e_pow_reverse[i] * com_sec.s[i];
  }

  com_sec_51.t = FrZero();
  for (int64_t l = 0; l < m * 2 - 1; ++l) {
    com_sec_51.t += com_ext_sec.tl[l] * e_pow[l];
  }

#ifdef _DEBUG
  auto check_a =
      PcComputeCommitmentG(input_51.x_g_offset, input_51.x, com_sec_51.r);
  assert(check_a == com_pub_51.a);
  auto check_b =
      PcComputeCommitmentG(input_51.y_g_offset, input_51.y, com_sec_51.s);
  assert(check_b == com_pub_51.b);
  auto check_c =
      PcComputeCommitmentG(input_51.z_g_offset, input_51.z, com_sec_51.t);
  assert(check_c == com_pub_51.c);
#endif

  sec51b::RomProve(rom_proof.rom_proof_51, seed, input_51, com_pub_51,
                   com_sec_51);
}

inline bool RomVerify(RomProof const& rom_proof, h256_t const& common_seed,
                      VerifierInput const& input) {
  auto m = rom_proof.m();
  assert(PcBase::kGSize >= rom_proof.n());

  auto const& com_pub = input.com_pub;
  auto const& com_ext_pub = rom_proof.com_ext_pub;

  auto seed = common_seed;
  UpdateSeed(seed, com_pub, com_ext_pub);
  Fr e = H256ToFr(seed);

  std::vector<Fr> e_pow;
  std::vector<Fr> e_pow_reverse;
  details::ComputePowOfE(e, m, e_pow, e_pow_reverse);

  std::cout << Tick::GetIndentString() << "2 multiexp(" << m << ")\n";
  std::cout << Tick::GetIndentString() << "multiexp(" << 2 * m - 1 << ")\n";

  sec51b::CommitmentPub com_pub_51;
  std::array<parallel::Task, 3> tasks;
  tasks[0] = [&com_pub_51, &com_pub, e_pow]() {
    com_pub_51.a = MultiExpBdlo12(com_pub.a, e_pow);
  };
  tasks[1] = [&com_pub_51, &com_pub, e_pow_reverse]() {
    com_pub_51.b = MultiExpBdlo12(com_pub.b, e_pow_reverse);
  };
  tasks[2] = [&com_pub_51, &com_ext_pub, e_pow]() {
    com_pub_51.c = MultiExpBdlo12(com_ext_pub.cl, e_pow);
  };
  parallel::Invoke(tasks);

  sec51b::VerifierInput input_51(input.t, com_pub_51, input.x_g_offset,
                                 input.y_g_offset, input.z_g_offset);
  return sec51b::RomVerify(rom_proof.rom_proof_51, seed, input_51);
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

  std::vector<Fr> t(n);
  FrRand(t.data(), t.size());

  h256_t common_seed = misc::RandH256();

  int64_t x_g_offset = 10;
  int64_t y_g_offset = 550;
  int64_t z_g_offset = -1;

  std::vector<std::vector<Fr>> yt(m);
  Fr z = FrZero();
  for (int64_t i = 0; i < m; ++i) {
    yt[i] = details::HadamardProduct(y[i], t);
    z += InnerProduct(x[i], yt[i]);
  }

  ProverInput prover_input(x, y, t, yt, z, x_g_offset, y_g_offset, z_g_offset);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, common_seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(t, com_pub, x_g_offset, y_g_offset, z_g_offset);
  bool success = RomVerify(rom_proof, common_seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace groth09::sec52b