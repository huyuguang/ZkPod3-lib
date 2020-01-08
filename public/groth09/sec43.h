#pragma once

#include "./details.h"
#include "./sec53.h"

// x, y, z: secret matrix<Fr>, size =m*n
// open: a1=com(x1)...am=com(xm)
// open b1=com(y1)...bm=com(ym)
// open c1=com(z1)...cm=com(zm)
// prove: z=x o y
// proof size: 2*log(m)+6 G1, 3n+5 Fr
// prove cost: 2*log(m)*mulexp(n)
// verify cost: 2*mulexp(n)
namespace groth09::sec43 {

struct CommitmentPub {
  std::vector<G1> a;  // a.size = m
  std::vector<G1> b;  // b.size = m
  std::vector<G1> c;  // c.size = m
  void Align() {
    int64_t old_m = a.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      static const G1 g0 = G1Zero();
      a.resize(new_m);
      std::fill(a.begin() + old_m, a.end(), g0);
      b.resize(new_m);
      std::fill(b.begin() + old_m, b.end(), g0);
      c.resize(new_m);
      std::fill(c.begin() + old_m, c.end(), g0);
    }
  }
};

struct CommitmentSec {
  std::vector<Fr> r;  // r.size = m
  std::vector<Fr> s;  // s.size = m
  std::vector<Fr> t;  // t.size = m
  void Align() {
    int64_t old_m = r.size();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (new_m > old_m) {
      static const Fr f0 = FrZero();
      r.resize(new_m);
      std::fill(r.begin() + old_m, r.end(), f0);
      s.resize(new_m);
      std::fill(s.begin() + old_m, s.end(), f0);
      t.resize(new_m);
      std::fill(t.begin() + old_m, t.end(), f0);
    }
  }
};

struct RomProof {
  G1 c;
  sec53::RomProof proof_53;      // 2*log(m)+4 G1, 2n+3 Fr
  hyrax::a2::RomProof proof_a2;  // 2 G1, n+2 Fr
  int64_t n() const { return proof_53.n(); }
  int64_t m() const { return proof_53.m(); }
};

inline bool operator==(RomProof const& left, RomProof const& right) {
  return left.c == right.c && left.proof_53 == right.proof_53 &&
         left.proof_a2 == right.proof_a2;
}
inline bool operator!=(RomProof const& left, RomProof const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, RomProof const &t) {
  ar &YAS_OBJECT_NVP("43.rp", ("c", t.c), ("53p", t.proof_53),
                     ("a2p", t.proof_a2));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, RomProof &t) {
  ar &YAS_OBJECT_NVP("43.rp", ("c", t.c), ("53p", t.proof_53),
                     ("a2p", t.proof_a2));
}

class ProverInput {
 private:
  std::vector<std::vector<Fr>> x_;  // m*n
  std::vector<std::vector<Fr>> y_;
  std::vector<std::vector<Fr>> z_;

 public:
  int64_t m() const { return x_.size(); }
  int64_t n() const { return x_[0].size(); }
  std::vector<std::vector<Fr>> const& x() const { return x_; }
  std::vector<Fr> const& x(size_t i) const { return x_[i]; }
  std::vector<std::vector<Fr>> const& y() const { return y_; }
  std::vector<Fr> const& y(size_t i) const { return y_[i]; }
  std::vector<std::vector<Fr>> const& z() const { return z_; }
  std::vector<Fr> const& z(size_t i) const { return z_[i]; }
  void Take(std::vector<std::vector<Fr>>& x, std::vector<std::vector<Fr>>& y,
            std::vector<std::vector<Fr>>& z) {
    x = std::move(x_);
    y = std::move(y_);
    z = std::move(z_);
  }

 public:
  ProverInput(std::vector<std::vector<Fr>> x, std::vector<std::vector<Fr>> y,
              std::vector<std::vector<Fr>> z)
      : x_(std::move(x)), y_(std::move(y)), z_(std::move(z)) {
    // Tick tick(__FUNCTION__);
    assert(!x_.empty());
    assert(x_.size() == y_.size());
    assert(x_.size() == z_.size());
    for (size_t i = 0; i < x_.size(); ++i) {
      assert(x_[i].size() == x_[0].size());
      assert(y_[i].size() == x_[0].size());
      assert(z_[i].size() == x_[0].size());
    }
  }

  // pad some trivial value
  void Align() {
    // Tick tick(__FUNCTION__);
    int64_t old_m = m();
    int64_t new_m = (int64_t)misc::Pow2UB(old_m);
    if (old_m == new_m) return;

    Fr f0 = FrZero();
    x_.resize(new_m);
    y_.resize(new_m);
    z_.resize(new_m);
    for (int64_t i = old_m; i < new_m; ++i) {
      auto& xi = x_[i];
      xi.resize(n());
      std::fill(xi.begin(), xi.end(), f0);
      auto& yi = y_[i];
      yi.resize(n());
      std::fill(yi.begin(), yi.end(), f0);
      auto& zi = z_[i];
      zi.resize(n());
      std::fill(zi.begin(), zi.end(), f0);
    }
  }
};

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                       ProverInput const& input) {
  // Tick tick(__FUNCTION__);
  auto const m = input.m();
  com_sec.r.resize(m);
  FrRand(com_sec.r.data(), m);

  com_sec.s.resize(m);
  FrRand(com_sec.s.data(), m);

  com_sec.t.resize(m);
  FrRand(com_sec.t.data(), m);

  com_pub.a.resize(m);
  com_pub.b.resize(m);
  com_pub.c.resize(m);

  auto parallel_f = [&com_sec, &com_pub, &input](int64_t i) mutable {
    auto const& r = com_sec.r[i];
    auto const& s = com_sec.s[i];
    auto const& t = com_sec.t[i];
    auto const& x = input.x(i);
    auto const& y = input.y(i);
    auto const& z = input.z(i);
    com_pub.a[i] = PcComputeCommitment(x, r);
    com_pub.b[i] = PcComputeCommitment(y, s);
    com_pub.c[i] = PcComputeCommitment(z, t);
  };
  parallel::For(m, parallel_f);
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const& com_pub) {
  // Tick tick(__FUNCTION__);
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, com_pub.a);
  HashUpdate(hash, com_pub.b);
  HashUpdate(hash, com_pub.c);
  hash.Final(seed.data());
}

inline void ComputeChallengeKT(h256_t const& seed, std::vector<Fr>& k,
                               std::vector<Fr>& t) {
  ComputeFst(seed, "gro09::sec43::k", k);
  ComputeFst(seed, "gro09::sec43::t", t);
}

// pad some trivial value
inline void AlignData(ProverInput& input, CommitmentPub& com_pub,
                      CommitmentSec& com_sec) {
  // Tick tick(__FUNCTION__);
  input.Align();
  com_pub.Align();
  com_sec.Align();
}

inline void RomProve(RomProof& rom_proof, h256_t const& rom_seed,
                     ProverInput input, CommitmentPub com_pub,
                     CommitmentSec com_sec) {
  // Tick tick(__FUNCTION__);
  auto m = input.m();
  auto n = input.n();

  auto seed = rom_seed;
  UpdateSeed(seed, com_pub);
  std::vector<Fr> k(m);
  std::vector<Fr> t(n);
  ComputeChallengeKT(seed, k, t);

  Fr input_53_z;
  G1 com_pub_53_c;
  Fr com_sec_53_t;

  std::vector<std::vector<Fr>> input_x;
  std::vector<std::vector<Fr>> input_y;
  std::vector<std::vector<Fr>> input_z;
  input.Take(input_x, input_y, input_z);

  {
    // Tick tick53(" prepare for sec53");

    sec53::CommitmentSec com_sec_53;
    sec53::CommitmentPub com_pub_53;
    std::vector<std::vector<Fr>> x_53(m);

    auto parallel_f = [&input_x, &k](int64_t i) { input_x[i] *= k[i]; };
    parallel::For(m, parallel_f, m < 1024);

    sec53::ProverInput input_53(std::move(input_x), std::move(input_y), &t);
    input_53_z = input_53.z();

    com_sec_53.r.resize(m);
    com_pub_53.a.resize(m);
    auto parallel_f2 = [&com_sec, &com_pub, &com_sec_53, &com_pub_53,
                        &k](int64_t i) {
      com_sec_53.r[i] = com_sec.r[i] * k[i];
      com_pub_53.a[i] = com_pub.a[i] * k[i];
    };
    parallel::For(m, parallel_f2, m < 16 * 1024);

    com_sec_53.s = com_sec.s;
    com_sec_53.t = FrRand();
    com_sec_53_t = com_sec_53.t;

    com_pub_53.b = com_pub.b;
    com_pub_53.c = PcComputeCommitment(input_53_z, com_sec_53.t);
    rom_proof.c = com_pub_53.c;  // verifier can not compute c by com_pub.c
    com_pub_53_c = com_pub_53.c;

    sec53::RomProve(rom_proof.proof_53, seed, std::move(input_53),
                    std::move(com_pub_53), std::move(com_sec_53));
  }

  {
    // Tick tick53("prepare for hyrax");
    hyrax::a2::CommitmentPub com_pub_hy;
    hyrax::a2::CommitmentSec com_sec_hy;

    std::vector<Fr> x_hy(n);

    std::fill(x_hy.begin(), x_hy.end(), FrZero());

    auto parallel_f = [&input_z, &x_hy, &k, m](int64_t j) mutable {
      for (int64_t i = 0; i < m; ++i) {
        x_hy[j] += input_z[i][j] * k[i];
      }
    };
    parallel::For(n, parallel_f, n < 16 * 1024);

    hyrax::a2::ProverInput input_hy(x_hy, t, input_53_z);

    com_sec_hy.r_xi = InnerProduct(com_sec.t, k);
    com_sec_hy.r_tau = com_sec_53_t;

    com_pub_hy.tau = com_pub_53_c;

    // do not need to compute the com1_pub_hy.xi in release build
    assert(input_hy.y == InnerProduct(input_hy.x, input_hy.a));

    com_pub_hy.xi = MultiExpBdlo12(com_pub.c, k);
#ifdef _DEBUG
    auto check_xi = PcComputeCommitment(input_hy.x, com_sec_hy.r_xi);
    assert(check_xi == com_pub_hy.xi);
#endif

    hyrax::a2::RomProve(rom_proof.proof_a2, seed, std::move(input_hy),
                        std::move(com_pub_hy), std::move(com_sec_hy));
  }
}

struct VerifierInput {
  VerifierInput(CommitmentPub const& com_pub) : com_pub(com_pub) {}
  CommitmentPub const& com_pub;
};

inline bool RomVerify(RomProof const& rom_proof, h256_t const& rom_seed,
                      VerifierInput const& input) {
  // Tick tick(__FUNCTION__);
  auto m = rom_proof.m();
  auto n = rom_proof.n();

  auto const& com_pub = input.com_pub;
  auto seed = rom_seed;
  UpdateSeed(seed, com_pub);
  std::vector<Fr> k(m);
  std::vector<Fr> t(n);
  ComputeChallengeKT(seed, k, t);

  std::array<parallel::Task, 2> tasks;
  bool ret_53 = false;
  tasks[0] = [&ret_53, &rom_proof, &input, m, &com_pub, &k, &t,
              &seed]() mutable {
    sec53::CommitmentPub com_pub_53;
    com_pub_53.c = rom_proof.c;
    com_pub_53.b = input.com_pub.b;
    com_pub_53.a.resize(m);
    auto parallel_f = [&com_pub_53, &com_pub, &k](int64_t i) {
      com_pub_53.a[i] = com_pub.a[i] * k[i];
    };
    parallel::For(m, parallel_f, m < 1024);

    sec53::VerifierInput input_53(&t, com_pub_53);
    ret_53 = sec53::RomVerify(rom_proof.proof_53, seed, input_53);
  };

  bool ret_a2 = false;
  tasks[1] = [&ret_a2, &com_pub, &rom_proof, &t, &k, &seed]() {
    hyrax::a2::CommitmentPub com_pub_hy(MultiExpBdlo12(com_pub.c, k),
                                        rom_proof.c);
    hyrax::a2::VerifierInput input_hy(t, com_pub_hy);
    ret_a2 = hyrax::a2::RomVerify(rom_proof.proof_a2, seed, input_hy);
  };

  parallel::Invoke(tasks);

  if (!ret_53 || !ret_a2) {
    std::cout << "ret_53: " << ret_53 << ", ret_a2: " << ret_a2 << "\n";
    assert(false);
  }
  return ret_53 && ret_a2;
}

inline bool TestRom(int64_t m, int64_t n) {
  std::cout << "old_m=" << m << ", n=" << n << "\n";

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

  std::vector<std::vector<Fr>> z(m);
  for (int64_t i = 0; i < m; ++i) {
    z[i] = details::HadamardProduct(x[i], y[i]);
  }

  h256_t common_seed = misc::RandH256();

  ProverInput prover_input(x, y, z);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  AlignData(prover_input, com_pub, com_sec);

  RomProof rom_proof;
  RomProve(rom_proof, common_seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(com_pub);
  return RomVerify(rom_proof, common_seed, verifier_input);
}
}  // namespace groth09::sec43
