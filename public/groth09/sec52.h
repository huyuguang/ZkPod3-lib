#pragma once

#include <map>
#include <memory>

#include "utils/fst.h"
#include "./details.h"
#include "./sec51.h"

// t is a public vector which size = n
// x, y are secret matric which size = m*n
// commit {x1}, {y1}...., z
// prove z = <x1,y1 o t> + <x2,y2 o t>...
namespace groth09::sec52 {

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
  VerifierInput(std::vector<Fr> const* t, CommitmentPub const& com_pub)
      : t(t), com_pub(com_pub) {}
  std::vector<Fr> const* t;  // n, optional, can be null
  CommitmentPub const& com_pub;
  int64_t m() const { return com_pub.m(); }
  bool CheckFormat(int64_t check_n) const {
    if (t && t->size() != (uint64_t)check_n) return false;
    return com_pub.CheckFormat();
  }
};

class ProverInput {
 private:
  std::vector<std::vector<Fr>> yt_data_;

 private:
  std::vector<std::vector<Fr>> const* x_;   // m*n
  std::vector<std::vector<Fr>> const* y_;   // m*n
  std::vector<Fr> const* t_;                // n, optional, can be null
  std::vector<std::vector<Fr>> const* yt_;  // m*n
  Fr z_;

 private:
  std::vector<Fr> sum_xy_;  // m * 2 - 1

 public:
  std::vector<Fr> const* t() const { return t_; }
  Fr const& z() const { return z_; }
  std::vector<std::vector<Fr>> const& x() const { return *x_; }
  std::vector<std::vector<Fr>> const& y() const { return *y_; }
  std::vector<std::vector<Fr>> const& yt() const { return *yt_; }
  std::vector<Fr> const& x(size_t i) const { return x()[i]; }
  std::vector<Fr> const& y(size_t i) const { return y()[i]; }
  std::vector<Fr> const& yt(size_t i) const { return yt()[i]; }
  std::vector<Fr> const& sum_xy() const { return sum_xy_; }
  Fr const& sum_xy(size_t i) const { return sum_xy_[i]; }
  int64_t m() const { return x().size(); }
  int64_t n() const { return x(0).size(); }

 public:
  ProverInput(std::vector<std::vector<Fr>> const* x,
              std::vector<std::vector<Fr>> const* y, std::vector<Fr> const* t)
      : x_(x), y_(y), t_(t) {
    Tick tick(__FUNCTION__);
    if (!CheckFormat(false, false)) throw std::invalid_argument("");

    if (t_) {
      yt_data_.resize(m());
      for (int64_t i = 0; i < m(); ++i) {
        auto& yti = yt_data_[i];
        auto const& yi = (*y_)[i];
        details::HadamardProduct(yti, yi, *t_);
      }
      yt_ = &yt_data_;
    } else {
      yt_ = y_;
    }

    std::vector<Fr> temp_z(x->size());
    auto parallel_f = [this, &temp_z](int64_t i) {
      temp_z[i] = InnerProduct((*x_)[i], (*yt_)[i]);
    };
    parallel::For(m(), parallel_f);

    z_ = parallel::Accumulate(temp_z.begin(), temp_z.end(), FrZero());

    BuildSumXY();
  }

  ProverInput(std::vector<std::vector<Fr>> const* x,
              std::vector<std::vector<Fr>> const* y, std::vector<Fr> const* t,
              std::vector<std::vector<Fr>> const* yt, Fr const& z)
      : x_(x), y_(y), t_(t), yt_(yt), z_(z) {
    if (!CheckFormat(true, true)) throw std::invalid_argument("");
    BuildSumXY();
  }

  bool CheckFormat(bool check_yt, bool check_z) const {
#ifndef _DEBUG
    check_yt = false;
    check_z = false;
#endif

    if (!x_ || !y_ || x_->size() != y_->size() || x_->empty()) {
      assert(false);
      return false;
    }

    if (t_ && t_->size() != (size_t)n()) {
      assert(false);
      return false;
    }

    for (int64_t i = 0; i < m(); ++i) {
      if (x(i).size() != (size_t)n() || y(i).size() != (size_t)n()) {
        assert(false);
        return false;
      }
    }

    if (check_yt) {
      if (!t_ && yt_ != y_) {
        assert(false);
        return false;
      }

      if (t_) {
        for (int64_t i = 0; i < m(); ++i) {
          std::vector<Fr> check_value(n());
          details::HadamardProduct(check_value, y(i), *t_);
          if (yt(i) != check_value) {
            assert(false);
            return false;
          }
        }
      }
    }

    if (check_z) {
      std::vector<Fr> temp_z(m());
      auto parallel_f = [this, &temp_z](int64_t i) {
        temp_z[i] = InnerProduct(x(i), yt(i));
      };
      parallel::For(m(), parallel_f);

      auto check_value =
          parallel::Accumulate(temp_z.begin(), temp_z.end(), FrZero());
      if (check_value != z_) {
        assert(false);
        return false;
      }
    }

    return true;
  }

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
      ret += InnerProduct(x(i), yt(j));
    }
    return ret;
  }

  void BuildSumXY() {
    Tick tick(__FUNCTION__);
    sum_xy_.resize(m() * 2 - 1);

    auto parallel_f = [this](int64_t l) {
      sum_xy_[l] = ComputeSumOfXY(l);
    };
    parallel::For(m() * 2 - 1, parallel_f);
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
  sec51::RomProof rom_proof_51;
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
  auto const n = input.n();
  com_sec.r.resize(m);
  FrRand(com_sec.r.data(), m);

  com_sec.s.resize(m);
  FrRand(com_sec.s.data(), m);

  com_sec.t = FrRand();

  com_pub.a.resize(m);
  com_pub.b.resize(m);

  auto parallel_f = [&input, &com_pub, &com_sec](int64_t i) {
    com_pub.a[i] = PcComputeCommitment(input.x(i), com_sec.r[i]);
    com_pub.b[i] = PcComputeCommitment(input.y(i), com_sec.s[i]);
  };
  parallel::For(m, parallel_f);

  com_pub.c = PcComputeCommitment(input.z(), com_sec.t);
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
    auto const& sum_xy = input.sum_xy(l);
    com_ext_pub.cl[l] = PcComputeCommitment(sum_xy, com_ext_sec.tl[l]);
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
      x[j] += e_pow[i] * input.x(i)[j];
    }
  };
  parallel::For(n, parallel_fx);

  std::vector<Fr> y(n);
  std::fill(y.begin(), y.end(), FrZero());

  auto parallel_fy = [m, &y, &e_pow_reverse, &input](int64_t j) {
    for (int64_t i = 0; i < m; ++i) {
      y[j] += e_pow_reverse[i] * input.y(i)[j];
    }
  };
  parallel::For(n, parallel_fy);

  sec51::ProverInput input_51(&x, &y, input.t());

  assert(input_51.z() == InnerProduct(input.sum_xy(), e_pow));

  sec51::CommitmentPub com_pub_51;
  sec51::CommitmentSec com_sec_51;
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
  auto check_a = PcComputeCommitment(input_51.x(), com_sec_51.r);
  assert(check_a == com_pub_51.a);
  auto check_b = PcComputeCommitment(input_51.y(), com_sec_51.s);
  assert(check_b == com_pub_51.b);
#endif

  sec51::RomProve(rom_proof.rom_proof_51, seed, input_51, com_pub_51,
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

  sec51::CommitmentPub com_pub_51;
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

  sec51::VerifierInput input_51(input.t, com_pub_51);
  return sec51::RomVerify(rom_proof.rom_proof_51, seed, input_51);
}

inline bool TestRom(int64_t m, int64_t n) {
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

  ProverInput prover_input(&x, &y, &t);
  CommitmentPub com_pub;
  CommitmentSec com_sec;
  ComputeCom(com_pub, com_sec, prover_input);

  RomProof rom_proof;
  RomProve(rom_proof, common_seed, prover_input, com_pub, com_sec);

  VerifierInput verifier_input(&t, com_pub);
  return RomVerify(rom_proof, common_seed, verifier_input);
}
}  // namespace groth09::sec52