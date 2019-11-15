#pragma once

#include <numeric>
#include "groth09/details.h"

// t is a public vector
// commit {x},{y},z
// prove z = <x, y o t>
namespace groth09::sec51 {

class ProverInput {
 private:
  std::vector<Fr> yt_data_;  // size = n or empty
  std::vector<Fr> const* x_;   // size = n
  std::vector<Fr> const* y_;   // size = n
  std::vector<Fr> const* t_;   // size = n or empty
  std::vector<Fr> const* yt_;  // if (!yt_data_.empty()), yt_=&yt_data_  
  Fr z_;

 public:
  std::vector<Fr> const* t() const { return t_; }
  Fr const& z() const { return z_; }
  std::vector<Fr> const& x() const { return *x_; }
  std::vector<Fr> const& y() const { return *y_; }
  std::vector<Fr> const& yt() const { return *yt_; }
  Fr const& x(size_t i) const { return x()[i]; }
  Fr const& y(size_t i) const { return y()[i]; }
  Fr const& yt(size_t i) const { return yt()[i]; }
  int64_t n() const { return x().size(); }

 public:
  ProverInput(std::vector<Fr> const* x,
              std::vector<Fr> const* y,
              std::vector<Fr> const* t)
      : x_(x), y_(y), t_(t) {
    if (!CheckFormat(false, false)) throw std::invalid_argument("");
    if (t_) {
      yt_data_.resize(n());
      details::HadamardProduct(yt_data_, *y_, *t_);
      yt_ = &yt_data_;
    } else {
      yt_ = y_;
    }
    z_ = InnerProduct(*x_, *yt_);
  }

  ProverInput(std::vector<Fr> const* x,
              std::vector<Fr> const* y,
              std::vector<Fr> const* t,
              std::vector<Fr> const* yt, Fr const& iz)
      : x_(x),
        y_(y),
        t_(t),
        yt_(yt),
        z_(iz) {
    if (!CheckFormat(true, true)) throw std::invalid_argument("");
  }

  ProverInput(ProverInput const& o)
      : yt_data_(o.yt_data_),
        x_(o.x_),
        y_(o.y_),
        t_(o.t_),
        yt_(o.yt_),        
        z_(o.z_) {
    if (!yt_data_.empty()) {
      yt_ = &yt_data_;
    }
  }

  ProverInput(ProverInput&& o)
      : yt_data_(std::move(o.yt_data_)),
        x_(o.x_),
        y_(o.y_),
        t_(o.t_),
        yt_(o.yt_),
        z_(std::move(o.z_)) {
    if (!yt_data_.empty()) {
      yt_ = &yt_data_;
    }
  }

  ProverInput& operator=(ProverInput const& o) {
    yt_data_ = o.yt_data_;
    x_ = o.x_;
    y_ = o.y_;
    t_ = o.t_;
    yt_ = o.yt_;
    z_ = o.z_;
    if (!yt_data_.empty()) {
      yt_ = &yt_data_;
    }
    return *this;
  }

  ProverInput& operator=(ProverInput&& o) {
    yt_data_ = std::move(o.yt_data_);
    x_ = o.x_;
    y_ = o.y_;
    t_ = o.t_;
    yt_ = o.yt_;
    z_ = o.z_;
    if (!yt_data_.empty()) {
      yt_ = &yt_data_;
    }
    return *this;
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

    if (check_yt) {
      if (!t_ && yt_ != y_) {
        assert(false);
        return false;
      }

      if (t_) {
        std::vector<Fr> check_value(n());
        details::HadamardProduct(check_value, y(), *t_);
        if (yt() != check_value) {
          assert(false);
          return false;
        }
      }
    }

    if (check_z) {
      Fr check_value = InnerProduct(x(), yt());
      if (check_value != z_) {
        assert(false);
        return false;
      }
    }

    return true;
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

struct CommitmentExtSec {
  std::vector<Fr> dx;
  std::vector<Fr> dy;
  std::vector<Fr> dyt;  // dy o t
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

inline bool operator==(Proof const& left,
                       Proof const& right) {
  return left.fx == right.fx && left.fy == right.fy && left.rx == right.rx &&
         left.sy == right.sy && left.tz == right.tz;
}

inline bool operator!=(Proof const& left,
                       Proof const& right) {
  return !(left == right);
}

struct RomProof {
  CommitmentExtPub com_ext_pub;
  Proof proof;
  bool CheckFormat() const {
    return true;  // TODO:
  }
  int64_t n() const { return proof.n(); }
};

inline bool operator==(RomProof const& left,
                       RomProof const& right) {
  return left.com_ext_pub == right.com_ext_pub && left.proof == right.proof;
}

inline bool operator!=(RomProof const& left,
                       RomProof const& right) {
  return !(left == right);
}

inline void ComputeCom(CommitmentPub& com_pub, CommitmentSec& com_sec,
                ProverInput const& input) {
  //Tick tick(__FUNCTION__);
  com_sec.r = FrRand();
  com_sec.s = FrRand();
  com_sec.t = FrRand();

  using details::ComputeCommitment;

  //std::cout << Tick::GetIndentString() << "2*multiexp(" << input.n() << ")\n";

  std::vector<parallel::Task> tasks(3);
  tasks[0] = [&com_pub,&input,&com_sec]() mutable {
    com_pub.a = ComputeCommitment(input.x(), com_sec.r);
  };
  tasks[1] = [&com_pub, &input, &com_sec]() mutable {
    com_pub.b = ComputeCommitment(input.y(), com_sec.s);
  };
  tasks[2] = [&com_pub, &input, &com_sec]() mutable {
    com_pub.c = ComputeCommitment(input.z(), com_sec.t);
  };
  parallel::Invoke(tasks);

//
////#ifdef MULTICORE
////#pragma omp parallel sections
////#endif
//  {
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_pub.a = ComputeCommitment(input.x(), com_sec.r);
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_pub.b = ComputeCommitment(input.y(), com_sec.s);
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_pub.c = ComputeCommitment(input.z(), com_sec.t);
//  }
}

inline void ComputeComExt(CommitmentExtPub& com_ext_pub, CommitmentExtSec& com_ext_sec,
                   ProverInput const& input) {
  //Tick tick(__FUNCTION__);
  auto const n = input.n();
  com_ext_sec.dx.resize(n);
  FrRand(com_ext_sec.dx.data(), n);
  com_ext_sec.dy.resize(n);
  FrRand(com_ext_sec.dy.data(), n);
  com_ext_sec.dyt = com_ext_sec.dy;
  if (input.t()) {
    com_ext_sec.dyt.resize(n);
    details::HadamardProduct(com_ext_sec.dyt, com_ext_sec.dy, *input.t());
  }
  com_ext_sec.dz = InnerProduct(com_ext_sec.dx, com_ext_sec.dyt);

  com_ext_sec.rd = FrRand();
  com_ext_sec.sd = FrRand();
  com_ext_sec.t1 = FrRand();
  com_ext_sec.t0 = FrRand();

  using details::ComputeCommitment;
  Fr xdy_dxy = InnerProduct(input.x(), com_ext_sec.dyt) +
               InnerProduct(com_ext_sec.dx, input.yt());

  //std::cout << Tick::GetIndentString() << "2*multiexp(" << input.n() << ")\n";

  std::vector<parallel::Task> tasks(3);
  tasks[0] = [&com_ext_pub,&com_ext_sec]() mutable {
    com_ext_pub.ad =
        ComputeCommitment(com_ext_sec.dx, com_ext_sec.rd);
  };
  tasks[1] = [&com_ext_pub,&com_ext_sec]() mutable {
    com_ext_pub.bd =
        ComputeCommitment(com_ext_sec.dy, com_ext_sec.sd);
  };
  tasks[2] = [&com_ext_pub,&com_ext_sec,&xdy_dxy]() mutable {
    com_ext_pub.c1 = ComputeCommitment(xdy_dxy, com_ext_sec.t1);
    com_ext_pub.c0 =
        ComputeCommitment(com_ext_sec.dz, com_ext_sec.t0);
  };
  parallel::Invoke(tasks);

////#ifdef MULTICORE
////#pragma omp parallel sections
////#endif
//  {
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_ext_pub.ad =
//        ComputeCommitment(com_ext_sec.dx, com_ext_sec.rd);
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_ext_pub.bd =
//        ComputeCommitment(com_ext_sec.dy, com_ext_sec.sd);
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_ext_pub.c1 = ComputeCommitment(xdy_dxy, com_ext_sec.t1);
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    com_ext_pub.c0 =
//        ComputeCommitment(com_ext_sec.dz, com_ext_sec.t0);
//  }
}

inline void ComputeProof(Proof& proof, ProverInput const& input,
                  CommitmentSec const& com_sec,
                  CommitmentExtSec const& com_ext_sec, Fr const& challenge) {
  //Tick tick(__FUNCTION__);
  auto n = input.n();
  proof.fx.resize(n);
  proof.fy.resize(n);
  //for (int64_t i = 0; i < n; ++i) {
  //  // fx = e * x + dx
  //  proof.fx[i] = challenge * input.x(i) + com_ext_sec.dx[i];
  //} 
  //for (int64_t i = 0; i < n; ++i) {
  //  // fy = e * y + dy
  //  proof.fy[i] = challenge * input.y(i) + com_ext_sec.dy[i];
  //}
  auto parallel_f = [&proof, &challenge, &input, &com_ext_sec](int64_t i) {
    // fx = e * x + dx
    proof.fx[i] = challenge * input.x(i) + com_ext_sec.dx[i];
    // fy = e * y + dy
    proof.fy[i] = challenge * input.y(i) + com_ext_sec.dy[i];
  };
  parallel::For(n, parallel_f, n < 16 * 1024);

  proof.rx = challenge * com_sec.r + com_ext_sec.rd;
  proof.sy = challenge * com_sec.s + com_ext_sec.sd;
  proof.tz = challenge * challenge * com_sec.t + challenge * com_ext_sec.t1 +
             com_ext_sec.t0;
}

inline void UpdateSeed(h256_t& seed, CommitmentPub const com_pub,
                       CommitmentExtPub const& com_ext_pub) {
  using details::HashUpdate;
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
  VerifierInput(std::vector<Fr> const* t, CommitmentPub const& com_pub)
      : t(t), com_pub(com_pub) {}
  std::vector<Fr> const* t;  // n, optional, can be null
  CommitmentPub const& com_pub;
};

inline bool VerifyInternal(VerifierInput const& input,
                    Fr const& challenge, CommitmentExtPub const& com_ext_pub,
                    Proof const& proof) {
  using details::ComputeCommitment;

  auto const n = proof.fx.size();
  assert(n == proof.fy.size());

  //std::cout << Tick::GetIndentString() << "multiexp(" << n << ")\n";

  bool ret1 = false;
  bool ret2 = false;
  auto const& com_pub = input.com_pub;

  std::vector<parallel::Task> tasks(2);
  tasks[0] = [&ret1, &com_pub, &challenge, &com_ext_pub, &proof]() {
    Fr alpha = FrRand();
    // (a^e * a_d)^alpha * b^e * b_d
    G1 left = (com_pub.a * challenge + com_ext_pub.ad) * alpha +
              (com_pub.b * challenge + com_ext_pub.bd);
    // com(alpha * fx + fy,alpha * rx + sy)
    std::vector<Fr> alpha_fx_fy;
    details::VectorMul(alpha_fx_fy, proof.fx, alpha);
    details::VectorAdd(alpha_fx_fy, alpha_fx_fy, proof.fy);
    Fr alpha_rx_sy = alpha * proof.rx + proof.sy;
    G1 right = ComputeCommitment(alpha_fx_fy, alpha_rx_sy);
    ret1 = left == right;
    assert(ret1);
  };

  tasks[1] = [&ret2, &com_pub, &challenge, &com_ext_pub, &proof, &input,
              n]() {
    // c^(e^2) * c_1^e * c_0 == com(f_x * f_y , t_z)
    Fr e2_square = challenge * challenge;
    G1 left =
        com_pub.c * e2_square + com_ext_pub.c1 * challenge + com_ext_pub.c0;
    Fr fz = [&input, &proof, n]() {
      if (input.t) {
        std::vector<Fr> proof_fyt(n);
        details::HadamardProduct(proof_fyt, proof.fy, *input.t);
        return InnerProduct(proof.fx, proof_fyt);
      } else {
        return InnerProduct(proof.fx, proof.fy);
      }
    }();

    G1 right = ComputeCommitment(fz, proof.tz);
    ret2 = left == right;
    assert(ret2);
  };

  parallel::Invoke(tasks);

////#ifdef MULTICORE
////#pragma omp parallel sections
////#endif
//  {
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    {
//      Fr alpha = FrRand();
//      // (a^e * a_d)^alpha * b^e * b_d
//      G1 left = (com_pub.a * challenge + com_ext_pub.ad) * alpha +
//                (com_pub.b * challenge + com_ext_pub.bd);
//      // com(alpha * fx + fy,alpha * rx + sy)
//      std::vector<Fr> alpha_fx_fy;
//      details::VectorMul(alpha_fx_fy, proof.fx, alpha);
//      details::VectorAdd(alpha_fx_fy, alpha_fx_fy, proof.fy);
//      Fr alpha_rx_sy = alpha * proof.rx + proof.sy;
//      G1 right = ComputeCommitment(alpha_fx_fy, alpha_rx_sy);
//      ret1 = left == right;
//      assert(ret1);
//    }
//
////#ifdef MULTICORE
////#pragma omp section
////#endif
//    {
//      // c^(e^2) * c_1^e * c_0 == com(f_x * f_y , t_z)
//      Fr e2_square = challenge * challenge;
//      G1 left =
//          com_pub.c * e2_square + com_ext_pub.c1 * challenge + com_ext_pub.c0;
//      Fr fz = [&input, &proof, n]() {
//        if (input.t) {
//          std::vector<Fr> proof_fyt(n);
//          details::HadamardProduct(proof_fyt, proof.fy, *input.t);
//          return InnerProduct(proof.fx, proof_fyt);
//        } else {
//          return InnerProduct(proof.fx, proof.fy);
//        }
//      }();
//
//      G1 right = ComputeCommitment(fz, proof.tz);
//      ret2 = left == right;
//      assert(ret2);
//    }
//  }
  assert(ret1 && ret2);
  return ret1 && ret2;
}

inline void RomProve(RomProof& rom_proof,
                     h256_t const& common_seed, ProverInput const& input,
                     CommitmentPub const& com_pub,
                     CommitmentSec const& com_sec) {
  //Tick tick(__FUNCTION__);

  assert(PdsPub::kGSize >= input.n());

  CommitmentExtSec com_ext_sec;
  ComputeComExt(rom_proof.com_ext_pub, com_ext_sec, input);

  auto seed = common_seed;
  UpdateSeed(seed, com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);
  
  ComputeProof(rom_proof.proof, input, com_sec, com_ext_sec, challenge);
}

inline bool RomVerify(RomProof const& rom_proof, 
               h256_t const& common_seed, VerifierInput const& input) {
  //Tick tick(__FUNCTION__);
  assert(PdsPub::kGSize >= rom_proof.n());

  auto seed = common_seed;
  UpdateSeed(seed, input.com_pub, rom_proof.com_ext_pub);
  Fr challenge = H256ToFr(seed);

  return VerifyInternal(input, challenge, rom_proof.com_ext_pub,
                        rom_proof.proof);
}

inline bool TestRom(int64_t n) {
  std::cout << "n=" << n << "\n";
  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> y(n);
  FrRand(y.data(), n);
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
}  // namespace groth09::sec51
