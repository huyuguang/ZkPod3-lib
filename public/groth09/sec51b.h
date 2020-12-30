#pragma once

#include <numeric>

#include "groth09/details.h"
#include "utils/fst.h"

// t: public vector<Fr>, size = n
// X, Y: secret vector<Fr>, size = n
// z: secret Fr
// open: com(gx,X), com(gy,Y), com(gz,z)
// prove: z = <X, Y o t>
// proof size: (2n+3) Fr and 4 G1
// prove cost: 2*mulexp(n)
// verify cost: gx != gy? 2*mulexp(n) : mulexp(n)
namespace groth09 {
struct Sec51b {
  struct ProveInput {
    std::vector<Fr> const& x;   // size = n
    std::vector<Fr> const& y;   // size = n
    std::vector<Fr> const& t;   // size = n
    std::vector<Fr> const& yt;  // yt = y o t
    Fr const z;                 // z = <x, y o t>
    GetRefG1 const get_gx;
    GetRefG1 const get_gy;
    G1 const gz;
    bool const gx_equal_gy;

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
          gz(gz),
          gx_equal_gy(get_gx(0) == get_gy(0) &&
                      get_gx(x.size() - 1) == get_gy(y.size() - 1)) {
      assert(!x.empty());
      assert(x.size() == y.size());
      assert(x.size() == t.size());
      assert(x.size() == yt.size());
      assert(yt == HadamardProduct(y, t));
      assert(z == InnerProduct(x, yt));
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
    bool operator==(CommitmentExtPub const& right) const {
      return ad == right.ad && bd == right.bd && c1 == right.c1 &&
             c0 == right.c0;
    }

    bool operator!=(CommitmentExtPub const& right) const {
      return !(*this == right);
    }
  };

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

  struct SubProof {
    std::vector<Fr> fx;  // fx.size = n
    std::vector<Fr> fy;  // fy.size = n
    Fr rx;
    Fr sy;
    Fr tz;
    int64_t n() const { return (int64_t)fx.size(); }
    bool operator==(SubProof const& right) const {
      return fx == right.fx && fy == right.fy && rx == right.rx &&
             sy == right.sy && tz == right.tz;
    }

    bool operator!=(SubProof const& right) const { return !(*this == right); }
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;  // 4 G1
    SubProof sub_proof;            // (2n+3) Fr
    bool CheckFormat() const {
      return true;  // TODO:
    }
    int64_t n() const { return sub_proof.n(); }
    bool operator==(Proof const& right) const {
      return com_ext_pub == right.com_ext_pub && sub_proof == right.sub_proof;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }
  };

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

  static void ComputeComExt(CommitmentExtPub& com_ext_pub,
                            CommitmentExtSec& com_ext_sec,
                            ProveInput const& input) {
    // Tick tick(__FN__);
    auto const n = input.n();
    com_ext_sec.dx.resize(n);
    FrRand(com_ext_sec.dx.data(), n);
    com_ext_sec.dy.resize(n);
    FrRand(com_ext_sec.dy.data(), n);
    com_ext_sec.dyt.resize(n);
    HadamardProduct(com_ext_sec.dyt, com_ext_sec.dy, input.t);
    com_ext_sec.dz = InnerProduct(com_ext_sec.dx, com_ext_sec.dyt);

    com_ext_sec.rd = FrRand();
    com_ext_sec.sd = FrRand();
    com_ext_sec.t1 = FrRand();
    com_ext_sec.t0 = FrRand();

    Fr xdy_dxy = InnerProduct(input.x, com_ext_sec.dyt) +
                 InnerProduct(com_ext_sec.dx, input.yt);

    std::array<parallel::VoidTask, 3> tasks;
    tasks[0] = [&com_ext_pub, &com_ext_sec, &input]() {
      com_ext_pub.ad =
          pc::ComputeCom(input.get_gx, com_ext_sec.dx, com_ext_sec.rd);
    };
    tasks[1] = [&com_ext_pub, &com_ext_sec, &input]() {
      com_ext_pub.bd =
          pc::ComputeCom(input.get_gy, com_ext_sec.dy, com_ext_sec.sd);
    };
    tasks[2] = [&com_ext_pub, &com_ext_sec, &xdy_dxy, &input]() {
      com_ext_pub.c1 = pc::ComputeCom(input.gz, xdy_dxy, com_ext_sec.t1);
      com_ext_pub.c0 = pc::ComputeCom(input.gz, com_ext_sec.dz, com_ext_sec.t0);
    };
    parallel::Invoke(tasks);
  }

  static void ComputeSubProof(SubProof& sub_proof, ProveInput const& input,
                              CommitmentSec const& com_sec,
                              CommitmentExtSec const& com_ext_sec,
                              Fr const& challenge) {
    // Tick tick(__FN__);
    auto n = input.n();
    sub_proof.fx.resize(n);
    sub_proof.fy.resize(n);
    auto parallel_f = [&sub_proof, &challenge, &input,
                       &com_ext_sec](int64_t i) {
      // fx = e * x + dx
      sub_proof.fx[i] = challenge * input.x[i] + com_ext_sec.dx[i];
      // fy = e * y + dy
      sub_proof.fy[i] = challenge * input.y[i] + com_ext_sec.dy[i];
    };
    parallel::For(n, parallel_f, n < 16 * 1024);

    sub_proof.rx = challenge * com_sec.r + com_ext_sec.rd;
    sub_proof.sy = challenge * com_sec.s + com_ext_sec.sd;
    sub_proof.tz = challenge * challenge * com_sec.t +
                   challenge * com_ext_sec.t1 + com_ext_sec.t0;
  }

  static void UpdateSeed(h256_t& seed, CommitmentPub const com_pub,
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

  struct VerifyInput {
    VerifyInput(std::vector<Fr> const& t, CommitmentPub const& com_pub,
                GetRefG1 const& get_gx, GetRefG1 const& get_gy, G1 const& gz)
        : t(t),
          com_pub(com_pub),
          get_gx(get_gx),
          get_gy(get_gy),
          gz(gz),
          gx_equal_gy(get_gx(0) == get_gy(0)) {}
    std::vector<Fr> const& t;
    CommitmentPub const& com_pub;
    GetRefG1 const get_gx;
    GetRefG1 const get_gy;
    G1 const gz;
    bool const gx_equal_gy;
  };

  static bool VerifyInternal(VerifyInput const& input, Fr const& challenge,
                             CommitmentExtPub const& com_ext_pub,
                             SubProof const& sub_proof) {
    auto const n = sub_proof.fx.size();
    assert(n == sub_proof.fy.size());

    auto const& com_pub = input.com_pub;
    auto const& e = challenge;
    std::vector<int64_t> rets;
    std::vector<parallel::VoidTask> tasks;
    if (input.gx_equal_gy) {
      // if x_g_offset == y_g_offset, we can combine the verification
      tasks.resize(2);
      rets.resize(2);
      tasks[0] = [&rets, &com_pub, &e, &com_ext_pub, &sub_proof, &input]() {
        auto const& get_g = input.get_gx;
        Fr alpha = FrRand();
        // (a^e * a_d)^alpha * b^e * b_d
        G1 left = (com_pub.a * e + com_ext_pub.ad) * alpha +
                  (com_pub.b * e + com_ext_pub.bd);
        // com(alpha * fx + fy,alpha * rx + sy)
        std::vector<Fr> alpha_fx_fy = sub_proof.fx * alpha + sub_proof.fy;
        Fr alpha_rx_sy = alpha * sub_proof.rx + sub_proof.sy;
        G1 right = pc::ComputeCom(get_g, alpha_fx_fy, alpha_rx_sy);
        rets[0] = left == right;
        assert(rets[0]);
      };
    } else {
      tasks.resize(3);
      rets.resize(3);
      tasks[0] = [&rets, &com_pub, &e, &com_ext_pub, &sub_proof, &input]() {
        // a^e * a_d
        G1 left = (com_pub.a * e + com_ext_pub.ad);
        // com(fx,rx)
        G1 right = pc::ComputeCom(input.get_gx, sub_proof.fx, sub_proof.rx);
        rets[0] = left == right;
        assert(rets[0]);
      };
      tasks[1] = [&rets, &com_pub, &e, &com_ext_pub, &sub_proof, &input]() {
        // b^e * b_d
        G1 left = com_pub.b * e + com_ext_pub.bd;
        // com(fy,sy)
        G1 right = pc::ComputeCom(input.get_gy, sub_proof.fy, sub_proof.sy);
        rets[1] = left == right;
        assert(rets[1]);
      };
    }

    tasks.back() = [&rets, &com_pub, &e, &com_ext_pub, &sub_proof, &input,
                    n]() {
      // c^(e^2) * c_1^e * c_0 == com(f_x * f_y , t_z)
      Fr e2_square = e * e;
      G1 left = com_pub.c * e2_square + com_ext_pub.c1 * e + com_ext_pub.c0;

      std::vector<Fr> proof_fyt(n);
      HadamardProduct(proof_fyt, sub_proof.fy, input.t);
      Fr fz = InnerProduct(sub_proof.fx, proof_fyt);
      G1 right = pc::ComputeCom(input.gz, fz, sub_proof.tz);
      rets.back() = left == right;
      assert(rets.back());
    };

    parallel::Invoke(tasks);

    auto is_true = [](int64_t const& r) { return !!r; };
    auto all_success = std::all_of(rets.begin(), rets.end(), is_true);
    assert(all_success);
    return all_success;
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input,
                    CommitmentPub const& com_pub,
                    CommitmentSec const& com_sec) {
    // Tick tick(__FN__);
    CommitmentExtSec com_ext_sec;
    ComputeComExt(proof.com_ext_pub, com_ext_sec, input);

    UpdateSeed(seed, com_pub, proof.com_ext_pub);
    Fr challenge = H256ToFr(seed);

    ComputeSubProof(proof.sub_proof, input, com_sec, com_ext_sec, challenge);
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    // Tick tick(__FN__);
    UpdateSeed(seed, input.com_pub, proof.com_ext_pub);
    Fr challenge = H256ToFr(seed);

    return VerifyInternal(input, challenge, proof.com_ext_pub, proof.sub_proof);
  }

  static bool Test(int64_t n);
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::SubProof const& t) {
  ar& YAS_OBJECT_NVP("51.sp", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::SubProof& t) {
  ar& YAS_OBJECT_NVP("51.sp", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::Proof const& t) {
  ar& YAS_OBJECT_NVP("51.pf", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Sec51b::Proof& t) {
  ar& YAS_OBJECT_NVP("51.pf", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

bool Sec51b::Test(int64_t n) {
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
  std::vector<Fr> yt(n);

  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  GetRefG1 get_gy = [y_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[y_g_offset + i];
  };

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
