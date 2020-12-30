#pragma once

#include "./details.h"

// recursive version of a2
// a: public vector<Fr>, size = n
// x: secret vector<Fr>, size = n
// y: secret Fr, = <x,a>
// open: com(gx,x), com(gy,y)
// prove: y=<x,a>
// proof size: (2+2log(n))G1 + 2Fr
namespace hyrax {
struct A3 {
  struct ProveInput {
    ProveInput(std::string const& tag, std::vector<Fr> const& x,
               std::vector<Fr> const& a, Fr const& y, GetRefG1 const& get_gx,
               G1 const& gy)
        : tag(tag), x(x), a(a), y(y), get_gx(get_gx), gy(gy) {
      assert(x.size() == a.size() && !a.empty());
      assert(y == InnerProduct(x, a));
    }
    int64_t n() const { return (int64_t)x.size(); }
    std::string to_string() const { return tag + ": " + std::to_string(n()); }

    std::string tag;
    std::vector<Fr> const& x;  // x.size = n // TODO: change to move
    std::vector<Fr> const& a;  // a.size = n
    Fr const y;                // y = <x, a>
    GetRefG1 const get_gx;
    G1 const gy;
  };

  struct CommitmentPub {
    CommitmentPub() {}
    CommitmentPub(G1 const& xi, G1 const& tau) : xi(xi), tau(tau) {}
    G1 xi;   // com(x,r_xi)
    G1 tau;  // com(y,r_tau)
    bool operator==(CommitmentPub const& right) const {
      return xi == right.xi && tau == right.tau;
    }

    bool operator!=(CommitmentPub const& right) const {
      return !(*this == right);
    }
  };

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
    bool operator==(CommitmentExtPub const& right) const {
      return delta == right.delta && beta == right.beta &&
             gamma_neg_1 == right.gamma_neg_1 &&
             gamma_pos_1 == right.gamma_pos_1;
    }

    bool operator!=(CommitmentExtPub const& right) const {
      return !(*this == right);
    }
  };

  struct CommitmentExtSec {
    // recursive rounds
    std::vector<Fr> r_gamma_neg_1;  // size=log(n)
    std::vector<Fr> r_gamma_pos_1;
    // finnal round
    Fr d;
    Fr r_beta;
    Fr r_delta;
  };

  struct SubProof {
    Fr z1;
    Fr z2;
    bool operator==(SubProof const& right) const {
      return z1 == right.z1 && z2 == right.z2;
    }

    bool operator!=(SubProof const& right) const { return !(*this == right); }
  };

  struct Proof {
    CommitmentExtPub com_ext_pub;
    SubProof sub_proof;
    int64_t aligned_n() const { return 1LL << com_ext_pub.gamma_neg_1.size(); }
    bool operator==(Proof const& right) const {
      return com_ext_pub == right.com_ext_pub && sub_proof == right.sub_proof;
    }

    bool operator!=(Proof const& right) const { return !(*this == right); }
    bool CheckFormat() const { return true; }
  };

  struct VerifyInput {
    VerifyInput(std::string const& tag, std::vector<Fr> const& a,
                CommitmentPub const& com_pub, GetRefG1 const& get_gx,
                G1 const& gy)
        : tag(tag), a(a), com_pub(com_pub), get_gx(get_gx), gy(gy) {}
    std::string tag;
    std::vector<Fr> const& a;  // a.size = n
    CommitmentPub const& com_pub;
    GetRefG1 const get_gx;
    G1 const gy;
    int64_t n() const { return (int64_t)a.size(); }
    std::string to_string() const { return tag + ": " + std::to_string(n()); }
  };

  static void UpdateSeed(h256_t& seed, std::vector<Fr> const& a,
                         CommitmentPub const& com_pub) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, a);
    HashUpdate(hash, com_pub.xi);
    HashUpdate(hash, com_pub.tau);
    hash.Final(seed.data());
  }

  static void UpdateSeed(h256_t& seed, G1 const& a, G1 const& b) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, a);
    HashUpdate(hash, b);
    hash.Final(seed.data());
  }

  // NOTE: maybe can optimize?
  static std::vector<G1> FuncO(std::vector<G1> const& g1, Fr const& k,
                               std::vector<G1> const& g2, Fr const& l) {
    assert(g1.size() == g2.size());
    std::vector<G1> g(g1.size());

    auto f = [&g, &g1, &g2, &k, &l](int64_t i) {
      g[i] = g1[i] * k + g2[i] * l;
    };
    parallel::For(g.size(), f, g.size() < 10240);

    return g;
  }

  template <typename T>
  static void Divide(std::vector<T> const& t, std::vector<T>& t1,
                     std::vector<T>& t2, T const& t0) {
    auto n = t.size();
    auto half = misc::Pow2UB(n) / 2;
    t1.resize(half);
    t2.resize(half);
    std::copy(t.begin(), t.begin() + half, t1.begin());
    std::copy(t.begin() + half, t.end(), t2.begin());
    std::fill(t2.begin() + (n - half), t2.end(), t0);
  }

  static void ComputeCom(CommitmentPub& com_pub, CommitmentSec const& com_sec,
                         ProveInput const& input) {
    Tick tick(__FN__);
    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&com_pub, &input, &com_sec]() {
      com_pub.xi = pc::ComputeCom(input.get_gx, input.x, com_sec.r_xi);
    };
    tasks[1] = [&com_pub, &input, &com_sec]() {
      com_pub.tau = pc::ComputeCom(input.gy, input.y, com_sec.r_tau);
    };
    parallel::Invoke(tasks, true);
  }

  // TODO: change to move
  static void Prove(Proof& proof, h256_t seed, ProveInput input,
                    CommitmentPub com_pub, CommitmentSec com_sec) {
    Tick tick(__FN__, input.to_string());
    UpdateSeed(seed, input.a, com_pub);

    auto x = input.x;
    auto a = input.a;
    auto y = input.y;
    int64_t n = x.size();
    int64_t round = (int64_t)misc::Log2UB(n);
    auto gx = pc::CopyG(input.get_gx, n);
    gx.resize(misc::Pow2UB(n), G1Zero());
    // std::fill(gx.begin() + n, gx.end(), G1Zero());
    auto const& h = pc::PcH();
    proof.com_ext_pub.gamma_neg_1.resize(round);
    proof.com_ext_pub.gamma_pos_1.resize(round);
    CommitmentExtSec com_ext_sec;
    com_ext_sec.r_gamma_neg_1.resize(round);
    com_ext_sec.r_gamma_pos_1.resize(round);
    auto r_gamma = com_sec.r_tau + com_sec.r_xi;
    CommitmentExtPub& com_ext_pub = proof.com_ext_pub;
    auto const& gy = input.gy;
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
      Fr x1_a2;
      Fr x2_a1;

      std::array<parallel::VoidTask, 2> tasks;
      tasks[0] = [&x1, &a2, &h, &r_gamma_neg_1, &gy, &g2, &gamma_neg_1,
                  &x1_a2]() {
        x1_a2 = InnerProduct(x1, a2);
        gamma_neg_1 = h * r_gamma_neg_1 + gy * x1_a2;
        gamma_neg_1 += MultiExpBdlo12(g2, x1);
      };
      tasks[1] = [&x2, &a1, &h, &r_gamma_pos_1, &gy, &g1, &gamma_pos_1,
                  &x2_a1]() {
        x2_a1 = InnerProduct(x2, a1);
        gamma_pos_1 = h * r_gamma_pos_1 + gy * x2_a1;
        gamma_pos_1 += MultiExpBdlo12(g1, x2);
      };
      parallel::Invoke(tasks, g2.size() < 10240);

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
    // std::cout << "gx[0]: " << gx[0] << "\n";
    // std::cout << "a[0]: " << a[0] << "\n";
    // std::cout << "gamma: " << gamma << "\n";

    // final round
    com_ext_sec.d = FrRand();
    com_ext_sec.r_beta = FrRand();
    com_ext_sec.r_delta = FrRand();
    com_ext_pub.delta = gx[0] * com_ext_sec.d + h * com_ext_sec.r_delta;
    com_ext_pub.beta = gy * com_ext_sec.d + h * com_ext_sec.r_beta;
    UpdateSeed(seed, com_ext_pub.delta, com_ext_pub.beta);
    Fr c = H256ToFr(seed);
    // std::cout << c << "\n";
    proof.sub_proof.z1 = com_ext_sec.d + c * y;
    proof.sub_proof.z2 =
        a[0] * (c * r_gamma + com_ext_sec.r_beta) + com_ext_sec.r_delta;
  }

  static void BuildS(std::vector<Fr>& s, std::vector<Fr> const& c,
                     std::vector<Fr> const& d, std::vector<Fr> const& cc) {
    Tick tick(__FN__);
    assert(c.size() == d.size());
    auto round = c.size();
    assert(s.size() <= (1ULL << round));

    for (size_t k = 0; k < s.size(); ++k) {
      if (k == 0) {
        s[k] = FrOne();
        for (size_t i = 0; i < round; ++i) {
          s[k] *= d[i];
        }
      } else {
        std::bitset<64> bits(k);
        for (size_t i = 0; i < round; ++i) {
          if (bits[round - i - 1]) {
            bits[round - i - 1] = 0;
            auto kk = bits.to_ulong();
            s[k] = s[kk] * cc[i];
            break;
          }
        }
      }

      //{
      //  Fr check = FrOne();
      //  std::bitset<64> bits(k);
      //  for (size_t i = 0; i < round; ++i) {
      //    if (bits[round - i - 1]) {
      //      check *= c[i];
      //    } else {
      //      check *= d[i];
      //    }
      //  }
      //  assert(s[k] == check);
      //}
    }
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FN__, input.to_string());
    auto n = input.n();
    if (!n || (int64_t)misc::Pow2UB(n) != proof.aligned_n()) {
      return false;
    }
    CommitmentPub const& com_pub = input.com_pub;
    CommitmentExtPub const& com_ext_pub = proof.com_ext_pub;
    int64_t round = (int64_t)misc::Log2UB(n);

    UpdateSeed(seed, input.a, com_pub);

    std::vector<Fr> vec_c(round);
    std::vector<Fr> vec_d(round);
    std::vector<Fr> vec_cc(round);
    std::vector<Fr> vec_dd(round);

    // recursive round
    for (int64_t loop = 0; loop < round; ++loop) {
      auto const& gamma_neg_1 = com_ext_pub.gamma_neg_1[loop];
      auto const& gamma_pos_1 = com_ext_pub.gamma_pos_1[loop];
      UpdateSeed(seed, gamma_neg_1, gamma_pos_1);
      vec_c[loop] = H256ToFr(seed);
      vec_cc[loop] = vec_c[loop] * vec_c[loop];
    }
    vec_d = vec_c;
    FrInv(vec_d);
    vec_dd = vec_cc;
    FrInv(vec_dd);

    std::vector<Fr> s(n);
    BuildS(s, vec_c, vec_d, vec_cc);

    G1 gx = MultiExpBdlo12<G1>(input.get_gx, s, s.size());
    Fr a = InnerProduct(input.a, s);

    G1 gamma = com_pub.xi + com_pub.tau;
    for (int64_t loop = 0; loop < round; ++loop) {
      auto const& gamma_neg_1 = com_ext_pub.gamma_neg_1[loop];
      auto const& gamma_pos_1 = com_ext_pub.gamma_pos_1[loop];
      gamma += gamma_neg_1 * vec_cc[loop] + gamma_pos_1 * vec_dd[loop];
    }

    // final round
    UpdateSeed(seed, com_ext_pub.delta, com_ext_pub.beta);
    Fr c = H256ToFr(seed);
    auto const& h = pc::PcH();
    auto const& gy = input.gy;
    auto const& sub_proof = proof.sub_proof;
    G1 left = (gamma * c + com_ext_pub.beta) * a + com_ext_pub.delta;
    G1 right = (gx + gy * a) * sub_proof.z1 + h * sub_proof.z2;
    if (left != right) {
      assert(false);
      return false;
    }
    return true;
  }

  static bool Test(int64_t n);
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, A3::CommitmentPub const& t) {
  ar& YAS_OBJECT_NVP("a3.cp", ("xi", t.xi), ("tau", t.tau));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A3::CommitmentPub& t) {
  ar& YAS_OBJECT_NVP("a3.cp", ("xi", t.xi), ("tau", t.tau));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A3::CommitmentExtPub const& t) {
  ar& YAS_OBJECT_NVP("a3.cep", ("gn1", t.gamma_neg_1), ("gp1", t.gamma_pos_1),
                     ("delta", t.delta), ("beta", t.beta));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A3::CommitmentExtPub& t) {
  ar& YAS_OBJECT_NVP("a3.cep", ("gn1", t.gamma_neg_1), ("gp1", t.gamma_pos_1),
                     ("delta", t.delta), ("beta", t.beta));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A3::SubProof const& t) {
  ar& YAS_OBJECT_NVP("a3.sp", ("z1", t.z1), ("z2", t.z2));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A3::SubProof& t) {
  ar& YAS_OBJECT_NVP("a3.sp", ("z1", t.z1), ("z2", t.z2));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, A3::Proof const& t) {
  ar& YAS_OBJECT_NVP("a3.rp", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, A3::Proof& t) {
  ar& YAS_OBJECT_NVP("a3.rp", ("c", t.com_ext_pub), ("p", t.sub_proof));
}

bool A3::Test(int64_t n) {
  Tick tick(__FN__);
  std::cout << "n = " << n << "\n";

  std::vector<Fr> x(n);
  FrRand(x.data(), n);
  std::vector<Fr> a(n);
  FrRand(a.data(), n);

  h256_t seed = misc::RandH256();

  int64_t x_g_offset = 20;
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  G1 gy = pc::PcU();
  auto y = InnerProduct(x, a);
  ProveInput prove_input("test", x, a, y, get_gx, pc::PcU());

  CommitmentPub com_pub;
  CommitmentSec com_sec(FrRand(), FrRand());
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

  VerifyInput verify_input("test", a, com_pub, get_gx, gy);
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace hyrax
