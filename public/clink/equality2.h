#pragma once

#include "./details.h"
#include "hyrax/hyrax.h"

// batch proof of equality
//

namespace clink {
struct Equality2 {
  static void UpdateSeed(h256_t& seed, std::vector<G1> const& cv,
                         std::vector<G1> const& cu,
                         std::vector<size_t> const& v_size) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, cv);
    HashUpdate(hash, cu);
    HashUpdate(hash, v_size);
    hash.Final(seed.data());
  }

  static void Prove(hyrax::A3::Proof& proof, h256_t seed,
                    std::vector<std::vector<Fr>> const& v,
                    std::vector<G1> const& cv, std::vector<Fr> const& rv,
                    GetRefG1 const& get_gv, std::vector<G1> const& cu,
                    std::vector<Fr> const& ru, GetRefG1 const& get_gu) {
    CHECK(v.size() == cv.size() && v.size() == cu.size() &&
              v.size() == rv.size() && v.size() == ru.size(),
          "");

    auto m = v.size();
    size_t n = 0;
    std::vector<size_t> v_size(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
      v_size[i] = v[i].size();
      n = std::max<size_t>(n, v_size[i]);
    }

    UpdateSeed(seed, cv, cu, v_size);

    std::vector<Fr> r(m);
    std::vector<Fr> s(n);
    ComputeFst(seed, "equality2 r", r);
    ComputeFst(seed, "equality2 r", s);

    std::vector<Fr> V(n, FrZero());
    for (size_t i = 0; i < v.size(); ++i) {
      auto vi = v[i];
      vi.resize(n, FrZero());
      V += vi * r[i];
    }

    G1 cw = MultiExpBdlo12(cv, r) + MultiExpBdlo12(cu, r);
    Fr rw = InnerProduct(rv, r) + InnerProduct(ru, r);

    std::vector<Fr> S(n * 2);
    for (size_t i = 0; i < n; ++i) {
      S[i] = s[i];
      S[i + n] = -s[i];
    }

    std::vector<Fr> W(n * 2);
    for (size_t i = 0; i < n; ++i) {
      W[i] = V[i];
      W[i + n] = V[i];
    }

    auto get_g = [&get_gv, &get_gu, n](int64_t i) -> G1 const& {
      if (i < (int64_t)n)
        return get_gv(i);
      else
        return get_gu(i - n);
    };
    hyrax::A3::ProveInput input("equality2", W, S, FrZero(), get_g, G1Zero());
    hyrax::A3::CommitmentPub com_pub(cw, G1Zero());
    hyrax::A3::CommitmentSec com_sec(rw, FrZero());
    hyrax::A3::Prove(proof, seed, input, com_pub, com_sec);
  }

  static bool Verify(hyrax::A3::Proof const& proof, h256_t seed,
                     std::vector<size_t> const& v_size,
                     std::vector<G1> const& cv, GetRefG1 const& get_gv,
                     std::vector<G1> const& cu, GetRefG1 const& get_gu) {
    CHECK(cu.size() == cv.size() && v_size.size() == cv.size(), "");

    auto m = v_size.size();
    size_t n = 0;
    for (size_t i = 0; i < v_size.size(); ++i) {
      n = std::max<size_t>(n, v_size[i]);
    }

    UpdateSeed(seed, cv, cu, v_size);

    std::vector<Fr> r(m);
    std::vector<Fr> s(n);
    ComputeFst(seed, "equality2 r", r);
    ComputeFst(seed, "equality2 r", s);

    std::vector<Fr> S(n * 2);
    for (size_t i = 0; i < n; ++i) {
      S[i] = s[i];
      S[i + n] = -s[i];
    }

    G1 cw = MultiExpBdlo12(cv, r) + MultiExpBdlo12(cu, r);

    auto get_g = [&get_gv, &get_gu, n](int64_t i) -> G1 const& {
      if (i < (int64_t)n) {
        return get_gv(i);
      } else {
        return get_gu(i - n);
      }
    };

    hyrax::A3::CommitmentPub com_pub(cw, G1Zero());
    hyrax::A3::VerifyInput input("equality2", S, com_pub, get_g, G1Zero());
    return hyrax::A3::Verify(proof, seed, input);
  }

  static bool Test() {
    Tick tick(__FN__);
    auto seed = misc::RandH256();
    std::vector<std::vector<Fr>> x(10);
    for (size_t i = 0; i < x.size(); ++i) {
      x[i].resize(5 + i);
      FrRand(x[i]);
    }

    std::vector<G1> gv(x.size() + 5);
    G1Rand(gv);

    std::vector<G1> gu(x.size() + 5);
    G1Rand(gu);

    std::vector<G1> cv(x.size());
    std::vector<Fr> rv(x.size());
    FrRand(rv);

    auto get_gv = [&gv](int64_t i) -> G1 const& { return gv[i]; };
    for (size_t i = 0; i < x.size(); ++i) {
      cv[i] = pc::ComputeCom(get_gv, x[i], rv[i]);
    }

    std::vector<G1> cu(x.size());
    std::vector<Fr> ru(x.size());
    FrRand(ru);

    auto get_gu = [&gu](int64_t i) -> G1 const& { return gu[i]; };
    for (size_t i = 0; i < x.size(); ++i) {
      cu[i] = pc::ComputeCom(get_gu, x[i], ru[i]);
    }

    hyrax::A3::Proof proof;
    Prove(proof, seed, x, cv, rv, get_gv, cu, ru, get_gu);

    // size_t m = x.size();
    // size_t n = 0;
    std::vector<size_t> v_size(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      v_size[i] = x[i].size();
      // n = std::max<size_t>(n, v_size[i]);
    }
    bool success = Verify(proof, seed, v_size, cv, get_gv, cu, get_gu);
    std::cout << Tick::GetIndentString() << success << "\n\n\n\n\n\n";
    return success;
  }
};

}  // namespace clink