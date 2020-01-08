#pragma once

#include "./details.h"
#include "./types.h"
#include "./notary.h"

namespace pod {
inline bool VerifyAndSign(VerifyOutput& output, h256_t seed, int64_t n,
                          int64_t s, GetCom const& get_com,
                          ProvedData const& proved_data) {
  Tick _tick_(__FUNCTION__);
  if ((int64_t)proved_data.k.size() != (n + 1)) {
    assert(false);
    return false;
  }
  if ((int64_t)proved_data.em.size() != n * (s + 1)) {
    assert(false);
    return false;
  }
  if ((int64_t)proved_data.vw.size() != (s + 1)) {
    assert(false);
    return false;
  }

  // w
  details::UpdateSeed(seed, proved_data.k);
  output.w.resize(n + 1);
  ComputeFst(seed, "pod", output.w);
  
  std::array<parallel::Task, 2> tasks;
  std::array<int64_t, 2> ret{0, 0};

  // check consistency of the encrypted m, vm and k.
  tasks[0] = [&ret, &get_com, &proved_data, &output, n, s]() {
    //Tick tick(__FUNCTION__);
    auto const& k = proved_data.k;
    G1 left1 = MultiExpBdlo12<G1>(get_com, output.w, n);
    left1 = parallel::Accumulate(k.begin(), k.begin() + n, left1);
    auto get_g = [s](int64_t ij) -> G1 const& {
      return PcU(ij % (s + 1));
    };
    // MultiExp(n*(s+1))! 70% of the time here.
    G1 right1 = MultiExpBdlo12<G1>(get_g, proved_data.em, n * (s + 1));
    if (left1 != right1) return;

    G1 left2 = MultiExpBdlo12(proved_data.k, output.w);
    G1 right2 = MultiExpBdlo12<G1>(PcU, proved_data.vw, s + 1);
    if (left2 != right2) return;
    ret[0] = 1;
  };

  // check vrs
  vrs::VerifyOutput vrs_output;
  tasks[1] = [&ret, n, s, &seed, &proved_data, &output, &vrs_output]() {
    ret[1] = details::CheckVrs(n, s, seed, proved_data, output, vrs_output);
  };

  parallel::Invoke(tasks);

  if (!ret[0] || !ret[1]) {
    assert(false);
    return false;
  }

  // generate receipt
  output.receipt.g = vrs_output.g;
  output.receipt.h = vrs_output.h;
  output.receipt.seed0_com = vrs_output.key_com;
  return true;
}

inline bool DecryptData(int64_t n, int64_t s,
                        std::vector<Fr> const& encrypted_m,
                        Secret const& secret, VerifyOutput const& output,
                        std::vector<Fr>& decrypted_m) {
  Tick _tick_(__FUNCTION__);
  // verify the secret opened by prover
  if (!VerifySecret(output.receipt, secret)) {
    assert(false);
    return false;
  }

  // compute v
  std::vector<Fr> v(output.plain.size());
  auto parallel_f = [&secret, &output, &v](int64_t i) {
    v[i] = vrs::Mimc5Enc(output.plain[i], secret.seed0);
  };
  parallel::For((int64_t)output.plain.size(), parallel_f);

  #ifdef _DEBUG
    Fr check_sigma_vw = FrZero();
    for (size_t i = 0; i < v.size(); ++i) {
      check_sigma_vw += v[i] * output.w[i / (s + 1)];
    }
    assert(check_sigma_vw == output.sigma_vw);
  #endif

  // decrypt m
  std::vector<Fr> inv_w = output.w;
  FrInv(inv_w.data(), inv_w.size());
  decrypted_m.resize(n * s);
  auto parallel_f2 = [&v, &inv_w, &decrypted_m, &encrypted_m, s](int64_t i) {
    for (int64_t j = 0; j < s; ++j) {
      auto d_ij = i * s + j;
      auto e_ij = i * (s + 1) + j + 1;
      decrypted_m[d_ij] = (encrypted_m[e_ij] - v[e_ij]) * inv_w[i];
    }
  };
  parallel::For(n, parallel_f2);
  return true;
}
}  // namespace pod