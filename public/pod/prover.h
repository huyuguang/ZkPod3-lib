#pragma once

#include "./details.h"
#include "./types.h"
#include "vrs/vrs.h"

namespace pod {

template<typename VrsScheme>
void EncryptAndProve(ProveOutput<VrsScheme>& output, h256_t seed,
                            CommitedData const& commited_data,
                            std::string cache_dir = "") {
  Tick _tick_(__FUNCTION__);
  auto n = commited_data.n;
  auto s = commited_data.s;

  if (cache_dir.empty()) {
    char const* data_dir = std::getenv("options:data_dir");
    cache_dir = data_dir ? data_dir : ".";
    cache_dir += "/vrs_cache/";
    cache_dir += VrsScheme::type();
  }

  output.auto_cache.reset(new vrs::AutoCacheFile<VrsScheme>(cache_dir, (n + 1) * (s + 1)));
  auto cache = output.auto_cache->LoadAndUpgrade();

  std::vector<std::vector<G1>> cached_var_coms;
  std::vector<std::vector<Fr>> cached_var_coms_r;
  if (cache) {
    output.proved_data.vrs_plain_seed = cache->seed;
    output.secret.seed0 = cache->key;
    output.secret.seed0_com_r = cache->key_com_r;
    cached_var_coms = std::move(cache->var_coms);
    cached_var_coms_r = std::move(cache->var_coms_r);
  } else {
    output.proved_data.vrs_plain_seed = misc::RandH256();
    output.secret.seed0 = FrRand();
    output.secret.seed0_com_r = FrRand();
  }

  // generate vrs plain by vrs_plain_seed
  std::vector<Fr> plain((n + 1) * (s + 1));
  auto parallel_f_p = [&plain, &output](uint64_t i) {
    plain[i] = vrs::GeneratePlain(output.proved_data.vrs_plain_seed, i);
  };
  parallel::For(plain.size(), parallel_f_p);

  // v = mimc_enc(plain, seed0)
  std::vector<Fr> v((n + 1) * (s + 1));
  auto parallel_f_v = [&v, &plain, &output](uint64_t i) {
    v[i] = vrs::Mimc5Enc(plain[i], output.secret.seed0);
  };
  parallel::For(v.size(), parallel_f_v);

  // k=com(v)
  output.proved_data.k.resize(n + 1);
  auto parallel_f_k = [&v, &output, s](uint64_t i) {
    Fr const* vi0 = &v[i * (s + 1)];
    output.proved_data.k[i] = MultiExpBdlo12<G1>(PcHG, vi0, s + 1);
    output.proved_data.k[i].normalize();
  };
  parallel::For(n + 1, parallel_f_k);

  // w
  details::UpdateSeed(seed, output.proved_data.k);
  std::vector<Fr> w(n + 1);
  ComputeFst(seed, "pod", w);

  // mij' = vij + mij * wi
  output.proved_data.em.resize(n * (s + 1));
  auto parallel_f_m = [&output, &commited_data, s, &v, &w](uint64_t i) {
    auto get_rm = [&commited_data](int64_t i, int64_t j)->Fr const& {
      assert(j < (commited_data.s + 1) && i < commited_data.n);
      if (!j) return commited_data.get_r(i);
      return commited_data.get_m(i, j - 1);
    };

    for (int64_t j = 0; j < s + 1; ++j) {
      auto ij = i * (s + 1) + j;
      output.proved_data.em[ij] = v[ij] + w[i] * get_rm(i, j);
    }
  };
  parallel::For(n, parallel_f_m);  

  // vw
  output.proved_data.vw.resize(s + 1);
  auto parallel_f_vw = [&output, n, s, &v, &w](int64_t j) {
    output.proved_data.vw[j] = FrZero();
    for (int64_t i = 0; i <= n; ++i) {
      output.proved_data.vw[j] += v[i * (s + 1) + j] * w[i];
    }
  };
  parallel::For(s + 1, parallel_f_vw);

  // use seed as fst_seed to prove vrs
  vrs::PublicInput public_input(plain.size(),
                                [&plain](int64_t i) { return plain[i]; });

  output.proved_data.vw_com_r = FrRand();
  vrs::SecretInput secret_input(output.secret.seed0, output.secret.seed0_com_r,
                                output.proved_data.vw_com_r);
  vrs::LargeProverLowRam<VrsScheme> prover(public_input, secret_input,
                                std::move(cached_var_coms),
                                std::move(cached_var_coms_r));
  auto get_w = [&w, s](int64_t i) { return w[i / (s + 1)]; };
  vrs::ProveOutput vrs_output;
  prover.Prove(seed, get_w, output.proved_data.vrs_proofs, vrs_output);

  // generate wanted receipt
  output.receipt.h = vrs_output.h;
  output.receipt.g = vrs_output.g;
  output.receipt.seed0_com = vrs_output.key_com;

#ifdef _DEBUG
  assert(v == prover.v());
  assert(parallel::Accumulate(output.proved_data.vw.begin(),
                              output.proved_data.vw.end(),
                              FrZero()) == prover.vw());
#endif
}
}  // namespace pod