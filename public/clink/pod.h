#pragma once

#include "./vrs/vrs_cache.h"
#include "./vrs/vrs_large.h"
#include "debug/flags.h"
#include "ecc/parallel_multiexp.h"

namespace clink {
template <typename Scheme, typename Policy>
struct Pod {
  using R1cs = typename clink::ParallelR1cs<Policy>;
  using HyraxA = typename Policy::HyraxA;
  using Vrs = VrsLarge<Scheme, Policy>;
  using VrsProveInput = typename Vrs::ProveInput;
  using VrsProof = typename Vrs::Proof;
  using VrsProveOutput = typename Vrs::ProveOutput;
  using VrsVerifyInput = typename Vrs::VerifyInput;
  using VrsVerifyOutput = typename Vrs::VerifyOutput;
  using AutoCacheFile = typename VrsCache<Scheme>::AutoFile;
  using AutoCacheFileUPtr = std::unique_ptr<AutoCacheFile>;

  typedef std::function<Fr const&(int64_t i, int64_t j)> GetMatrix;
  typedef std::function<Fr const&(int64_t i)> GetR;
  typedef std::function<G1 const&(int64_t i)> GetCom;

  static G1 const& kVwG() { return pc::PcU(); }

  // com is base on g_offset = 0
  struct CommitedData {
    int64_t n = 0;
    int64_t s = 0;
    GetMatrix get_m;
    GetR get_r;
    GetCom get_com;
  };

  struct ProvedData {
    std::vector<G1> k;   // size = n+1, k=com(v)
    std::vector<Fr> em;  // size = n*(s+1), encrypted r&m
    std::vector<Fr> vw;  // s + 1
    h256_t vrs_plain_seed;
    Fr vw_com_r;
    std::vector<VrsProof> vrs_proofs;
    bool operator==(ProvedData const& b) const {
      return k == b.k && em == b.em && vw == b.vw &&
             vrs_plain_seed == b.vrs_plain_seed && vw_com_r == b.vw_com_r &&
             vrs_proofs == b.vrs_proofs;
    }

    bool operator!=(ProvedData const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("pod.pd", ("k", k), ("em", em), ("vw", vw),
                         ("seed", vrs_plain_seed), ("r", vw_com_r),
                         ("p", vrs_proofs));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("pod.pd", ("k", k), ("em", em), ("vw", vw),
                         ("seed", vrs_plain_seed), ("r", vw_com_r),
                         ("p", vrs_proofs));
    }
  };

  struct Receipt {
    G1 h;
    G1 g;
    G1 key_com;
    bool operator==(Receipt const& b) const {
      return h == b.h && g == b.g && key_com == b.key_com;
    }
    bool operator!=(Receipt const& b) const { return !(*this == b); }
  };

  struct Secret {
    Fr key;
    Fr key_com_r;
  };

  struct ProveOutput {
    ProvedData proved_data;
    Receipt receipt;
    Secret secret;
    AutoCacheFileUPtr cache;
  };

  struct VerifyOutput {
    Receipt receipt;
    std::vector<Fr> plain;
    std::vector<Fr> w;
    Fr sigma_vw;
  };

  static void EncryptAndProve(ProveOutput& output, h256_t seed,
                              CommitedData const& commited_data,
                              std::string const& data_dir) {
    Tick _tick_(__FN__);
    auto n = commited_data.n;
    auto s = commited_data.s;

    std::string cache_dir;
    if (!data_dir.empty() && !debug::flags::disable_vrs_cache) {
      cache_dir = data_dir + "/vrs_cache/" + Scheme::type();
    }

    output.cache.reset(new AutoCacheFile(cache_dir, (n + 1) * (s + 1)));
    auto cache = output.cache->LoadAndUpgrade();

    std::vector<std::vector<G1>> cached_var_coms;
    std::vector<std::vector<Fr>> cached_var_coms_r;
    if (cache) {
      output.proved_data.vrs_plain_seed = cache->seed;
      output.secret.key = cache->key;
      output.secret.key_com_r = cache->key_com_r;
      cached_var_coms = std::move(cache->var_coms);
      cached_var_coms_r = std::move(cache->var_coms_r);
    } else {
      output.proved_data.vrs_plain_seed = misc::RandH256();
      output.secret.key = FrRand();
      output.secret.key_com_r = FrRand();
    }

    // generate vrs plain by vrs_plain_seed
    std::vector<Fr> plain((n + 1) * (s + 1));
    auto parallel_f_p = [&plain, &output](uint64_t i) {
      VrsPub<Scheme>::GeneratePlain(&plain[i],
                                    output.proved_data.vrs_plain_seed, i);
    };
    parallel::For(plain.size(), parallel_f_p);

    auto v = GenerateV(plain, output.secret.key);

    // k=com(v)
    output.proved_data.k.resize(n + 1);
    auto parallel_f_k = [&v, &output, s](uint64_t i) {
      Fr const* vi0 = &v[i * (s + 1)];
      output.proved_data.k[i] = MultiExpBdlo12<G1>(pc::PcHG, vi0, s + 1);
      output.proved_data.k[i].normalize();
    };
    parallel::For(n + 1, parallel_f_k);

    // w
    UpdateSeed(seed, output.proved_data.k);
    std::vector<Fr> w(n + 1);
    ComputeFst(seed, "pod", w);

    // mij' = vij + mij
    output.proved_data.em.resize(n * (s + 1));
    auto parallel_f_m = [&output, &commited_data, s, &v](uint64_t i) {
      auto get_rm = [&commited_data](int64_t i, int64_t j) -> Fr const& {
        assert(j < (commited_data.s + 1) && i < commited_data.n);
        if (!j) return commited_data.get_r(i);
        return commited_data.get_m(i, j - 1);
      };

      for (int64_t j = 0; j < s + 1; ++j) {
        auto ij = i * (s + 1) + j;
        output.proved_data.em[ij] = v[ij] + get_rm(i, j);
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

    // vrs
    VrsProve(s, seed, plain, w, output, std::move(cached_var_coms),
             std::move(cached_var_coms_r));
  }

  static bool VerifyAndSign(VerifyOutput& output, h256_t seed, int64_t n,
                            int64_t s, GetCom const& get_com,
                            ProvedData const& proved_data) {
    Tick _tick_(__FN__);
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
    UpdateSeed(seed, proved_data.k);
    output.w.resize(n + 1);
    ComputeFst(seed, "pod", output.w);

    std::array<parallel::VoidTask, 2> tasks;
    std::array<int64_t, 2> ret{0, 0};

    // check consistency of the encrypted m, vm and k.
    tasks[0] = [&ret, &get_com, &proved_data, &output, n, s]() {
      Tick tick("Check consistency of the encrypted m, vm and k.");
      G1 left1 = G1Zero();
      for (int64_t i = 0; i < n; ++i) {
        left1 += get_com(i);
      }
      auto const& k = proved_data.k;
      for (int64_t i = 0; i < n; ++i) {
        left1 += k[i];
      }

      // MultiExp(n*(s+1))! 70% of the time here.
      auto get_g = [s](int64_t ij) -> G1 const& {
        return pc::PcHG(ij % (s + 1));
      };
      G1 right1 =
          ParallelMultiExpBdlo12<G1>(get_g, proved_data.em, n * (s + 1));
      if (left1 != right1) return;

      G1 left2 = MultiExpBdlo12(proved_data.k, output.w);
      G1 right2 = MultiExpBdlo12<G1>(pc::PcHG, proved_data.vw, s + 1);
      if (left2 != right2) return;
      ret[0] = 1;
    };

    // check vrs
    VrsVerifyOutput vrs_output;
    tasks[1] = [&ret, n, s, &seed, &proved_data, &output, &vrs_output]() {
      ret[1] = VrsVerify(n, s, seed, proved_data, output, vrs_output);
    };

    parallel::Invoke(tasks);

    if (!ret[0] || !ret[1]) {
      assert(false);
      return false;
    }

    // generate receipt
    output.receipt.g = vrs_output.g;
    output.receipt.h = vrs_output.h;
    output.receipt.key_com = vrs_output.key_com;
    return true;
  }

  static std::vector<Fr> GenerateV(std::vector<Fr> const& plain,
                                   Fr const& key) {
    Tick _tick_(__FN__);
    // ex: v = mimc_enc(plain, key)
    std::vector<Fr> v(plain.size());
    auto parallel_f = [&v, &plain, &key](uint64_t i) {
      v[i] = Scheme::Generate(plain[i], key);
    };
    parallel::For(v.size(), parallel_f);
    return v;
  }

  static bool DecryptData(int64_t n, int64_t s,
                          std::vector<Fr> const& encrypted_m,
                          Secret const& secret, VerifyOutput const& output,
                          std::vector<Fr>& decrypted_m) {
    Tick _tick_(__FN__);
    // verify the secret opened by prover
    if (!VerifySecret(output.receipt, secret)) {
      assert(false);
      return false;
    }

    // compute v
    auto v = GenerateV(output.plain, secret.key);

#ifdef _DEBUG
    Fr check_sigma_vw = FrZero();
    for (size_t i = 0; i < v.size(); ++i) {
      check_sigma_vw += v[i] * output.w[i / (s + 1)];
    }
    assert(check_sigma_vw == output.sigma_vw);
#endif

    // decrypt m
    decrypted_m.resize(n * s);
    auto parallel_f2 = [&v, &decrypted_m, &encrypted_m, s](int64_t i) {
      for (int64_t j = 0; j < s; ++j) {
        auto d_ij = i * s + j;
        auto e_ij = i * (s + 1) + j + 1;
        decrypted_m[d_ij] = encrypted_m[e_ij] - v[e_ij];
      }
    };
    parallel::For(n, parallel_f2);
    return true;
  }

  static bool VerifySecret(Receipt const& receipt, Secret const& secret) {
    return VrsPub<Scheme>::VerifySecret(receipt.h, receipt.g, receipt.key_com,
                                        secret.key_com_r, secret.key);
  }

  static bool Test(int64_t n, int64_t s, std::string const& cache_dir);

 private:
  static void UpdateSeed(h256_t& seed, std::vector<G1> const& k) {
    // Tick tick(__FN__);
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, k);
    hash.Final(seed.data());
  }

  static bool CheckCommitedData(CommitedData const& data) {
    if (!data.n || !data.s || data.s > pc::Base::GSize()) return false;

    bool all_success = false;
    auto parallel_f = [&data](int64_t i) {
      auto get_x = [&data, i](int64_t j) -> Fr const& {
        return data.get_m(i, j);
      };
      auto const& r = data.get_r(i);
      return pc::ComputeCom(data.s, get_x, r) == data.get_com(i);
    };
    parallel::For(&all_success, data.n, parallel_f);
    return all_success;
  }

  static void VrsProve(int64_t s, h256_t const& seed,
                       std::vector<Fr> const& plain, std::vector<Fr> const& w,
                       ProveOutput& output,
                       std::vector<std::vector<G1>>&& cached_var_coms,
                       std::vector<std::vector<Fr>>&& cached_var_coms_r) {
    Tick tick(__FN__);
    auto get_p = [&plain](int64_t i) -> Fr const& { return plain[i]; };
    auto get_w = [&w, s](int64_t i) -> Fr const& { return w[i / (s + 1)]; };
    output.proved_data.vw_com_r = FrRand();
    VrsProveInput vrs_prove_input(plain.size(), output.secret.key,
                                  output.secret.key_com_r, std::move(get_p),
                                  std::move(get_w), output.proved_data.vw_com_r,
                                  kVwG());
    VrsProveOutput vrs_output;
    Vrs::Prove(output.proved_data.vrs_proofs, vrs_output, seed,
               std::move(vrs_prove_input), std::move(cached_var_coms),
               std::move(cached_var_coms_r));

    // generate wanted receipt
    output.receipt.h = vrs_output.h;
    output.receipt.g = vrs_output.g;
    output.receipt.key_com = vrs_output.key_com;
  }

  static bool VrsVerify(int64_t n, int64_t s, h256_t const& seed,
                        ProvedData const& proved_data, VerifyOutput& output,
                        VrsVerifyOutput& vrs_output) {
    output.plain.resize((n + 1) * (s + 1));
    auto parallel_f = [&output, &proved_data](int64_t i) {
      VrsPub<Scheme>::GeneratePlain(&output.plain[i],
                                    proved_data.vrs_plain_seed, i);
    };
    parallel::For((int64_t)output.plain.size(), parallel_f);

    auto get_p = [&output](int64_t i) -> Fr const& { return output.plain[i]; };
    auto get_w = [&output, s](int64_t i) -> Fr const& {
      return output.w[i / (s + 1)];
    };

    VrsVerifyInput vrs_input(output.plain.size(), std::move(get_p),
                             std::move(get_w), kVwG());

    if (!Vrs::Verify(vrs_output, proved_data.vrs_proofs, seed,
                     std::move(vrs_input))) {
      assert(false);
      return false;
    }

    auto const& vw = proved_data.vw;
    output.sigma_vw = parallel::Accumulate(vw.begin(), vw.end(), FrZero());
    G1 check_vw_com =
        pc::ComputeCom(kVwG(), output.sigma_vw, proved_data.vw_com_r);
    if (Vrs::SumProofsComVw(proved_data.vrs_proofs) != check_vw_com) {
      assert(false);
      return false;
    }
    return true;
  }
};

template <typename Scheme, typename Policy>
bool Pod<Scheme, Policy>::Test(int64_t n, int64_t s,
                               std::string const& data_dir) {
  std::vector<Fr> m(n * s);
  FrRand(m.data(), m.size());
  auto get_m = [s, &m](int64_t i, int64_t j) -> Fr const& {
    return m[i * s + j];
  };

  std::vector<Fr> r(n);
  FrRand(r.data(), r.size());
  auto get_r = [&r](int64_t i) -> Fr const& { return r[i]; };

  std::vector<G1> com(n);
  auto parallel_f = [&m, &com, &r, s](int64_t i) {
    com[i] = pc::ComputeCom(s, m.data() + i * s, r[i]);
    com[i].normalize();
  };
  parallel::For(n, parallel_f);
  auto get_com = [&com](int64_t i) -> G1 const& { return com[i]; };

  CommitedData commited_data;
  commited_data.n = n;
  commited_data.s = s;
  commited_data.get_m = get_m;

  commited_data.get_r = get_r;
  commited_data.get_com = get_com;
  assert(CheckCommitedData(commited_data));

  h256_t seed = misc::RandH256();

  Tick tick(__FN__);
  std::cout << "n = " << n << ", s = " << s << "\n";

  // prove
  ProveOutput prove_output;
  EncryptAndProve(prove_output, seed, commited_data, data_dir);

  // prover send prove_output.proved_data to verifier
#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(prove_output.proved_data);
  std::cout << "proof with encrypted data size: " << os.get_shared_buffer().size
            << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  ProvedData proved_data2;
  ia.serialize(proved_data2);
  if (prove_output.proved_data != proved_data2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  // verifier verify the proved_data, if ok, sign verify_output.receipt and send
  // to prover
  VerifyOutput verify_output;
  if (!VerifyAndSign(verify_output, seed, n, s, get_com,
                     prove_output.proved_data)) {
    assert(false);
    return false;
  }

  // prover verify the signed receipt
  if (prove_output.receipt != verify_output.receipt) {
    assert(false);
    return false;
  }

  // prover indicate the signed receipt and the secret to notary(blockchain)
  if (prove_output.cache) {
    prove_output.cache->SetLeaked();
  }

  // notary verify the consistent of the receipt and the secret, if ok,
  // automaticly transfer money
  if (!VerifySecret(prove_output.receipt, prove_output.secret)) {
    assert(false);
    return false;
  }

  // verifier decrypt the encrypted data by secret
  std::vector<Fr> decrypted_m;
  if (!DecryptData(n, s, prove_output.proved_data.em, prove_output.secret,
                   verify_output, decrypted_m)) {
    assert(false);
    return false;
  }

  if (m != decrypted_m) {
    assert(false);
    return false;
  }

  std::cout << "success\n";
  return true;
}
}  // namespace clink