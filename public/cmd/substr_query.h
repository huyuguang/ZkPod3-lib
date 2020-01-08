#pragma once

#include "pc_utils/pc_utils.h"
#include "pod/pod.h"

// substr query and pod

namespace cmd::substr_query {

struct Proof {
  pod::ProvedData pod_proved_data;
  std::vector<pc_utils::substrpack::Proof> sp_proofs;
};

inline bool operator==(Proof const& a, Proof const& b) {
  return a.pod_proved_data == b.pod_proved_data && a.sp_proofs == b.sp_proofs;
}

inline bool operator!=(Proof const& a, Proof const& b) { return !(a == b); }

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof const& t) {
  ar& YAS_OBJECT_NVP("sq.p", ("pod", t.pod_proved_data), ("sp", t.sp_proofs));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof& t) {
  ar& YAS_OBJECT_NVP("sq.p", ("pod", t.pod_proved_data), ("sp", t.sp_proofs));
}

struct ProveOutput {
  pod::ProveOutput pod_output;
  std::vector<pc_utils::substrpack::ProveOutput> sp_outputs;

  Proof BuildProof() const {
    Proof ret;
    ret.pod_proved_data = pod_output.proved_data;
    ret.sp_proofs.resize(sp_outputs.size());
    for (auto i = 0; i < sp_outputs.size(); ++i) {
      ret.sp_proofs[i] = sp_outputs[i].BuildProof();
    }
    return ret;
  }
};

inline void UpdateSeed(h256_t& seed, int64_t n, int64_t s,
                       std::string const& key) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, n);
  HashUpdate(hash, s);
  HashUpdate(hash, key);
  hash.Final(seed.data());
}

inline void ProveLine(ProveOutput& output, h256_t seed, std::string const& k,
                      int64_t i, pod::CommitedData const& data_x) {
  auto& sp_output = output.sp_outputs[i];

  std::vector<Fr> x(data_x.s);
  for (int64_t j = 0; j < data_x.s; ++j) {
    x[j] = data_x.get_m(i, j);
  }

  G1 com_x = data_x.get_com(i);
  Fr com_x_r = data_x.get_r(i);

  pc_utils::substrpack::Prove(sp_output, seed, x, k, com_x, com_x_r);
}

inline void Prove(ProveOutput& output, h256_t seed, std::string const& key,
                  pod::CommitedData const& data_x,
                  std::string vrs_cache_dir = "") {
  Tick tick(__FUNCTION__);
  int64_t n = data_x.n;
  int64_t s = data_x.s;

  UpdateSeed(seed, n, s, key);

  output.sp_outputs.resize(n);
  auto parallel_f = [&output, &key, &data_x, &seed](int64_t i) {
    ProveLine(output, seed, key, i, data_x);
  };
  parallel::For(data_x.n, parallel_f);

  for (auto const& i : output.sp_outputs) {
    (void)i;
    assert((int64_t)i.y.size() == s);
    assert((int64_t)i.pack_y.size() == (s + 252LL) / 253LL);
  }

  // pod y to bob
  auto pod_n = n;
  auto pod_s = (s + 252) / 253;
  pod::CommitedData data_y;
  data_y.n = pod_n;
  data_y.s = pod_s;
  data_y.get_m = [&output](int64_t i, int64_t j) -> Fr const& {
    return output.sp_outputs[i].pack_y[j];
  };
  data_y.get_com = [&output](int64_t i) -> G1 const& {
    return output.sp_outputs[i].com_pack_y;
  };
  data_y.get_r = [&output](int64_t i) -> Fr const& {
    return output.sp_outputs[i].com_pack_y_r;
  };

  pod::EncryptAndProve(output.pod_output, seed, data_y, vrs_cache_dir);
}

inline bool Verify(Proof const& proof, h256_t seed, std::string const& key,
                   int64_t s, std::vector<G1> const& com_x,
                   pod::VerifyOutput& output) {
  Tick tick(__FUNCTION__);
  int64_t n = (int64_t)com_x.size();

  UpdateSeed(seed, n, s, key);

  if ((int64_t)proof.sp_proofs.size() != n) {
    assert(false);
    return false;
  }
  for (auto i = 0; i < n; ++i) {
    auto const& com_w = proof.sp_proofs[i].substr_proof.com_w;
    if (com_w.empty() || com_w[0] != com_x[i]) {
      assert(false);
      return false;
    }
  }

  bool all_success = false;
  auto parallel_f = [&proof, &seed, s, &key](int64_t i) {
    auto const& sp_proof = proof.sp_proofs[i];
    return pc_utils::substrpack::Verify(sp_proof, seed, s, key);
  };
  parallel::For(&all_success, n, parallel_f);
  if (!all_success) {
    assert(false);
    return false;
  }

  auto get_com = [&proof](int64_t i) -> G1 const& {
    return proof.sp_proofs[i].com_pack_y;
  };
  auto pod_n = n;
  auto pod_s = (s + 252) / 253;
  if (!pod::VerifyAndSign(output, seed, pod_n, pod_s, get_com,
                          proof.pod_proved_data)) {
    assert(false);
    return false;
  }
  return true;
}

inline bool DecryptData(int64_t n, int64_t s,
                        pod::ProvedData const& proved_data,
                        pod::Secret const& secret,
                        pod::VerifyOutput const& verify_output,
                        std::vector<boost::dynamic_bitset<uint8_t>>& rets) {
  Tick tick(__FUNCTION__);
  int64_t pack_s = (s + 252) / 253;
  std::vector<Fr> m;
  if (!pod::DecryptData(n, pack_s, proved_data.em, secret, verify_output, m)) {
    assert(false);
    return false;
  }

  assert((int64_t)m.size() == n * pack_s);

  rets.resize(n);

  auto parallel_f = [&rets, &m, pack_s, s](int64_t i) {
    rets[i] = FrsToBitset(m.data() + i * pack_s, pack_s);
    rets[i].resize(s);
  };
  parallel::For(n, parallel_f);

  return true;
}

inline bool Test() {
  auto seed = misc::RandH256();
  std::string key = "cde";
  int64_t s = 32 * 1000;
  int64_t n = 24;
  std::vector<std::vector<Fr>> x(n);
  auto parallel_f = [&x, s](int64_t i) {
    auto& xi = x[i];
    xi.resize(s);
    for (auto j = 0; j < s; ++j) {
      uint8_t buf[33];
      int len = rand() % 32;
      for (int k = 0; k < len; ++k) {
        buf[k] = (uint8_t)rand();
      }
      buf[len] = 0;
      xi[j] = PackStrToFr((char*)buf);
    }
  };
  parallel::For(n, parallel_f);

  x[0][0] = PackStrToFr("cdefg123");
  x[0][3] = PackStrToFr("abcdefg");
  x[1][9] = PackStrToFr("123abfgcde");

  std::vector<boost::dynamic_bitset<uint8_t>> check_rets(n);
  for (auto i = 0; i < n; ++i) {
    check_rets[i].resize(s);
    for (auto j = 0; j < s; ++j) {
      auto str = UnPackStrFromFr(x[i][j]);
      check_rets[i][j] = str.find(key) != std::string::npos;
    }
  }

  std::vector<Fr> com_x_r(n);
  FrRand(com_x_r);

  std::vector<G1> com_x(n);
  auto parallel_f2 = [&x, &com_x, &com_x_r](int64_t i) {
    com_x[i] = PcComputeCommitment(x[i], com_x_r[i]);
  };
  parallel::For(n, parallel_f);

  Tick tick(__FUNCTION__);

  pod::CommitedData data_x;
  data_x.n = n;
  data_x.s = s;
  data_x.get_com = [&com_x](int64_t i) -> G1 const& { return com_x[i]; };
  data_x.get_r = [&com_x_r](int64_t i) -> Fr const& { return com_x_r[i]; };
  data_x.get_m = [&x](int64_t i, int64_t j) -> Fr const& { return x[i][j]; };

  ProveOutput prove_output;
  Prove(prove_output, seed, key, data_x);
  Proof proof = prove_output.BuildProof();

  pod::VerifyOutput verify_output;
  if (!Verify(proof, seed, key, s, com_x, verify_output)) {
    assert(false);
    return false;
  }

  // prover verify the signed receipt
  if (prove_output.pod_output.receipt != verify_output.receipt) {
    assert(false);
    return false;
  }

  prove_output.pod_output.auto_cache->SetLeaked();

  std::vector<boost::dynamic_bitset<uint8_t>> rets;
  if (!DecryptData(n, s, proof.pod_proved_data, prove_output.pod_output.secret,
                   verify_output, rets)) {
    assert(false);
    return false;
  }

  if (check_rets != rets) {
    assert(false);
    std::cout << __LINE__ << ", oops\n";
  }

  std::cout << __FUNCTION__ << ": success\n";
  return true;
}
}  // namespace cmd::substr_query