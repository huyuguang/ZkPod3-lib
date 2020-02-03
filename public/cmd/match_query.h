#pragma once

#include "pc_utils/pc_utils.h"
#include "pod/pod.h"

// match query and pod

namespace cmd {

template <typename Policy>
struct MatchQuery {
  struct Proof {
    pod::ProvedData<Policy> pod_proved_data;
    std::vector<typename pc_utils::MatchPack<Policy>::Proof> mp_proofs;
    bool operator==(Proof const& b) {
      return pod_proved_data == b.pod_proved_data && mp_proofs == b.mp_proofs;
    }

    bool operator!=(Proof const& b) { return !(*this == b); }
  };

  struct ProveOutput {
    pod::ProveOutput<vrs::Mimc5Scheme, Policy> pod_output;
    std::vector<typename pc_utils::MatchPack<Policy>::ProveOutput> mp_outputs;

    Proof BuildProof() const {
      Proof ret;
      ret.pod_proved_data = pod_output.proved_data;
      ret.mp_proofs.resize(mp_outputs.size());
      for (size_t i = 0; i < mp_outputs.size(); ++i) {
        ret.mp_proofs[i] = mp_outputs[i].BuildProof();
      }
      return ret;
    }
  };

  static void UpdateSeed(h256_t& seed, int64_t n, int64_t s, Fr const& key) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, n);
    HashUpdate(hash, s);
    HashUpdate(hash, key);
    hash.Final(seed.data());
  }

  struct ProverInput {
    ProverInput(Fr const& key, pod::CommitedData const& data_x,
                int64_t x_g_offset, std::string const& vrs_cache_dir = "")
        : key(key),
          data_x(data_x),
          x_g_offset(x_g_offset),
          vrs_cache_dir(vrs_cache_dir),
          n(data_x.n),
          s(data_x.s) {}
    Fr const& key;
    pod::CommitedData const& data_x;
    int64_t const x_g_offset;
    int64_t const py_g_offset = 0;  // must be 0 because of pod
    std::string vrs_cache_dir;
    int64_t const n;
    int64_t const s;
  };

  static void ProveLine(ProveOutput& output, h256_t seed,
                        ProverInput const& input, int64_t i) {
    auto& sp_output = output.mp_outputs[i];

    std::vector<Fr> x(input.s);
    for (int64_t j = 0; j < input.s; ++j) {
      x[j] = input.data_x.get_m(i, j);
    }

    G1 com_x = input.data_x.get_com(i);
    Fr com_x_r = input.data_x.get_r(i);

    pc_utils::MatchPack<Policy>::ProverInput m_input(
        input.key, x, com_x, com_x_r, input.x_g_offset, input.py_g_offset);
    pc_utils::MatchPack<Policy>::Prove(sp_output, seed, m_input);
  }

  static void Prove(ProveOutput& output, h256_t seed,
                    ProverInput const& input) {
    Tick tick(__FUNCTION__);
    int64_t n = input.n;
    int64_t s = input.s;

    UpdateSeed(seed, n, s, input.key);

    output.mp_outputs.resize(n);
    auto parallel_f = [&output, &input, &seed](int64_t i) {
      ProveLine(output, seed, input, i);
    };
    parallel::For(n, parallel_f);

    for (auto const& i : output.mp_outputs) {
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
      return output.mp_outputs[i].pack_y[j];
    };
    data_y.get_com = [&output](int64_t i) -> G1 const& {
      return output.mp_outputs[i].com_pack_y;
    };
    data_y.get_r = [&output](int64_t i) -> Fr const& {
      return output.mp_outputs[i].com_pack_y_r;
    };

    pod::EncryptAndProve<vrs::Mimc5Scheme>(output.pod_output, seed, data_y,
                                           input.vrs_cache_dir);
  }

  struct VerifierInput {
    VerifierInput(Fr const& key, int64_t s, std::vector<G1> const& com_x,
                  int64_t x_g_offset)
        : key(key),
          s(s),
          com_x(com_x),
          x_g_offset(x_g_offset),
          n((int64_t)com_x.size()) {}
    Fr const& key;
    int64_t const s;
    std::vector<G1> const& com_x;
    int64_t const x_g_offset;
    int64_t const py_g_offset = 0;
    int64_t const n;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifierInput const& input, pod::VerifyOutput& output) {
    Tick tick(__FUNCTION__);
    int64_t n = input.n;
    int64_t s = input.s;

    UpdateSeed(seed, n, s, input.key);

    if ((int64_t)proof.mp_proofs.size() != n) {
      assert(false);
      return false;
    }
    for (auto i = 0; i < n; ++i) {
      auto const& com_w = proof.mp_proofs[i].match_proof.com_w;
      if (com_w.empty() || com_w[0] != input.com_x[i]) {
        assert(false);
        return false;
      }
    }

    bool all_success = false;
    auto parallel_f = [&proof, &seed, s, &input](int64_t i) {
      auto const& sp_proof = proof.mp_proofs[i];
      pc_utils::MatchPack<Policy>::VerifierInput m_input(s, input.key, input.x_g_offset,
                                                 input.py_g_offset);
      return pc_utils::MatchPack<Policy>::Verify(sp_proof, seed, m_input);
    };
    parallel::For(&all_success, n, parallel_f);

    if (!all_success) {
      assert(false);
      return false;
    }

    auto get_com = [&proof](int64_t i) -> G1 const& {
      return proof.mp_proofs[i].com_pack_y;
    };
    auto pod_n = n;
    auto pod_s = (s + 252) / 253;
    if (!pod::VerifyAndSign<vrs::Mimc5Scheme>(output, seed, pod_n, pod_s,
                                              get_com, proof.pod_proved_data)) {
      assert(false);
      return false;
    }
    return true;
  }

  static bool DecryptData(int64_t n, int64_t s,
                          pod::ProvedData<Policy> const& proved_data,
                          pod::Secret const& secret,
                          pod::VerifyOutput const& verify_output,
                          std::vector<boost::dynamic_bitset<uint8_t>>& rets) {
    Tick tick(__FUNCTION__);
    int64_t pack_s = (s + 252) / 253;
    std::vector<Fr> m;
    if (!pod::DecryptData(n, pack_s, proved_data.em, secret, verify_output,
                          m)) {
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

  static bool Test();
};

// save to bin
template <typename Ar, typename Policy>
void serialize(Ar& ar, typename MatchQuery<Policy>::Proof const& t) {
  ar& YAS_OBJECT_NVP("mq.p", ("pod", t.pod_proved_data), ("sp", t.mp_proofs));
}

// load from bin
template <typename Ar, typename Policy>
void serialize(Ar& ar, typename MatchQuery<Policy>::Proof& t) {
  ar& YAS_OBJECT_NVP("mq.p", ("pod", t.pod_proved_data), ("sp", t.mp_proofs));
}

template <typename Policy>
bool MatchQuery<Policy>::Test() {
  auto seed = misc::RandH256();
  int64_t x_g_offset = 10;
  Fr key = FrRand();
  int64_t s = 100;
  int64_t n = 24;
  std::vector<std::vector<Fr>> x(n);
  auto parallel_f = [&x, s, &key](int64_t i) {
    auto& xi = x[i];
    xi.resize(s);
    FrRand(xi);
    xi[rand() % s] = key;
  };
  parallel::For(n, parallel_f);

  std::vector<boost::dynamic_bitset<uint8_t>> check_rets(n);
  for (auto i = 0; i < n; ++i) {
    check_rets[i].resize(s);
    for (auto j = 0; j < s; ++j) {
      check_rets[i][j] = x[i][j] == key;
    }
  }

  std::vector<Fr> com_x_r(n);
  FrRand(com_x_r);

  std::vector<G1> com_x(n);
  auto parallel_f2 = [&x, &com_x, &com_x_r, x_g_offset](int64_t i) {
    com_x[i] = PcComputeCommitmentG(x_g_offset, x[i], com_x_r[i]);
  };
  parallel::For(n, parallel_f2);

  Tick tick(__FUNCTION__);

  pod::CommitedData data_x;
  data_x.n = n;
  data_x.s = s;
  data_x.get_com = [&com_x](int64_t i) -> G1 const& { return com_x[i]; };
  data_x.get_r = [&com_x_r](int64_t i) -> Fr const& { return com_x_r[i]; };
  data_x.get_m = [&x](int64_t i, int64_t j) -> Fr const& { return x[i][j]; };

  ProverInput prover_input(key, data_x, x_g_offset);

  ProveOutput prove_output;
  Prove(prove_output, seed, prover_input);
  Proof proof = prove_output.BuildProof();

  VerifierInput verifier_input(key, s, com_x, x_g_offset);
  pod::VerifyOutput verify_output;
  if (!Verify(proof, seed, verifier_input, verify_output)) {
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

  bool success = check_rets == rets;
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n";
  return success;
}
}  // namespace cmd