#pragma once

#include "clink/clink.h"

// match query and pod

namespace cmd {

template <typename Policy>
struct MatchQuery {
  using Pod = typename clink::Pod<clink::VrsMimc5Scheme, Policy>;
  using MatchPack = typename clink::MatchPack<Policy>;

  struct Proof {
    typename Pod::ProvedData pod_proved_data;
    std::vector<typename MatchPack::Proof> mp_proofs;
    bool operator==(Proof const& b) const {
      return pod_proved_data == b.pod_proved_data && mp_proofs == b.mp_proofs;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("mq.p", ("pod", pod_proved_data), ("sp", mp_proofs));
    }

    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("mq.p", ("pod", pod_proved_data), ("sp", mp_proofs));
    }
  };

  struct ProveOutput {
    typename Pod::ProveOutput pod_output;
    std::vector<typename MatchPack::ProveOutput> mp_outputs;

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

  struct ProveInput {
    ProveInput(Fr const& key, typename Pod::CommitedData const& data_x,
               GetRefG1 const& get_gx, std::string const& data_dir)
        : key(key),
          data_x(data_x),
          get_gx(get_gx),
          data_dir(data_dir),
          n(data_x.n),
          s(data_x.s) {}
    Fr const& key;
    typename Pod::CommitedData const& data_x;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gpy = pc::kGetRefG1;  // must be 0 because of pod
    std::string data_dir;
    int64_t const n;
    int64_t const s;
  };

  static void ProveLine(ProveOutput& output, h256_t seed,
                        ProveInput const& input, int64_t i) {
    auto& sp_output = output.mp_outputs[i];

    std::vector<Fr> x(input.s);
    for (int64_t j = 0; j < input.s; ++j) {
      x[j] = input.data_x.get_m(i, j);
    }

    G1 com_x = input.data_x.get_com(i);
    Fr com_x_r = input.data_x.get_r(i);

    typename clink::MatchPack<Policy>::ProveInput m_input(
        input.key, x, com_x, com_x_r, input.get_gx, input.get_gpy);
    clink::MatchPack<Policy>::Prove(sp_output, seed, m_input);
  }

  static void Prove(ProveOutput& output, h256_t seed, ProveInput const& input) {
    Tick tick(__FN__);
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
    typename Pod::CommitedData data_y;
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

    Pod::EncryptAndProve(output.pod_output, seed, data_y, input.data_dir);
  }

  struct VerifyInput {
    VerifyInput(Fr const& key, int64_t s, std::vector<G1> const& com_x,
                GetRefG1 const& get_gx)
        : key(key),
          s(s),
          com_x(com_x),
          get_gx(get_gx),
          n((int64_t)com_x.size()) {}
    Fr const& key;
    int64_t const s;
    std::vector<G1> const& com_x;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gpy = pc::kGetRefG1;
    int64_t const n;
  };

  static bool Verify(Proof const& proof, h256_t seed, VerifyInput const& input,
                     typename Pod::VerifyOutput& output) {
    Tick tick(__FN__);
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
      typename clink::MatchPack<Policy>::VerifyInput m_input(
          s, input.key, input.get_gx, input.get_gpy);
      return clink::MatchPack<Policy>::Verify(sp_proof, seed, m_input);
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
    if (!Pod::VerifyAndSign(output, seed, pod_n, pod_s, get_com,
                            proof.pod_proved_data)) {
      assert(false);
      return false;
    }
    return true;
  }

  static bool DecryptData(int64_t n, int64_t s,
                          typename Pod::ProvedData const& proved_data,
                          typename Pod::Secret const& secret,
                          typename Pod::VerifyOutput const& verify_output,
                          std::vector<boost::dynamic_bitset<uint8_t>>& rets) {
    Tick tick(__FN__);
    int64_t pack_s = (s + 252) / 253;
    std::vector<Fr> m;
    if (!Pod::DecryptData(n, pack_s, proved_data.em, secret, verify_output,
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

  static bool Test(int64_t n, int64_t s, std::string const& key,
                   std::string const& data_dir);
};

template <typename Policy>
bool MatchQuery<Policy>::Test(int64_t n, int64_t s, std::string const& key,
                              std::string const& data_dir) {
  if (key.size() > 31) {
    std::cout << "invalid parameter: k.size() must <= 31.\n";
    return false;
  }
  if (s >= pc::Base::GSize() / 2) {
    std::cout << "invalid parameter: s must < " << pc::Base::GSize() / 2
              << "\n";
    return false;
  }

  auto seed = misc::RandH256();
  int64_t x_g_offset = 10;
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };

  Fr fr_key = PackStrToFr(key.c_str());
  std::vector<std::vector<Fr>> x(n);
  auto parallel_f = [&x, s](int64_t i) {
    auto& xi = x[i];
    xi.resize(s);
    FrRand(xi);
  };
  parallel::For(n, parallel_f);

  x[rand() % n][rand() % s] = fr_key;
  x[rand() % n][rand() % s] = fr_key;

  std::vector<boost::dynamic_bitset<uint8_t>> check_rets(n);
  for (auto i = 0; i < n; ++i) {
    check_rets[i].resize(s);
    for (auto j = 0; j < s; ++j) {
      check_rets[i][j] = x[i][j] == fr_key;
    }
  }

  std::vector<Fr> com_x_r(n);
  FrRand(com_x_r);

  std::vector<G1> com_x(n);
  auto parallel_f2 = [&x, &com_x, &com_x_r, &get_gx](int64_t i) {
    com_x[i] = pc::ComputeCom(get_gx, x[i], com_x_r[i]);
  };
  parallel::For(n, parallel_f2);

  Tick tick(__FN__);

  typename Pod::CommitedData data_x;
  data_x.n = n;
  data_x.s = s;
  data_x.get_com = [&com_x](int64_t i) -> G1 const& { return com_x[i]; };
  data_x.get_r = [&com_x_r](int64_t i) -> Fr const& { return com_x_r[i]; };
  data_x.get_m = [&x](int64_t i, int64_t j) -> Fr const& { return x[i][j]; };

  ProveInput prove_input(fr_key, data_x, get_gx, data_dir);

  ProveOutput prove_output;
  Prove(prove_output, seed, prove_input);
  Proof proof = prove_output.BuildProof();

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

  VerifyInput verify_input(fr_key, s, com_x, get_gx);
  typename Pod::VerifyOutput verify_output;
  if (!Verify(proof, seed, verify_input, verify_output)) {
    assert(false);
    return false;
  }

  // prover verify the signed receipt
  if (prove_output.pod_output.receipt != verify_output.receipt) {
    assert(false);
    return false;
  }

  prove_output.pod_output.cache->SetLeaked();

  std::vector<boost::dynamic_bitset<uint8_t>> rets;
  if (!DecryptData(n, s, proof.pod_proved_data, prove_output.pod_output.secret,
                   verify_output, rets)) {
    assert(false);
    return false;
  }

  bool success = check_rets == rets;
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace cmd