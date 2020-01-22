#pragma once

#include "./pod.h"

namespace pod {
template<typename Scheme>
bool TestBasic(int64_t n, int64_t s) {
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
    int64_t g_offset = 0;
    com[i] = PcComputeCommitmentG(g_offset, s, m.data() + i * s, r[i]);
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
  assert(details::CheckCommitedData(commited_data));

  h256_t seed = misc::RandH256();

  Tick tick(__FUNCTION__);
  std::cout << "n = " << n << ", s = " << s << "\n";

  // prove
  ProveOutput<Scheme> prove_output;
  EncryptAndProve<Scheme>(prove_output, seed, commited_data);

  // prover send proveoutput.proved_data to verifier

  // verifier verify the proved_data, if ok, sign verify_output.receipt and send
  // to prover
  VerifyOutput verify_output;
  if (!VerifyAndSign<Scheme>(verify_output, seed, n, s, get_com,
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

  prove_output.auto_cache->SetLeaked();

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

inline bool Test() {
#ifdef _DEBUG
  int64_t n = 3;
  int64_t s = 5;
#else
  int64_t n = 32 * 24 - 1;
  int64_t s = 1023;
  //1023;
#endif

  //return TestBasic<vrs::Sha256cScheme>(n, s);
  return TestBasic<vrs::Mimc5Scheme>(n, s);
}
}  // namespace pod