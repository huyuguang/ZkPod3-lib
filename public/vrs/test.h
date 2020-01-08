#include "./vrs.h"

namespace vrs {

inline void Test() {
  auto rom_seed = misc::RandH256();
  // int64_t offset = 0;
#ifdef _DEBUG
  int64_t count = 10;
#else
  int64_t count = 1024 * 32;
#endif

  h256_t vrs_plain_seed = misc::RandH256();
  auto get_p = [&vrs_plain_seed](int64_t i) {
    return GeneratePlain(vrs_plain_seed, i);
  };
  vrs::PublicInput public_input(count, std::move(get_p));

  vrs::SecretInput secret_input{FrRand(), FrRand(), FrRand()};

  vrs::Proof proof;
  vrs::ProveOutput prove_output;
  vrs::Prover prover(public_input, secret_input, std::vector<G1>(),
                     std::vector<Fr>());
  prover.Evaluate();

  std::vector<Fr> w(count);
  FrRand(w.data(), count);
  auto get_w = [&w](int64_t i) { return w[i]; };
  prover.Prove(rom_seed, get_w, proof, prove_output);

  auto com_vw =
      PcComputeCommitment(prover.vw(), secret_input.vw_com_r);
  assert(com_vw == proof.com_vw);

  vrs::VerifyOutput verify_output;
  vrs::Verifier verifier(public_input);
  bool ret = verifier.Verify(rom_seed, get_w, proof, verify_output);
  assert(ret);

  assert(prove_output.g == verify_output.g);
  assert(prove_output.h == verify_output.h);
  assert(prove_output.key_com == verify_output.key_com);

  ret = vrs::VerifySecret(prove_output.h, prove_output.g, prove_output.key_com,
                          secret_input.key_com_r, secret_input.key);
  assert(ret);
  std::cout << (ret ? "success" : "failed") << "\n";
}

inline void TestLarge() {
  auto rom_seed = misc::RandH256();
#ifdef _DEBUG
  int64_t count = 64;
#else
  int64_t count = 1024 * 128;
#endif
  h256_t vrs_plain_seed = misc::RandH256();
  auto get_p = [&vrs_plain_seed](int64_t i) {
    return GeneratePlain(vrs_plain_seed, i);
  };
  vrs::PublicInput public_input(count, std::move(get_p));

  vrs::SecretInput secret_input{FrRand(), FrRand(), FrRand()};

  std::vector<vrs::Proof> proofs;
  vrs::ProveOutput prove_output;
  vrs::LargeProver prover(public_input, secret_input,
                          std::vector<std::vector<G1>>(),
                          std::vector<std::vector<Fr>>());
  prover.Evaluate();

  std::vector<Fr> w(count);
  FrRand(w.data(), count);
  auto get_w = [&w](int64_t i) { return w[i]; };
  prover.Prove(rom_seed, get_w, proofs, prove_output);

  auto check_com_vw =
      PcComputeCommitment(prover.vw(), secret_input.vw_com_r);
  auto com_vw =
      std::accumulate(proofs.begin(), proofs.end(), G1Zero(),
                      [](G1 const& a, Proof const& b) { return a + b.com_vw; });
  assert(check_com_vw == com_vw);

  vrs::VerifyOutput verify_output;
  vrs::LargeVerifier verifier(public_input);
  bool ret = verifier.Verify(rom_seed, get_w, proofs, verify_output);
  assert(ret);

  assert(prove_output.g == verify_output.g);
  assert(prove_output.h == verify_output.h);
  assert(prove_output.key_com == verify_output.key_com);
  assert(com_vw == verifier.com_vw());

  ret = vrs::VerifySecret(prove_output.h, prove_output.g, prove_output.key_com,
                          secret_input.key_com_r, secret_input.key);
  assert(ret);

  std::cout << (ret ? "success" : "failed") << "\n";
}

inline void TestCache() {
  std::vector<bool> rets;
  // std::string output_file;
  // auto cache = vrs::CreateCache(2);
  // bool ret = vrs::SaveCache("vrs_cache", cache, output_file);
  // assert(ret);
  // rets.push_back(ret);

  // cache = vrs::CreateCache(4);
  // ret = vrs::SaveCache("vrs_cache", cache, output_file);
  // assert(ret);
  // rets.push_back(ret);

  // cache = vrs::CreateCache(8);
  // ret = vrs::SaveCache("vrs_cache", cache, output_file);
  // assert(ret);
  // rets.push_back(ret);

  // cache = vrs::CreateCache(100);
  // ret = vrs::SaveCache("vrs_cache", cache, output_file);
  // assert(ret);
  // rets.push_back(ret);

  auto pathname = vrs::SelectCacheFile("vrs_cache", 121);
  assert(!pathname.empty());

  Cache cache;
  bool ret = vrs::LoadCache(pathname, cache, true);
  assert(ret);
  rets.push_back(ret);

  vrs::UpgradeCache(cache, 121);

  ret = CheckCache(cache);
  assert(ret);
  rets.push_back(ret);

  std::cout << "summary:\n";
  for (auto i : rets) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}

}  // namespace vrs