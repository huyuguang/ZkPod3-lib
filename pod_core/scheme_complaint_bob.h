#pragma once

#include "misc/misc.h"
#include "ecc/ecc.h"
#include "parallel/parallel.h"
#include "public.h"
#include "scheme_complaint_notary.h"
#include "scheme_complaint_protocol.h"
#include "scheme/public_misc.h"

namespace scheme::complaint {

template <typename BobData>
class Bob {
 public:
  Bob(std::shared_ptr<BobData> b, h256_t const& self_id, h256_t const& peer_id,
      std::vector<Range> demands);

 public:
  void GetRequest(Request& request);
  bool OnResponse(Response response, Receipt& receipt);
  bool OnSecret(Secret const& secret);
  bool GenerateClaim(Claim& claim);
  bool SaveDecrypted(std::string const& file);

 private:
  void BuildMapping();
  bool CheckEncryptedM();
  bool CheckK(std::vector<Fr> const& v);
  void DecryptM(std::vector<Fr> const& v);
  void BuildClaim(uint64_t i, Claim& claim);

 private:
  std::shared_ptr<BobData> b_;
  h256_t const self_id_;
  h256_t const peer_id_;
  uint64_t const n_;
  uint64_t const s_;
  std::vector<Range> const demands_;
  uint64_t demands_count_ = 0;
  h256_t alice_nonce_;
  h256_t bob_nonce_;
  Receipt receipt_;

 private:
  std::vector<G1> k_;

 private:
  struct Mapping {
    uint64_t global_index;
  };
  std::vector<Mapping> mappings_;

 private:
  h256_t seed2_;
  std::vector<Fr> w_;  // size() is L
  h256_t k_mkl_root_;
  std::vector<Fr> decrypted_m_;
  std::vector<Fr> encrypted_m_;
  int64_t claim_i_ = -1;
};

template <typename BobData>
using BobPtr = std::shared_ptr<Bob<BobData>>;
}  // namespace scheme::complaint

#include "scheme_complaint_bob.inc"