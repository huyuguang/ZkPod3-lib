#pragma once

#include <memory>
#include <string>

#include "misc/misc.h"
#include "ecc/ecc.h"
#include "scheme/public_misc.h"
#include "scheme_ot_complaint_protocol.h"

namespace scheme::ot_complaint {

template <typename BobData>
class Bob {
 public:
  typedef std::shared_ptr<BobData> BobDataPtr;
  Bob(BobDataPtr b, h256_t const& self_id, h256_t const& peer_id,
      std::vector<Range> demands, std::vector<Range> phantoms);

 public:
  void GetNegoReqeust(NegoBRequest& request);
  bool OnNegoRequest(NegoARequest const& request, NegoAResponse& response);
  bool OnNegoResponse(NegoBResponse const& response);

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
  BobDataPtr b_;
  h256_t const self_id_;
  h256_t const peer_id_;
  uint64_t const n_;
  uint64_t const s_;
  std::vector<Range> const demands_;
  std::vector<Range> const phantoms_;
  uint64_t demands_count_ = 0;
  uint64_t phantoms_count_ = 0;

 private:
  Request request_;
  std::vector<G1> k_;      // sizeof() = L
  std::vector<G1> ot_ui_;  // sizeof() = K
  h256_t alice_nonce_;
  h256_t bob_nonce_;
  Receipt receipt_;

 private:
  struct Mapping {
    uint64_t phantom_offset;
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

 private:
  G1 ot_self_pk_;
  G2 ot_peer_pk_;
  G1 ot_sk_;
  Fr ot_beta_;
  Fr ot_rand_a_;
  Fr ot_rand_b_;
};

template <typename BobData>
using BobPtr = std::shared_ptr<Bob<BobData>>;

}  // namespace scheme::ot_complaint

#include "scheme_ot_complaint_bob.inc"