#pragma once

#include <memory>
#include <string>
#include <vector>

#include "ecc/ecc.h"
#include "scheme_ot_vrfq_protocol.h"
#include "scheme_table_alice_data.h"

namespace scheme::table::ot_vrfq {

class Alice {
 public:
  // The self_id and peer_id are useless now, just for later convenience.
  Alice(AliceDataPtr a, h256_t const& self_id, h256_t const& peer_id);

 public:
  void GetNegoReqeust(NegoARequest& request);
  bool OnNegoRequest(NegoBRequest const& request, NegoBResponse& response);
  bool OnNegoResponse(NegoAResponse const& response);

 public:
  bool OnRequest(Request const& request, Response& response);
  bool OnReceipt(Receipt const& receipt, Secret& secret);

 private:
  AliceDataPtr a_;
  h256_t const self_id_;
  h256_t const peer_id_;

 private:
  Fr r_;
  G1 g_exp_r_;

 private:
  G2 ot_self_pk_;
  G1 ot_peer_pk_;
  G2 ot_sk_;
  Fr ot_alpha_;
};

typedef std::shared_ptr<Alice> AlicePtr;
}  // namespace scheme::table::ot_vrfq