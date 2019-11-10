#pragma once

#include "basic_types_serialize.h"
#include "misc.h"
#include "scheme_atomic_swap_vc_protocol.h"

namespace scheme::atomic_swap_vc {
// save to bin
template <typename Ar>
void serialize(Ar &ar, Request const &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("p", t.demands));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Request &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("p", t.demands));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Response const &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m),
                     ("vw", t.vw), ("vrsps", t.vrs_plain_seed),
                     ("vwcr", t.vw_com_r), ("vrs_proofs", t.vrs_proofs));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Response &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m),
                     ("vw", t.vw), ("vrsps", t.vrs_plain_seed),
                     ("vwcr", t.vw_com_r), ("vrs_proofs", t.vrs_proofs));
}

// save to json
template <typename Ar>
void serialize(Ar &ar, Receipt const &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("h", t.h), ("g", t.g), ("c", t.seed0_com));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Receipt &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("h", t.h), ("g", t.g), ("c", t.seed0_com));
}

// save to json
template <typename Ar>
void serialize(Ar &ar, Secret const &t) {
  ar &YAS_OBJECT_NVP("Secret", ("s", t.seed0), ("r", t.seed0_r));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Secret &t) {
  ar &YAS_OBJECT_NVP("Secret", ("s", t.seed0), ("r", t.seed0_r));
}
}  // namespace scheme::atomic_swap_vc
