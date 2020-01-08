#pragma once

#include "misc/misc.h"
#include "scheme_atomic_swap_protocol.h"

namespace scheme::atomic_swap {
// save to json
template <typename Ar>
void serialize(Ar &ar, Request const &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("p", t.demands));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Request &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("p", t.demands));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Response const &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m),
                     ("vw", t.vw));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Response &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m),
                     ("vw", t.vw));
}

// save to json
template <typename Ar>
void serialize(Ar &ar, Receipt const &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("s", t.seed2), ("vw", t.sigma_vw),
                     ("c", t.count));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Receipt &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("s", t.seed2), ("vw", t.sigma_vw),
                     ("c", t.count));
}

// save to json
template <typename Ar>
void serialize(Ar &ar, Secret const &t) {
  ar &YAS_OBJECT_NVP("Secret", ("s", t.seed0));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Secret &t) {
  ar &YAS_OBJECT_NVP("Secret", ("s", t.seed0));
}
}  // namespace scheme::atomic_swap
