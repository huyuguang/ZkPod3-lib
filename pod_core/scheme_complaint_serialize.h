#pragma once

#include "misc/misc.h"
#include "scheme_complaint_protocol.h"

namespace scheme::complaint {
// save to bin
template <typename Ar>
void serialize(Ar &ar, Request const &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("d", t.demands));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Request &t) {
  ar &YAS_OBJECT_NVP("Request", ("b", t.bob_nonce), ("d", t.demands));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Response const &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Response &t) {
  ar &YAS_OBJECT_NVP("Response", ("a", t.alice_nonce), ("k", t.k), ("m", t.m));
}

// save to json
template <typename Ar>
void serialize(Ar &ar, Receipt const &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("s", t.seed2), ("k", t.k_mkl_root),
                     ("c", t.count));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Receipt &t) {
  ar &YAS_OBJECT_NVP("Receipt", ("s", t.seed2), ("k", t.k_mkl_root),
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

// save to json
template <typename Ar>
void serialize(Ar &ar, Claim const &t) {
  ar &YAS_OBJECT_NVP("Claim", ("i", t.i), ("k", t.ki), ("m", t.mkl_path));
}

// load from json
template <typename Ar>
void serialize(Ar &ar, Claim &t) {
  ar &YAS_OBJECT_NVP("Claim", ("i", t.i), ("k", t.ki), ("m", t.mkl_path));
}
}  // namespace scheme::complaint
