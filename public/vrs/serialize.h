#pragma once

#include "basic_types_serialize.h"
#include "groth09/serialize.h"
#include "vrs_types.h"

namespace vrs {
// save to bin
template <typename Ar>
void serialize(Ar &ar, Proof const &t) {
  ar &YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Proof &t) {
  ar &YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Cache const &t) {
  ar &YAS_OBJECT_NVP("cache", ("c", t.count), ("s", t.seed), ("k", t.key),
                     ("kcr", t.key_com_r), ("vc", t.var_coms),
                     ("vcr", t.var_coms_r));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Cache &t) {
  ar &YAS_OBJECT_NVP("cache", ("c", t.count), ("s", t.seed), ("k", t.key),
                     ("kcr", t.key_com_r), ("vc", t.var_coms),
                     ("vcr", t.var_coms_r));
}

}  // namespace vrs