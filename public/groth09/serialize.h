#pragma once

#include "basic_types_serialize.h"
#include "hyrax/serialize.h"
#include "sec43.h"
#include "sec51.h"
#include "sec53.h"

namespace groth09::sec51 {
// save to bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub const &t) {
  ar &YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub &t) {
  ar &YAS_OBJECT_NVP("51.cep", ("ad", t.ad), ("bd", t.bd), ("c1", t.c1),
                     ("c0", t.c0));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Proof const &t) {
  ar &YAS_OBJECT_NVP("51.pf", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Proof &t) {
  ar &YAS_OBJECT_NVP("51.pf", ("fx", t.fx), ("fy", t.fy), ("rx", t.rx),
                     ("sy", t.sy), ("tz", t.tz));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, RomProof const &t) {
  ar &YAS_OBJECT_NVP("51.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, RomProof &t) {
  ar &YAS_OBJECT_NVP("51.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

}  // namespace groth09::sec51

namespace groth09::sec53 {
// save to bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub const &t) {
  ar &YAS_OBJECT_NVP("53.cep", ("cl", t.cl), ("cu", t.cu));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub &t) {
  ar &YAS_OBJECT_NVP("53.cep", ("cl", t.cl), ("cu", t.cu));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, RomProof const &t) {
  ar &YAS_OBJECT_NVP("53.rp", ("c", t.com_ext_pub), ("r", t.rom_proof_51));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, RomProof &t) {
  ar &YAS_OBJECT_NVP("53.rp", ("c", t.com_ext_pub), ("r", t.rom_proof_51));
}

}  // namespace groth09::sec53

namespace groth09::sec43 {

// save to bin
template <typename Ar>
void serialize(Ar &ar, RomProof const &t) {
  ar &YAS_OBJECT_NVP("43.rp", ("c", t.c), ("53p", t.proof_53),
                     ("a2p", t.proof_a2));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, RomProof &t) {
  ar &YAS_OBJECT_NVP("43.rp", ("c", t.c), ("53p", t.proof_53),
                     ("a2p", t.proof_a2));
}

}  // namespace groth09::sec43