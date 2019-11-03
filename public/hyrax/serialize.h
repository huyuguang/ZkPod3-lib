#pragma once

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100)
#pragma warning(disable : 4101)
#pragma warning(disable : 4702)
#endif
#include <yas/binary_iarchive.hpp>
#include <yas/binary_oarchive.hpp>
#include <yas/file_streams.hpp>
#include <yas/json_iarchive.hpp>
#include <yas/json_oarchive.hpp>
#include <yas/mem_streams.hpp>
#include <yas/serialize.hpp>
#include <yas/std_types.hpp>
#include <yas/types/std/array.hpp>
#include <yas/types/std/pair.hpp>
#include <yas/types/std/string.hpp>
#include <yas/types/std/vector.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "a2.h"
#include "basic_types_serialize.h"

namespace hyrax::a2 {
// save to bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub const &t) {
  ar &YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, CommitmentExtPub &t) {
  ar &YAS_OBJECT_NVP("a2.cep", ("delta", t.delta), ("beta", t.beta));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, Proof const &t) {
  ar &YAS_OBJECT_NVP("a2.pf", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, Proof &t) {
  ar &YAS_OBJECT_NVP("a2.pf", ("z", t.z), ("z_delta", t.z_delta),
                     ("z_beta", t.z_beta));
}

// save to bin
template <typename Ar>
void serialize(Ar &ar, RomProof const &t) {
  ar &YAS_OBJECT_NVP("a2.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

// load from bin
template <typename Ar>
void serialize(Ar &ar, RomProof &t) {
  ar &YAS_OBJECT_NVP("a2.rp", ("c", t.com_ext_pub), ("p", t.proof));
}

}  // namespace hyrax::a2
