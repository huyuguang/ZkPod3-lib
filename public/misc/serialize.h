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

#include "./funcs.h"
#include "./types.h"

inline constexpr size_t YasBinF() {
  return yas::options::binary | yas::options::ebig | yas::options::compacted;
}

// save
template <typename Ar>
void serialize(Ar &ar, Range const &t) {
  if (ar.type() == yas::binary) {
    ar &t.start;
    ar &t.count;
  } else {
    assert(ar.type() == yas::json);
    std::string str = Range::to_string(t);
    ar &str;
  }
}

// load
template <typename Ar>
void serialize(Ar &ar, Range &t) {
  if (ar.type() == yas::binary) {
    ar &t.start;
    ar &t.count;
  } else {
    assert(ar.type() == yas::json);
    std::string str;
    ar &str;
    t = Range::from_string(str);  // throw
  }
}

// save
template <typename Ar>
void serialize(Ar &ar, h256_t const &t) {
  if (ar.type() == yas::binary) {
    ar &t.to_array();
  } else {
    assert(ar.type() == yas::json);
    auto str = misc::HexToStr(t);
    ar &str;
  }
}

// load
template <typename Ar>
void serialize(Ar &ar, h256_t &t) {
  if (ar.type() == yas::binary) {
    ar &t.to_array();
  } else {
    assert(ar.type() == yas::json);
    std::string str;
    ar &str;
    misc::HexStrToH256(str, t);
  }
}
