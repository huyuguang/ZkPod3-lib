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
//
// template<typename T>
// bool YasLoadJson(char const* data, size_t size, T& t) {
//  try {
//    yas::mem_istream is(data, size);
//    yas::json_iarchive<yas::mem_istream, YasBinF()> ia(is);
//    ia.serialize(t);
//    return true;
//  } catch (std::exception &e) {
//    std::cerr << e.what() << "\n";
//    return false;
//  }
//}

template <typename T>
bool YasLoadBin(std::string const &file, T &t) {
  try {
    yas::file_istream is(file.c_str());
    yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
    ia.serialize(t);
    return true;
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n";
    return false;
  }
}

template <typename T>
bool YasLoadBin(char const *data, size_t size, T &t) {
  try {
    yas::mem_istream is(data, size);
    yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
    ia.serialize(t);
    return true;
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n";
    return false;
  }
}

template <typename T>
bool YasSaveBin(std::string const &file, T const &t) {
  try {
    // boost::system::error_code dummy;
    fs::remove(file);
    yas::file_ostream os(file.c_str());
    yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
    oa.serialize(t);
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n";
    return false;
  }

  if (DEBUG_CHECK) {
    T check;
    CHECK(YasLoadBin(file, check), "");
    CHECK(t == check, "");
  }

  return true;
}

template <typename T>
bool YasSaveBin(yas::shared_buffer &data, T const &t) {
  try {
    // boost::system::error_code dummy;
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(t);
    data = os.get_shared_buffer();
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n";
    return false;
  }

  if (DEBUG_CHECK) {
    T check;
    CHECK(YasLoadBin(data.data.get(), data.size, check), "");
    CHECK(t == check, "");
  }

  return true;
}

template <typename T>
size_t YasGetBinLen(T const &t) {
  yas::shared_buffer data;
  if (!YasSaveBin(data, t)) return (size_t)-1;
  return data.size;
}