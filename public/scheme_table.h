#pragma once

#include <stdint.h>

#include <vector>

#include "ecc.h"
#include "scheme_misc.h"
#include "vrf.h"
#include "vrf_meta.h"

namespace scheme::table {
typedef std::vector<std::string> Record;

typedef std::vector<Record> Table;

enum Type {
  kCsv,  // now only support csv format
};

inline std::istream& operator>>(std::istream& in, Type& t) {
  std::string token;
  in >> token;
  if (token == "csv") {
    t = Type::kCsv;
  } else {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, Type const& t) {
  if (t == Type::kCsv) {
    os << "csv";
  } else {
    os.setstate(std::ios_base::failbit);
  }
  return os;
}

inline void UniqueRecords(Table& table,
                          std::vector<uint64_t> const& vrf_key_colnums) {
  Tick _tick_(__FUNCTION__);
  auto unique = [&table](uint64_t pos) {
    std::map<std::string, size_t> count;
    for (auto& record : table) {
      auto& key = record[pos];
      auto& c = count[key];
      key += "_" + std::to_string(c++);
    }
  };

  for (auto i : vrf_key_colnums) {
    unique(i);
  }
}

inline uint64_t GetRecordSize(Record const& record) {
  uint64_t len = record.size() * sizeof(uint32_t);
  for (auto const& i : record) len += i.size();
  return len;
}

inline uint64_t GetMaxRecordSize(Table const& table) {
  uint64_t max_record_len = 0;
  Record const* r = nullptr;
  for (auto const& i : table) {
    auto len = GetRecordSize(i);
    if (len > max_record_len) {
      max_record_len = len;
      r = &i;
    }
  }
  if (r) {
    std::cout << "max long record: " << max_record_len << "\n";
    for (auto const& i : *r) {
      std::cout << i << ";";
    }
  }
  return max_record_len;
}

inline Fr GetPadFr(uint32_t len) {
  uint8_t bin[31];
  misc::RandomBytes(bin, sizeof(bin));
  len = boost::endian::native_to_big(len);
  memcpy(bin, &len, sizeof(len));
  return BinToFr31(bin, bin + sizeof(bin));
}

inline bool BinToRecord(std::vector<uint8_t> const& bin, Record& record,
                        uint64_t max_item) {
  uint8_t const* p = bin.data();
  size_t left_len = bin.size();
  for (;;) {
    uint32_t item_len;
    if (left_len < sizeof(item_len)) return false;
    memcpy(&item_len, p, sizeof(item_len));
    item_len = boost::endian::big_to_native(item_len);
    p += sizeof(item_len);
    left_len -= sizeof(item_len);
    if (item_len > left_len) return false;
    record.resize(record.size() + 1);
    auto& r = record.back();
    r.assign((char*)p, item_len);
    p += item_len;
    left_len -= item_len;
    if (left_len == 0) break;
    if (record.size() >= max_item) break;
  }
  return true;
}

inline void RecordToBin(Record const& record, std::vector<uint8_t>& bin) {
  assert(bin.size() >= GetRecordSize(record));

  uint8_t* p = bin.data();
  for (auto& i : record) {
    uint32_t len = (uint32_t)i.size();
    len = boost::endian::native_to_big(len);
    memcpy(p, &len, sizeof(len));
    p += sizeof(len);
    memcpy(p, i.data(), i.size());
    p += i.size();
  }

  misc::RandomBytes(p, bin.data() + bin.size() - p);

#ifdef _DEBUG
  Record debug_record;
  if (!BinToRecord(bin, debug_record, record.size())) {
    assert(false);
  }
#endif

  assert((p - bin.data()) == (int64_t)GetRecordSize(record));
}

inline void DataToM(Table const& table, std::vector<uint64_t> columens_index,
                    uint64_t s, vrf::Sk<> const& vrf_sk, std::vector<Fr>& m) {
  Tick _tick_(__FUNCTION__);
  auto record_fr_num = s - 1 - columens_index.size();
  auto n = table.size();

  //#ifdef MULTICORE
  //#pragma omp parallel for
  //#endif
  //  for (int64_t i = 0; i < (int64_t)n; ++i) {
  //    std::vector<uint8_t> bin(31 * record_fr_num);
  //    auto const& record = table[i];
  //    auto record_size = GetRecordSize(record);
  //    auto offset = i * s;
  //    for (auto j : columens_index) {
  //      auto h = HashVrfKey(record[j], vrf_sk);
  //      m[offset++] = BinToFr31(h.data(), h.data() + 31);  // drop the last
  //      byte
  //    }
  //
  //    m[offset++] = GetPadFr((uint32_t)record_size);
  //
  //    RecordToBin(record, bin);
  //    for (uint64_t j = 0; j < record_fr_num; ++j) {
  //      uint8_t const* p = bin.data() + j * 31;
  //      m[offset++] = BinToFr31(p, p + 31);
  //    }
  //  }

  auto parallel_f = [record_fr_num, s, &m, &columens_index, &table,
                     &vrf_sk](int64_t i) mutable {
    std::vector<uint8_t> bin(31 * record_fr_num);
    auto const& record = table[i];
    auto record_size = GetRecordSize(record);
    auto offset = i * s;
    for (auto j : columens_index) {
      auto h = HashVrfKey(record[j], vrf_sk);
      m[offset++] = BinToFr31(h.data(), h.data() + 31);  // drop the last byte
    }

    m[offset++] = GetPadFr((uint32_t)record_size);

    RecordToBin(record, bin);
    for (uint64_t j = 0; j < record_fr_num; ++j) {
      uint8_t const* p = bin.data() + j * 31;
      m[offset++] = BinToFr31(p, p + 31);
    }
  };
  parallel::For((int64_t)n, parallel_f);
}

inline VrfKeyMeta const* GetKeyMetaByName(VrfMeta const& vrf_meta,
                                          std::string const& name) {
  auto const& names = vrf_meta.column_names;
  auto it_name = std::find(names.begin(), names.end(), name);
  if (it_name == names.end()) return nullptr;
  auto col_index = std::distance(names.begin(), it_name);

  for (uint64_t i = 0; i < vrf_meta.keys.size(); ++i) {
    auto const& vrf_key = vrf_meta.keys[i];
    if (vrf_meta.keys[i].column_index == (uint64_t)col_index) return &vrf_key;
  }
  return nullptr;
}

inline bool DecryptedRangeMToFile(std::string const& file, uint64_t s,
                                  VrfMeta const& vrf_meta,
                                  std::vector<Range> const& demands,
                                  std::vector<Fr> const& part_m) {
  boost::system::error_code err;
  fs::remove(file, err);

  std::ofstream out(file);
  if (!out) return false;

  std::string str;
  for (auto const& i : vrf_meta.column_names) {
    str += i + ",";
  }
  str.pop_back();
  out << str;
  out << "\n";

  uint64_t demands_count = 0;
  for (auto const& i : demands) demands_count += i.count;
  h256_t pad_fr;
  uint64_t prefix_count = vrf_meta.keys.size() + 1;
  uint64_t max_record_len = (s - prefix_count) * 31 + 1;
  std::vector<uint8_t> record_bin;
  for (size_t i = 0; i < demands_count; ++i) {
    // pad fr
    size_t j = vrf_meta.keys.size();
    FrToBin(part_m[i * s + j], pad_fr.data());
    uint32_t real_len;
    memcpy(&real_len, pad_fr.data(), sizeof(real_len));
    real_len = boost::endian::big_to_native(real_len);
    record_bin.resize(max_record_len);
    if (real_len >= record_bin.size()) {
      assert(false);
      return false;
    }

    for (j = j + 1; j < s; ++j) {
      auto row_index = j - prefix_count;
      auto prow = record_bin.data() + row_index * 31;
      FrToBin(part_m[i * s + j], prow);
      assert(prow[31] == 0);
    }

    record_bin.resize(real_len);

    Record record;
    record.reserve(vrf_meta.column_names.size());
    if (!BinToRecord(record_bin, record, vrf_meta.column_names.size())) {
      assert(false);
      return false;
    }
    assert(record.size() <= vrf_meta.column_names.size());

    str.clear();
    for (auto const& r : record) {
      str += r + ",";
    }
    if (!str.empty()) str.pop_back();

    out << str << "\n";
  }

  return true;
}

}  // namespace scheme::table