#pragma once

#include "public.h"
#include "serialize.h"
#include "vrs_mimc5_gadget.h"
#include "vrs_misc.h"

namespace vrs {

inline static const std::string kExtensionUsing = ".using";
inline static const std::string kExtensionUsed = ".used";

inline bool CheckVarComs(h256_t const& seed, Fr const& key, Fr const& key_com_r,
                         int64_t begin, int64_t end,
                         std::vector<G1> const& var_coms,
                         std::vector<Fr> const& var_coms_r) {
  using groth09::details::ComputeCommitment;
  static constexpr int64_t kPrimaryInputSize = 1;
  auto count = end - begin;
  std::vector<Fr> v(count);
  std::vector<std::vector<Fr>> values(count);
  libsnark::protoboard<Fr> pb;
  Mimc5Gadget gadget(pb);
  pb.set_input_sizes(kPrimaryInputSize);
  auto num_var = (int64_t)pb.num_variables();
  if ((int64_t)var_coms.size() != num_var) return false;
  if ((int64_t)var_coms_r.size() != num_var) return false;

  for (int64_t i = 0; i < count; ++i) {
    auto plain = GeneratePlain(seed, i + begin);
    gadget.Assign(plain, key);
    assert(pb.is_satisfied());
    v[i] = pb.val(gadget.result());
    values[i] = pb.full_variable_assignment();
  }

  for (int64_t i = 0; i < num_var; ++i) {
    auto& var_com = var_coms[i];
    auto& var_com_r = var_coms_r[i];
    std::vector<Fr> data(count);
    for (int64_t j = 0; j < count; ++j) {
      auto const& values_j = values[j];
      data[j] = values_j[i];
    }
    if (i < kPrimaryInputSize) {
      if (var_com_r != FrZero()) {
        assert(false);
        return false;
      }
    } else if (i == kPrimaryInputSize) {
      if (var_com_r != key_com_r) {
        assert(false);
        return false;
      }
    }
    auto check_var_com = ComputeCommitment(data, var_com_r);
    if (var_com != check_var_com) {
      assert(false);
      return false;
    }
  }
  return true;
}

inline bool CheckCache(Cache const& cache) {
  static constexpr int64_t kPrimaryInputSize = 1;
  if (!cache.count) return false;
  auto count = cache.count;
  auto items = SplitLargeTask(count);
  if (cache.var_coms.size() != items.size()) return false;
  if (cache.var_coms_r.size() != items.size()) return false;

  Fr check_key_com_r =
      std::accumulate(cache.var_coms_r.begin(), cache.var_coms_r.end(),
                      FrZero(), [](Fr const& a, std::vector<Fr> const& b) {
                        return a + b[kPrimaryInputSize];
                      });

  if (check_key_com_r != cache.key_com_r) return false;

  auto const& var_coms = cache.var_coms;
  auto const& var_coms_r = cache.var_coms_r;

  std::vector<int64_t> rets(items.size());
#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (int64_t i = 0; i < (int64_t)items.size(); ++i) {
    auto const& item = items[i];
    auto const& key_com_r = var_coms_r[i][kPrimaryInputSize];
    rets[i] = CheckVarComs(cache.seed, cache.key, key_com_r, item.first,
                           item.second, var_coms[i], var_coms_r[i]);
  }

  if (std::any_of(rets.begin(), rets.end(), [](int64_t r) { return !r; })) {
    assert(false);
    return false;
  }
  return true;
}

inline void ComputeVarComs(h256_t const& seed, Fr const& key,
                           Fr const& key_com_r, int64_t begin, int64_t end,
                           std::vector<G1>& var_coms,
                           std::vector<Fr>& var_coms_r) {
  using groth09::details::ComputeCommitment;
  static constexpr int64_t kPrimaryInputSize = 1;
  auto count = end - begin;
  std::vector<Fr> v(count);
  std::vector<std::vector<Fr>> values(count);
  libsnark::protoboard<Fr> pb;
  Mimc5Gadget gadget(pb);
  pb.set_input_sizes(kPrimaryInputSize);

  for (int64_t i = 0; i < count; ++i) {
    auto plain = GeneratePlain(seed, i + begin);
    gadget.Assign(plain, key);
    assert(pb.is_satisfied());
    v[i] = pb.val(gadget.result());
    values[i] = pb.full_variable_assignment();
  }

  auto num_var = (int64_t)pb.num_variables();
  var_coms.resize(num_var);
  var_coms_r.resize(num_var);
  for (int64_t i = 0; i < num_var; ++i) {
    auto& var_com = var_coms[i];
    auto& var_com_r = var_coms_r[i];
    std::vector<Fr> data(count);
    for (int64_t j = 0; j < count; ++j) {
      auto const& values_j = values[j];
      data[j] = values_j[i];
    }
    if (i < kPrimaryInputSize) {
      var_com_r = FrZero();
    } else if (i == kPrimaryInputSize) {
      var_com_r = key_com_r;
    } else {
      var_com_r = FrRand();
    }
    var_com = ComputeCommitment(data, var_com_r);
  }
}

inline void UpgradeVarComs(h256_t const& seed, Fr const& key, int64_t begin,
                           int64_t old_end, int64_t new_end,
                           std::vector<G1>& var_coms) {
  if (old_end == new_end) return;
  using groth09::details::ComputeCommitment;
  static constexpr int64_t kPrimaryInputSize = 1;
  auto min_end = std::min(old_end, new_end);
  auto max_end = std::max(old_end, new_end);
  bool is_grow = new_end > old_end;
  auto count = max_end - min_end;
  std::vector<std::vector<Fr>> values(count);
  libsnark::protoboard<Fr> pb;
  Mimc5Gadget gadget(pb);
  pb.set_input_sizes(kPrimaryInputSize);

  for (int64_t i = 0; i < count; ++i) {
    auto plain = GeneratePlain(seed, i + min_end);
    gadget.Assign(plain, key);
    assert(pb.is_satisfied());
    values[i] = pb.full_variable_assignment();
  }

  auto num_var = (int64_t)pb.num_variables();
  assert(num_var == (int64_t)var_coms.size());
  var_coms.resize(num_var);
  auto const* g = GetPdsPub().g().data() + min_end - begin;
  for (int64_t i = 0; i < num_var; ++i) {
    auto& var_com = var_coms[i];
    std::vector<Fr> data(count);
    for (int64_t j = 0; j < count; ++j) {
      auto const& values_j = values[j];
      data[j] = values_j[i];
    }
    auto delta = MultiExpBdlo12(g, data.data(), count);
    if (is_grow) {
      var_com += delta;
    } else {
      var_com -= delta;
    }
  }
}

inline Cache CreateCache(int64_t count) {
  Tick tick(__FUNCTION__);
  Cache cache;
  cache.count = count;
  cache.seed = misc::RandH256();
  cache.key = FrRand();

  auto items = SplitLargeTask(cache.count);

  std::vector<Fr> key_com_rs(items.size());
  for (auto& i : key_com_rs) i = FrRand();
  cache.key_com_r =
      std::accumulate(key_com_rs.begin(), key_com_rs.end(), FrZero());

  cache.var_coms.resize(items.size());
  cache.var_coms_r.resize(items.size());
  auto& var_coms = cache.var_coms;
  auto& var_coms_r = cache.var_coms_r;

#ifdef MULTICORE
#pragma omp parallel for
#endif
  for (int64_t i = 0; i < (int64_t)items.size(); ++i) {
    auto const& item = items[i];
    ComputeVarComs(cache.seed, cache.key, key_com_rs[i], item.first,
                   item.second, var_coms[i], var_coms_r[i]);
  }

  return cache;
}

inline bool LoadCache(std::string const& pathname, Cache& cache,
                      bool check_name) {
  Tick tick(__FUNCTION__);
  try {
    yas::file_istream is(pathname.c_str());
    yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
    // yas::json_iarchive<yas::file_istream> ia(is);
    ia.serialize(cache);
  } catch (std::exception& e) {
    std::cerr << __FUNCTION__ << ": " << e.what() << "\n";
    boost::system::error_code ec;
    fs::remove(pathname, ec);
    return false;
  }

  if (check_name) {
    auto base = fs::basename(pathname);
    auto check_base =
        std::to_string(cache.count) + "_" + misc::HexToStr(cache.seed);
    return base == check_base;
  }

  return true;
}

inline bool SaveCache(std::string const& dir, Cache const& cache,
                      std::string& output) {
  Tick tick(__FUNCTION__);
  std::string base_name =
      std::to_string(cache.count) + "_" + misc::HexToStr(cache.seed);
  std::string temp_name = base_name + ".tmp";
  std::string temp_path_name = dir + "/" + temp_name;
  output = dir + "/" + base_name;

  try {
    yas::file_ostream os(temp_path_name.c_str());
    yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
    // yas::json_oarchive<yas::file_ostream> oa(os);
    oa.serialize(cache);
  } catch (std::exception& e) {
    std::cerr << __FUNCTION__ << ": " << e.what() << "\n";
    fs::remove(temp_path_name);
    return false;
  }

  boost::system::error_code ec;
  fs::rename(temp_path_name, output, ec);
  if (ec) {
    fs::remove(temp_path_name);
    return false;
  }

#ifdef _DEBUG
  Cache check_cache;
  assert(LoadCache(output, check_cache, true) && check_cache == cache);
#endif
  return true;
}

inline std::string SelectCacheFile(std::string const& dir, int64_t count) {
  auto get_count_from_name = [](std::string const& name) -> int64_t {
    auto pos = name.find("_");
    if (pos == std::string::npos) return 0;
    if (name.size() != pos + 1 + 64) return 0;
    auto s = name.substr(0, pos);
    try {
      return std::stoull(s.c_str());
    } catch (std::exception&) {
      return 0;
    }
  };

  struct Item {
    Item(int64_t c, std::string n = "") : count(c), name(std::move(n)) {}
    int64_t count;
    std::string name;
    bool operator<(Item const& a) const { return count < a.count; }
    bool operator==(Item const& a) const { return count == a.count; }
  };

  std::vector<Item> files;
  auto range = boost::make_iterator_range(fs::directory_iterator(dir), {});
  for (auto& entry : range) {
    auto basename = fs::basename(entry);
    auto extension = fs::extension(entry);
    if (!extension.empty()) continue;
    auto this_count = get_count_from_name(basename);
    if (this_count == 0) continue;
    files.push_back(Item(this_count, std::move(basename)));
  }
  if (files.empty()) return "";
  std::sort(files.begin(), files.end());

  static std::mutex mutex;
  while (!files.empty()) {
    auto it = std::upper_bound(files.begin(), files.end(), Item(count));
    if (it == files.begin()) {
      // all files > count
      if (count < kMaxUnitPerZkp) {
        auto count2 = std::min(it->count, kMaxUnitPerZkp);
        if (count2 >= count * 2) return "";  // no benefit
      }
    } else if (it == files.end()) {
      // all files <= count
      --it;
    } else {
      // find the nearest one
      auto gap1 = it->count - count;
      assert(gap1 > 0);
      auto it2 = it - 1;
      auto gap2 = count - it2->count;
      assert(gap2 >= 0);
      if (gap1 >= gap2) it = it2;
    }

    // try to change the selected file extension
    auto ori_name = dir + "/" + it->name;
    auto using_name = ori_name + kExtensionUsing;
    // auto using_name = fs::change_extension(ori_name, "use");

    std::lock_guard<std::mutex> lock(mutex);
    try {
      fs::rename(ori_name, using_name);
      return using_name;
    } catch (std::exception& e) {
      std::cerr << __FUNCTION__ << ":" << e.what() << "\n";
      files.erase(it);
    }
  }
}

inline void ReturnCacheFile(std::string const& using_name) {
  assert(fs::extension(using_name) == kExtensionUsing);
  auto ori_name = fs::change_extension(using_name, "");
  try {
    fs::rename(using_name, ori_name);
  } catch (std::exception&) {
    // ignore the error
  }
}

inline void ExhaustCacheFile(std::string const& using_name) {
  assert(fs::extension(using_name) == kExtensionUsing);
  auto used_name = fs::change_extension(using_name, kExtensionUsed);
  try {
    fs::rename(using_name, used_name);
  } catch (std::exception&) {
    // ignore the error
  }
}

inline void UpgradeCache(Cache& cache, int64_t count) {
  static constexpr int64_t kPrimaryInputSize = 1;
  if (cache.count == count) return;
  auto old_count = cache.count;
  auto new_count = count;
  auto old_items = SplitLargeTask(old_count);
  auto new_items = SplitLargeTask(new_count);

  if (old_items.size() <= new_items.size()) {
    auto const& old_item = old_items.back();
    auto const& new_item = new_items[old_items.size() - 1];
    assert(old_item.first == new_item.first);
    auto& var_com = cache.var_coms.back();
    UpgradeVarComs(cache.seed, cache.key, old_item.first, old_item.second,
                   new_item.second, var_com);

    cache.var_coms.resize(new_items.size());
    cache.var_coms_r.resize(new_items.size());
    for (size_t i = old_items.size(); i < new_items.size(); ++i) {
      auto const& add_new_item = new_items[i];
      auto& add_var_com = cache.var_coms[i];
      auto& add_var_com_r = cache.var_coms_r[i];
      Fr key_com_r = FrRand();
      ComputeVarComs(cache.seed, cache.key, key_com_r, add_new_item.first,
                     add_new_item.second, add_var_com, add_var_com_r);
      cache.key_com_r += key_com_r;
    }
  } else {
    auto const& new_item = new_items.back();
    auto const& old_item = old_items[new_items.size() - 1];
    assert(old_item.first == new_item.first);
    auto& var_com = cache.var_coms[new_items.size() - 1];
    UpgradeVarComs(cache.seed, cache.key, old_item.first, old_item.second,
                   new_item.second, var_com);
    for (size_t i = new_items.size(); i < old_items.size(); ++i) {
      auto& add_var_com_r = cache.var_coms_r[i];
      auto const& key_com_r = add_var_com_r[kPrimaryInputSize];
      cache.key_com_r -= key_com_r;
    }
    cache.var_coms.resize(new_items.size());
    cache.var_coms_r.resize(new_items.size());
  }

  cache.count = count;

  assert(CheckCache(cache));
}

}  // namespace vrs