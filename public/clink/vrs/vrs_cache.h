#pragma once

#include "../details.h"
#include "./vrs_pub.h"

namespace clink {

template <typename Scheme>
struct VrsCache {
  inline static const std::string kExtensionUsing = ".using";
  inline static const std::string kExtensionUsed = ".used";

  struct FileData {
    std::string type;
    int64_t max_unit_per_zkp;
    int64_t count = 0;
    h256_t seed;  // plain seed
    Fr key;
    Fr key_com_r;
    std::vector<std::vector<G1>> var_coms;
    std::vector<std::vector<Fr>> var_coms_r;

    bool operator==(FileData const& b) const {
      return type == b.type && max_unit_per_zkp == b.max_unit_per_zkp &&
             count == b.count && seed == b.seed && key == b.key &&
             key_com_r == b.key_com_r && var_coms == b.var_coms &&
             var_coms_r == b.var_coms_r;
    }

    bool operator!=(FileData const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("cache", ("type", type), ("max", max_unit_per_zkp),
                         ("c", count), ("s", seed), ("k", key),
                         ("kcr", key_com_r), ("vc", var_coms),
                         ("vcr", var_coms_r));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("cache", ("type", type), ("max", max_unit_per_zkp),
                         ("c", count), ("s", seed), ("k", key),
                         ("kcr", key_com_r), ("vc", var_coms),
                         ("vcr", var_coms_r));
    }
  };

  typedef std::unique_ptr<FileData> FileDataUPtr;

  static bool CheckVarComs(h256_t const& seed, Fr const& key,
                           Fr const& key_com_r, int64_t begin, int64_t end,
                           std::vector<G1> const& var_coms,
                           std::vector<Fr> const& var_coms_r) {
    constexpr int64_t kPrimaryInputSize = Scheme::kPrimaryInputSize;
    auto count = end - begin;
    std::vector<Fr> v(count);
    std::vector<std::vector<Fr>> values(count);
    Scheme scheme;
    libsnark::protoboard<Fr> pb;
    auto gadget = scheme.CreateGadget(pb);
    auto num_var = (int64_t)pb.num_variables();
    if ((int64_t)var_coms.size() != num_var) {
      assert(false);
      return false;
    }
    if ((int64_t)var_coms_r.size() != num_var) {
      assert(false);
      return false;
    }

    for (int64_t i = 0; i < count; ++i) {
      auto plain = VrsPub<Scheme>::GeneratePlain(seed, i + begin);
      gadget->Assign(plain, key);
      assert(pb.is_satisfied());
      v[i] = pb.val(gadget->result());
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
      auto check_var_com = pc::ComputeCom(data, var_com_r);
      if (var_com != check_var_com) {
        assert(false);
        return false;
      }
    }
    return true;
  }

  static bool CheckFileData(FileData const& cache) {
    static constexpr int64_t kPrimaryInputSize = Scheme::kPrimaryInputSize;
    if (!cache.count) return false;
    auto count = cache.count;
    auto items = VrsPub<Scheme>::SplitLargeTask(count);
    if (cache.var_coms.size() != items.size()) return false;
    if (cache.var_coms_r.size() != items.size()) return false;

    Fr check_key_com_r = parallel::Accumulate(
        cache.var_coms_r.begin(), cache.var_coms_r.end(), FrZero(),
        [](Fr const& a, std::vector<Fr> const& b) {
          return a + b[kPrimaryInputSize];
        });

    if (check_key_com_r != cache.key_com_r) return false;

    auto const& var_coms = cache.var_coms;
    auto const& var_coms_r = cache.var_coms_r;

    bool all_success = false;
    auto f = [&var_coms_r, &items, &cache, &var_coms](int64_t i) {
      auto const& item = items[i];
      auto const& key_com_r = var_coms_r[i][kPrimaryInputSize];
      return CheckVarComs(cache.seed, cache.key, key_com_r, item.begin,
                          item.end, var_coms[i], var_coms_r[i]);
    };
    parallel::For(&all_success, items.size(), f);

    if (!all_success) {
      return false;
    }
    return true;
  }

  static void ComputeVarComs(h256_t const& seed, Fr const& key,
                             Fr const& key_com_r, int64_t begin, int64_t end,
                             std::vector<G1>& var_coms,
                             std::vector<Fr>& var_coms_r) {
    Tick tick(__FN__);
    int64_t const kPrimaryInputSize = Scheme::kPrimaryInputSize;
    auto count = end - begin;
    std::vector<Fr> v(count);
    std::vector<std::vector<Fr>> values(count);
    Scheme scheme;
    libsnark::protoboard<Fr> pb;
    auto gadget = scheme.CreateGadget(pb);

    for (int64_t i = 0; i < count; ++i) {
      auto plain = VrsPub<Scheme>::GeneratePlain(seed, i + begin);
      gadget->Assign(plain, key);
      assert(pb.is_satisfied());
      v[i] = pb.val(gadget->result());
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
      var_com = pc::ComputeCom(data, var_com_r);
    }
  }

  static void UpgradeVarComs(h256_t const& seed, Fr const& key, int64_t begin,
                             int64_t old_end, int64_t new_end,
                             std::vector<G1>& var_coms) {
    Tick tick(__FN__);
    if (old_end == new_end) return;
    // int64_t const kPrimaryInputSize = Scheme::kPrimaryInputSize;
    auto min_end = std::min(old_end, new_end);
    auto max_end = std::max(old_end, new_end);
    bool is_grow = new_end > old_end;
    auto count = max_end - min_end;
    std::vector<std::vector<Fr>> values(count);
    Scheme scheme;
    libsnark::protoboard<Fr> pb;
    auto gadget = scheme.CreateGadget(pb);

    for (int64_t i = 0; i < count; ++i) {
      auto plain = VrsPub<Scheme>::GeneratePlain(seed, i + min_end);
      gadget->Assign(plain, key);
      assert(pb.is_satisfied());
      values[i] = pb.full_variable_assignment();
    }

    auto num_var = (int64_t)pb.num_variables();
    assert(num_var == (int64_t)var_coms.size());
    var_coms.resize(num_var);
    auto const* g = pc::PcG() + min_end - begin;
    auto parallel_f = [&var_coms, &values, is_grow, count, g](uint64_t i) {
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
    };
    parallel::For(num_var, parallel_f);
  }

  static FileData CreateFileData(int64_t count) {
    Tick tick(__FN__);
    FileData cache;
    cache.count = count;
    cache.seed = misc::RandH256();
    cache.key = FrRand();
    cache.type = Scheme::type();
    cache.max_unit_per_zkp = Scheme::kMaxUnitPerZkp;

    auto items = VrsPub<Scheme>::SplitLargeTask(cache.count);

    std::vector<Fr> key_com_rs(items.size());
    for (auto& i : key_com_rs) i = FrRand();
    cache.key_com_r =
        parallel::Accumulate(key_com_rs.begin(), key_com_rs.end(), FrZero());

    cache.var_coms.resize(items.size());
    cache.var_coms_r.resize(items.size());
    auto& var_coms = cache.var_coms;
    auto& var_coms_r = cache.var_coms_r;

    auto f = [&items, &cache, &key_com_rs, &var_coms, &var_coms_r](int64_t i) {
      auto const& item = items[i];
      ComputeVarComs(cache.seed, cache.key, key_com_rs[i], item.begin, item.end,
                     var_coms[i], var_coms_r[i]);
    };
    parallel::For(items.size(), f);

    assert(CheckFileData(cache));
    return cache;
  }

  static bool LoadFile(std::string const& pathname, FileData& cache,
                       bool check_name) {
    Tick tick(__FN__);
    try {
      yas::file_istream is(pathname.c_str());
      yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
      // yas::json_iarchive<yas::file_istream> ia(is);
      ia.serialize(cache);
    } catch (std::exception& e) {
      std::cerr << __FN__ << ": " << e.what() << "\n";
      boost::system::error_code ec;
      fs::remove(pathname, ec);
      return false;
    }

    if (cache.max_unit_per_zkp != Scheme::kMaxUnitPerZkp ||
        cache.type != Scheme::type()) {
      std::cout << "invalid max_unit_per_zkp or type\n";
      assert(false);
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

  static bool SaveFileData(std::string const& cache_dir, FileData const& cache,
                           std::string& output) {
    Tick tick(__FN__);
    if (cache.max_unit_per_zkp != Scheme::kMaxUnitPerZkp ||
        cache.type != Scheme::type()) {
      std::cout << "invalid max_unit_per_zkp or type\n";
      assert(false);
      return false;
    }

    std::string base_name =
        std::to_string(cache.count) + "_" + misc::HexToStr(cache.seed);
    std::string temp_name = base_name + ".tmp";
    std::string temp_path_name = cache_dir + "/" + temp_name;
    output = cache_dir + "/" + base_name;

    try {
      yas::file_ostream os(temp_path_name.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      // yas::json_oarchive<yas::file_ostream> oa(os);
      oa.serialize(cache);
    } catch (std::exception& e) {
      std::cerr << __FN__ << ": " << e.what() << "\n";
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
    FileData check_cache;
    assert(LoadFile(output, check_cache, true) && check_cache == cache);
#endif
    return true;
  }

  static std::string SelectFile(std::string const& cache_dir, int64_t count) {
    boost::system::error_code ec;
    if (!fs::is_directory(cache_dir, ec)) return "";

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
    auto range =
        boost::make_iterator_range(fs::directory_iterator(cache_dir), {});
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

    auto const kMaxUnitPerZkp = Scheme::kMaxUnitPerZkp;
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
      auto ori_name = cache_dir + "/" + it->name;
      auto using_name = ori_name + kExtensionUsing;
      // auto using_name = fs::change_extension(ori_name, "use");

      try {
        fs::rename(ori_name, using_name);
        return using_name;
      } catch (std::exception& e) {
        std::cerr << __FN__ << ":" << e.what() << "\n";
        files.erase(it);
      }
    }
    return "";
  }

  static void ReturnFile(std::string const& using_name) {
    assert(fs::extension(using_name) == kExtensionUsing);
    auto ori_name = fs::change_extension(using_name, "");
    try {
      fs::rename(using_name, ori_name);
    } catch (std::exception&) {
      // ignore the error
    }
  }

  static void ExhaustFile(std::string const& using_name) {
    assert(fs::extension(using_name) == kExtensionUsing);
    auto used_name = fs::change_extension(using_name, kExtensionUsed);
    try {
      fs::rename(using_name, used_name);
    } catch (std::exception&) {
      // ignore the error
    }
  }

  static void UpgradeFileData(FileData& cache, int64_t count) {
    Tick tick(__FN__);
    static constexpr int64_t kPrimaryInputSize = 1;
    if (cache.count == count) return;
    auto old_count = cache.count;
    auto new_count = count;
    auto old_items = VrsPub<Scheme>::SplitLargeTask(old_count);
    auto new_items = VrsPub<Scheme>::SplitLargeTask(new_count);
    std::cout << "old_count: " << old_count
              << ", old_items: " << old_items.size() << "\n";
    std::cout << "new_count: " << new_count
              << ", new_items: " << new_items.size() << "\n";
    if (old_items.size() <= new_items.size()) {
      auto const& old_item = old_items.back();
      auto const& new_item = new_items[old_items.size() - 1];
      assert(old_item.begin == new_item.begin);
      auto& var_com = cache.var_coms.back();
      UpgradeVarComs(cache.seed, cache.key, old_item.begin, old_item.end,
                     new_item.end, var_com);

      std::cout << "add " << new_items.size() - old_items.size() << " items\n";
      cache.var_coms.resize(new_items.size());
      cache.var_coms_r.resize(new_items.size());
      std::vector<Fr> key_com_rs(new_items.size() - old_items.size());
      auto parallel_f = [&cache, &key_com_rs, &new_items,
                         old_item_size = old_items.size()](size_t i) mutable {
        auto const& add_new_item = new_items[i];
        auto& add_var_com = cache.var_coms[i];
        auto& add_var_com_r = cache.var_coms_r[i];
        auto& key_com_r = key_com_rs[i - old_item_size];
        key_com_r = FrRand();
        ComputeVarComs(cache.seed, cache.key, key_com_r, add_new_item.begin,
                       add_new_item.end, add_var_com, add_var_com_r);
      };
      parallel::For(old_items.size(), new_items.size(), parallel_f);
      cache.key_com_r = parallel::Accumulate(key_com_rs.begin(),
                                             key_com_rs.end(), cache.key_com_r);
    } else {
      auto const& new_item = new_items.back();
      auto const& old_item = old_items[new_items.size() - 1];
      assert(old_item.begin == new_item.begin);
      auto& var_com = cache.var_coms[new_items.size() - 1];
      UpgradeVarComs(cache.seed, cache.key, old_item.begin, old_item.end,
                     new_item.end, var_com);
      for (size_t i = new_items.size(); i < old_items.size(); ++i) {
        auto& add_var_com_r = cache.var_coms_r[i];
        auto const& key_com_r = add_var_com_r[kPrimaryInputSize];
        cache.key_com_r -= key_com_r;
      }
      cache.var_coms.resize(new_items.size());
      cache.var_coms_r.resize(new_items.size());
    }

    cache.count = count;

    assert(CheckFileData(cache));
  }

  static bool CreateAndSave(std::string const& data_dir, int64_t count) {
    std::string cache_dir = data_dir + "/vrs_cache/";
    cache_dir += Scheme::type();

    fs::create_directories(cache_dir);
    if (!fs::is_directory(cache_dir)) {
      std::cerr << "create directory failed: " << cache_dir << "\n";
      return false;
    }
    auto cache = CreateFileData(count);
    std::string cache_file;
    bool ret = SaveFileData(cache_dir, cache, cache_file);
    if (ret) {
      std::cout << "Create vrs cache Success: cache_file: " << cache_file
                << "\n";
    } else {
      assert(false);
      std::cout << "Save vrs cache failed: " << cache_file << "\n";
    }
    return ret;
  }

  class AutoFile : boost::noncopyable {
   public:
    AutoFile(std::string cache_dir, int64_t count)
        : cache_dir_(std::move(cache_dir)), count_(count) {}
    ~AutoFile() {
      if (cache_file_.empty()) return;
      if (!leaked_) {
        ReturnFile(cache_file_);
      } else {
        ExhaustFile(cache_file_);
      }
    }
    FileDataUPtr LoadAndUpgrade() {
      FileDataUPtr cache;
      if (cache_dir_.empty()) return cache;
      cache_file_ = SelectFile(cache_dir_, count_);
      if (!cache_file_.empty()) {
        cache.reset(new FileData());
        if (!LoadFile(cache_file_, *cache, false)) {
          std::cerr << "LoadFile failed\n";
          cache.reset();
          cache_file_.clear();
        } else {
          UpgradeFileData(*cache, count_);
          // CheckFileData(*cache);
        }
      }
      return cache;
    }

    void SetLeaked() { leaked_ = true; }

   private:
    std::string cache_dir_;
    int64_t count_;
    std::string cache_file_;
    bool leaked_{false};
  };
};
}  // namespace clink