#pragma once

#include <algorithm>
#include <array>
#include <boost/endian/conversion.hpp>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <numeric>
#include <vector>

#include "ecc.h"
#include "tick.h"

class PdsPub : boost::noncopyable {
 public:
#ifdef _DEBUG
  static inline int64_t const kGSize = 32;
#else
  static inline int64_t const kGSize = 1024 * 32;
#endif

  PdsPub(std::string const& file) {
    LoadInternal(file);
    BuildSigmaG();
  }

  PdsPub() {
    Create();
    BuildSigmaG();
  }

  G1 const& h() const& { return h_; }
  std::array<G1, kGSize> const& g() const { return g_; }

  bool Save(std::string const& file) {
    try {
      SaveInternal(file);
      return true;
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
  }

  G1 ComputeSigmaG(uint64_t count) const {
    if (count > kGSize) throw std::runtime_error("bad count");
    G1 ret;
    auto index = count / kSigmaGInterval;
    if (index == 0) {
      ret = std::accumulate(g_.data() + index * kSigmaGInterval,
                            g_.data() + count, G1Zero());
    } else {
      ret = sigma_g_[index - 1] +
            std::accumulate(g_.data() + index * kSigmaGInterval,
                            g_.data() + count, G1Zero());
    }
    assert(ret == std::accumulate(g_.data(), g_.data() + count, G1Zero()));
    return ret;
  }

 private:
  struct Header {
    uint64_t header_size;
    uint64_t g_size;
  };

  void Create() {
    Tick tick(__FUNCTION__);

    GenerateG1(0xffffffff, &h_);

    auto parallel_f = [this](int64_t i) mutable { GenerateG1(i, &g_[i]); };
    parallel::For((int64_t)kGSize, parallel_f);
  }

  void SaveInternal(std::string const& file) {
    Tick tick(__FUNCTION__);
    FILE* f = fopen(file.c_str(), "wb+");
    if (!f) throw std::runtime_error("Create file failed");

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    header.header_size = sizeof(Header);
    header.g_size = kGSize;

    if (!WriteHeader(f, header)) {
      throw std::runtime_error("Write header failed");
    }

    if (!WriteG1(f, h_)) {
      throw std::runtime_error("Write pedersen_base failed");
    }

    for (auto& i : g_) {
      if (!WriteG1(f, i)) {
        throw std::runtime_error("Write pedersen_base failed");
      }
    }
  }

  void LoadInternal(std::string const& file) {
    Tick tick(__FUNCTION__);
    FILE* f = fopen(file.c_str(), "rb");
    if (!f) throw std::runtime_error("Open failed");

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    if (!ReadHeader(f, header)) throw std::runtime_error("Read header failed");
    if (header.g_size != kGSize) throw std::runtime_error("Invalid data");

    if (!ReadG1(f, h_)) throw std::runtime_error("Read pedersen_base failed");

    for (auto& i : g_) {
      if (!ReadG1(f, i)) throw std::runtime_error("Read pedersen_base failed");
    }
  }

  void BuildSigmaG() {
    static_assert(kGSize % kSigmaGInterval == 0 && kGSize >= kSigmaGInterval,
                  "");
    sigma_g_.resize(kGSize / kSigmaGInterval);
    auto begin = g_.data();
    auto end = begin + kSigmaGInterval;
    sigma_g_[0] = std::accumulate(begin, end, G1Zero());
    for (size_t i = 1; i < sigma_g_.size(); ++i) {
      begin += kSigmaGInterval;
      end += kSigmaGInterval;
      sigma_g_[i] = std::accumulate(begin, end, sigma_g_[i - 1]);
    }
  }

  enum {
    kFpBinSize = 32,
    kG1BufSize = 1 + kFpBinSize * 2,
  };

  static void GenerateG1(uint64_t index, G1* g) {
    std::string seed = "pod_pedersen_base_" + std::to_string(index);
    MapToG1(seed, g);
    g->normalize();
  }

  template <typename T>
  static bool WriteUint(FILE* f, T v) {
    v = boost::endian::native_to_big(v);
    return fwrite(&v, sizeof(v), 1, f) == 1;
  }

  template <typename T>
  static bool ReadUint(FILE* f, T& v) {
    if (fread(&v, sizeof(v), 1, f) != 1) return false;
    v = boost::endian::big_to_native(v);
    return true;
  }

  static bool WriteHeader(FILE* f, Header const& v) {
    if (!WriteUint(f, v.header_size)) return false;
    if (!WriteUint(f, v.g_size)) return false;
    return true;
  }

  static bool ReadHeader(FILE* f, Header& v) {
    if (!ReadUint(f, v.header_size)) return false;
    if (v.header_size != sizeof(Header)) return false;
    if (!ReadUint(f, v.g_size)) return false;
    return true;
  }

  static bool WriteG1(FILE* f, G1 const& v) {
    uint8_t buf[kG1FlatBinSize];
    G1ToFlatBin(v, buf);
    return fwrite(buf, kG1FlatBinSize, 1, f) == 1;
  }

  static bool ReadG1(FILE* f, G1& v) {
    uint8_t buf[kG1BufSize];
    if (fread(buf, kG1BufSize, 1, f) != 1) return false;
    if (!FlatBinToG1(buf, &v)) return false;
    return true;
  }

 private:
  static inline size_t const kSigmaGInterval = 32;
  static_assert(kGSize % kSigmaGInterval == 0 && kGSize >= kSigmaGInterval, "");

  G1 h_;
  std::array<G1, kGSize> g_;
  std::vector<G1> sigma_g_;
};

inline PdsPub& GetPdsPub(std::string const& file = "") {
  static std::unique_ptr<PdsPub> _instance_(new PdsPub(file));
  return *_instance_;
}

inline bool operator==(PdsPub const& a, PdsPub const& b) {
  if (a.h() != b.h()) return false;
  if (a.g() != b.g()) return false;
  return true;
}

inline bool operator!=(PdsPub const& a, PdsPub const& b) { return !(a == b); }

inline bool OpenOrCreatePdsPub(std::string const& file) {
  auto Load = [](std::string const& file) {
    try {
      GetPdsPub(file);
      return true;
    } catch (std::exception&) {
      return false;
    }
  };

  if (Load(file)) return true;

  try {
    boost::system::error_code ec;
    boost::filesystem::remove(file, ec);

    std::unique_ptr<PdsPub> pub(new PdsPub);
    if (!pub->Save(file)) {
      std::cerr << "Save ecc pds file" << file << " failed.\n";
      return false;
    }

    std::cout << "Create pds pub file success.\n";
    if (!Load(file)) return false;
    assert(*pub == GetPdsPub());
    return true;
  } catch (std::exception& e) {
    std::cerr << "Create pds pub file exception: " << e.what() << "\n";
    return false;
  }
}

inline G1 ComputePdsSigmaG(uint64_t count) {
  auto const& pds_pub = GetPdsPub();
  return pds_pub.ComputeSigmaG(count);
}
//
// struct PdsBase {
//  PdsBase(G1 const& h, G1 const* gstart, int64_t count)
//      : h(h), gstart(gstart), gend(gstart + count) {
//  }
//  G1 const& h;
//  G1 const* gstart;
//  G1 const* gend;
//  int64_t size() const { return gend - gstart; }
//};

// std::unique_ptr<PdsBase> GetPdsBase(uint64_t count) const {
//  if (count > kPdsSize) throw std::runtime_error("bad count");
//  return std::unique_ptr<PdsBase>(new PdsBase(pds_h_, pds_g_.data(), count));
//}
