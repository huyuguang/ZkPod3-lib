#pragma once

#include <algorithm>
#include <array>
#include <boost/endian/conversion.hpp>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <numeric>
#include <vector>

#include "./funcs.h"
#include "./multiexp.h"
#include "./types.h"
#include "log/tick.h"

// pedersen commitment base H&G
class PcBase : boost::noncopyable {
 public:
#ifdef _DEBUG
  static inline int64_t const kGSize = 1024 * 128;
#else
  static inline int64_t const kGSize = 1024 * 128;
#endif

  PcBase(std::string const& file) {
    LoadInternal(file);
    BuildSigmaG();
  }

  PcBase() {
    Create();
    BuildSigmaG();
  }

  G1 const& u() const& { return u_; }
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

  G1 ComputeSigmaG(int64_t offset, uint64_t count) const {
    if (offset < 0 || offset >= kGSize || (offset + count) > kGSize) {
      throw std::invalid_argument(std::to_string(offset) + "," +
                                  std::to_string(count));
    }
    int64_t index1 = (offset + kSigmaGInterval - 1) / kSigmaGInterval;
    int64_t index2 = (offset + count) / kSigmaGInterval;
    G1 ret = std::accumulate(sigma_g_.begin() + index1,
                             sigma_g_.begin() + index2, G1Zero());
    ret = std::accumulate(g_.begin() + index2 * kSigmaGInterval,
                          g_.begin() + offset + count, ret);
    ret = std::accumulate(g_.begin() + offset,
                          g_.begin() + index1 * kSigmaGInterval, ret);

    assert(ret == std::accumulate(g_.begin() + offset,
                                  g_.begin() + offset + count, G1Zero()));
    return ret;
  }

  typedef std::function<Fr const&(int64_t i)> GetX;
  G1 ComputeCommitmentG(int64_t g_offset, int64_t n, GetX const& x, Fr const& r,
                        bool check_01) const {
    if (g_offset < 0 || (kGSize < g_offset + n)) {
      assert(false);
      throw std::runtime_error(std::to_string(g_offset) + " + " +
                               std::to_string(n) + " too large");
    }
    auto get_g = [this, g_offset](int64_t i) -> G1 const& {
      return i ? g_[g_offset + i - 1] : h_;
    };
    auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x(i - 1) : r; };
    return MultiExpBdlo12<G1>(get_g, get_f, n + 1, check_01);
  }

  G1 ComputeCommitmentG(int64_t g_offset, std::vector<Fr> const& x, Fr const& r,
                        bool check_01) const {
    if (g_offset < 0 || (kGSize < g_offset + (int64_t)x.size())) {
      assert(false);
      throw std::runtime_error(std::to_string(g_offset) + " + " +
                               std::to_string(x.size()) + " too large");
    }
    auto get_g = [this, g_offset](int64_t i) -> G1 const& {
      return i ? g_[g_offset + i - 1] : h_;
    };
    auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x[i - 1] : r; };
    return MultiExpBdlo12<G1>(get_g, get_f, x.size() + 1, check_01);
  }

  G1 ComputeCommitmentG(int64_t g_offset, Fr const& x, Fr const& r) const {
    if (g_offset == -1) return h_ * r + u_ * x;
    if (g_offset < 0 || g_offset >= kGSize) {
      assert(false);
      throw std::runtime_error(std::to_string(g_offset) + " too large");
    }
    return h_ * r + g_[g_offset] * x;
  }

  std::vector<G1> CopyG(size_t offset, size_t count) const {
    assert(offset + count <= kGSize);
    auto begin = g_.begin() + offset;
    auto end = begin + count;
    return std::vector<G1>(begin, end);
  }

 private:
  struct Header {
    uint64_t header_size;
    uint64_t g_size;
  };

  void Create() {
    Tick tick(__FUNCTION__);

    GenerateG1(0xffffffff, &h_);

    GenerateG1(0xfffffffe, &u_);

    auto parallel_f = [this](int64_t i) { GenerateG1(i, &g_[i]); };
    parallel::For(kGSize, parallel_f);
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

    if (!WriteG1(f, u_)) {
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

    if (!ReadG1(f, u_)) throw std::runtime_error("Read pedersen_base failed");

    for (auto& i : g_) {
      if (!ReadG1(f, i)) throw std::runtime_error("Read pedersen_base failed");
    }
  }

  void BuildSigmaG() {
    static_assert(kGSize % kSigmaGInterval == 0 && kGSize >= kSigmaGInterval,
                  "");
    sigma_g_.resize(kGSize / kSigmaGInterval);

    auto parallel_f = [this](int64_t i) {
      auto begin = g_.data() + i * kSigmaGInterval;
      auto end = begin + kSigmaGInterval;
      sigma_g_[i] = std::accumulate(begin, end, G1Zero());
    };
    parallel::For(sigma_g_.size(), parallel_f);
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
  G1 u_;
  G1 h_;
  std::array<G1, kGSize> g_;
  std::vector<G1> sigma_g_;
};

inline PcBase& GetPcBase(std::string const& file = "") {
  static std::unique_ptr<PcBase> _instance_(new PcBase(file));
  return *_instance_;
}

inline bool operator==(PcBase const& a, PcBase const& b) {
  if (a.h() != b.h()) return false;
  if (a.g() != b.g()) return false;
  return true;
}

inline bool operator!=(PcBase const& a, PcBase const& b) { return !(a == b); }

inline bool OpenOrCreatePdsPub(std::string const& file) {
  auto Load = [](std::string const& file) {
    try {
      GetPcBase(file);
      return true;
    } catch (std::exception&) {
      return false;
    }
  };

  if (Load(file)) return true;

  try {
    boost::system::error_code ec;
    boost::filesystem::remove(file, ec);

    std::unique_ptr<PcBase> base(new PcBase);
    if (!base->Save(file)) {
      std::cerr << "Save ecc pds file" << file << " failed.\n";
      return false;
    }

    std::cout << "Create pc base file success.\n";
    if (!Load(file)) return false;
    assert(*base == GetPcBase());
    return true;
  } catch (std::exception& e) {
    std::cerr << "Create pc base file exception: " << e.what() << "\n";
    return false;
  }
}

inline G1 PcComputeSigmaG(int64_t offset, uint64_t count) {
  auto const& base = GetPcBase();
  return base.ComputeSigmaG(offset, count);
}

inline G1 PcComputeCommitmentG(int64_t g_offset, int64_t n,
                               PcBase::GetX const& x, Fr const& r,
                               bool check_01 = false) {
  auto const& base = GetPcBase();
  return base.ComputeCommitmentG(g_offset, n, x, r, check_01);
}

inline G1 PcComputeCommitmentG(int64_t g_offset, int64_t n, Fr const* x,
                               Fr const& r, bool check_01 = false) {
  auto const& base = GetPcBase();
  auto get_x = [x](int64_t i) -> Fr const& { return x[i]; };
  return base.ComputeCommitmentG(g_offset, n, get_x, r, check_01);
}

inline G1 PcComputeCommitmentG(int64_t g_offset, std::vector<Fr> const& x,
                               Fr const& r, bool check_01 = false) {
  auto const& base = GetPcBase();
  return base.ComputeCommitmentG(g_offset, x, r, check_01);
}

inline G1 PcComputeCommitmentG(int64_t g_offset, Fr const& x, Fr const& r) {
  auto const& base = GetPcBase();
  return base.ComputeCommitmentG(g_offset, x, r);
}

inline G1 const& PcG(int64_t i) {
  auto const& base = GetPcBase();
  if (i == -1) {
    return base.u();
  } else {
    return base.g()[i];
  }
}

inline G1 const& PcH() {
  auto const& base = GetPcBase();
  return base.h();
}

inline G1 const& PcHG(int64_t i) {
  if (i == 0) return PcH();
  return PcG(i - 1);
}