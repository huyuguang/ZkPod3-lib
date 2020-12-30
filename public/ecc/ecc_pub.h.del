#pragma once

#include "./funcs.h"
#include "./multiexp.h"
#include "./types.h"
#include "log/tick.h"
#include "public.h"

// Order of the G1 can get through Fr::getModulo(), which return
// 21888242871839275222246405745257275088548364400416034343698204186575808495617
// It's a prime number, that means for any generator u = g^xx, order of the sub
// group is same.

class EccPub : boost::noncopyable {
 public:
  static inline size_t const kU1Size = 2050;
  static inline size_t const kU2Size = 2;
  static inline size_t const kU1WmSize = 64;
  static inline size_t const kU2WmSize = kU2Size;

 public:
  G1WM const& g1_wm() const { return g1_wm_; }
  G2WM const& g2_wm() const { return g2_wm_; }
  std::array<G1, kU1Size> const& u1() const { return u1_; }
  std::array<G1WM, kU1WmSize> const& u1_wm() const { return u1_wm_; }
  std::array<G2, kU2Size> const& u2() const { return u2_; }
  std::array<G2WM, kU2WmSize> const& u2_wm() const { return u2_wm_; }

  EccPub(std::string const& file) { LoadInternal(file); }

  EccPub() { Create(); }

  bool Save(std::string const& file) {
    try {
      SaveInternal(file);
      return true;
    } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      return false;
    }
  }

  G1 PowerG1(Fr const& f) const {
    G1 ret;
    g1_wm_.mul(ret, f);
    return ret;
  }

  G2 PowerG2(Fr const& f) const {
    G2 ret;
    g2_wm_.mul(ret, f);
    return ret;
  }

  G1 PowerU1(uint64_t u_index, Fr const& f) const {
    CHECK(u_index < kU1Size, std::to_string(u_index));
    if (u_index < kU1WmSize) {
      G1WM const& wm = u1_wm_[u_index];
      G1 ret;
      wm.mul(ret, f);
      return ret;
    } else {
      return u1_[u_index] * f;
    }
  }

  G2 PowerU2(uint64_t u_index, Fr const& f) const {
    CHECK(u_index < kU2Size, std::to_string(u_index));
    if (u_index < kU2WmSize) {
      G2WM const& wm = u2_wm_[u_index];
      G2 ret;
      wm.mul(ret, f);
      return ret;
    } else {
      return u2_[u_index] * f;
    }
  }

 private:
  struct Header {
    uint64_t header_size;
    uint64_t u1_size;
    uint64_t u2_size;
    uint64_t u1wm_size;
    uint64_t u2wm_size;
    uint64_t g1wm_len;
    uint64_t g2wm_len;
    uint64_t u1wm_len;
    uint64_t u2wm_len;
  };

  void Create() {
    Tick tick(__FN__);

    auto fr_bits = Fr::getBitSize();

    for (size_t i = 0; i < kU1Size; ++i) {
      std::string seed = "pod_u1_" + std::to_string(i);
      u1_[i] = MapToG1(seed);
      u1_[i].normalize();
      if (i < kU1WmSize) {
        u1_wm_[i].init(u1_[i], fr_bits, 4);  // use 4 is ok
      }
    }

    for (size_t i = 0; i < kU2Size; ++i) {
      std::string seed = "pod_u2_" + std::to_string(i);
      u2_[i] = MapToG2(seed);
      u2_[i].normalize();
      if (i < kU2WmSize) {
        u2_wm_[i].init(u2_[i], fr_bits, 4);  // use 4 is ok
      }
    }

    g1_wm_.init(G1One(), fr_bits, 8);  // use 8
    g2_wm_.init(G2One(), fr_bits, 8);
  }

  void SaveInternal(std::string const& file) {
    Tick tick(__FN__);
    FILE* f = fopen(file.c_str(), "wb+");
    CHECK(f, file);

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    header.header_size = sizeof(Header);
    header.u1_size = kU1Size;
    header.u2_size = kU2Size;
    header.u1wm_size = kU1WmSize;
    header.u2wm_size = kU2WmSize;
    header.g1wm_len = GetG1wmFlatLen(g1_wm_);
    header.g2wm_len = GetG2wmFlatLen(g2_wm_);
    header.u1wm_len = GetG1wmFlatLen(u1_wm_[0]);
    header.u2wm_len = GetG2wmFlatLen(u2_wm_[0]);

    CHECK(WriteHeader(f, header), "");

    CHECK(WriteG1wm(f, g1_wm_), "");

    CHECK(WriteG2wm(f, g2_wm_), "");

    for (auto& i : u1_) {
      CHECK(WriteG1(f, i), "u1");
    }

    for (auto& i : u1_wm_) {
      if (!WriteG1wm(f, i)) {
        throw std::runtime_error("Write u1_wm failed");
      }
    }

    for (auto& i : u2_) {
      if (!WriteG2(f, i)) {
        throw std::runtime_error("Write u2 failed");
      }
    }

    for (auto& i : u2_wm_) {
      if (!WriteG2wm(f, i)) {
        throw std::runtime_error("Write u2_wm failed");
      }
    }
  }

  void LoadInternal(std::string const& file) {
    Tick tick(__FN__);
    FILE* f = fopen(file.c_str(), "rb");
    if (!f) throw std::runtime_error("Open failed");

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    if (!ReadHeader(f, header)) throw std::runtime_error("Read header failed");
    if (header.u1_size != kU1Size || header.u2_size != kU2Size)
      throw std::runtime_error("Invalid data");
    if (!header.g1wm_len || !header.g2wm_len || !header.u1wm_len ||
        !header.u2wm_len)
      throw std::runtime_error("Invalid data");

    if (!ReadG1wm(f, header.g1wm_len, g1_wm_)) {
      throw std::runtime_error("Read g1_wm failed");
    }

    if (!ReadG2wm(f, header.g2wm_len, g2_wm_))
      throw std::runtime_error("Read g2_wm failed");

    for (auto& i : u1_) {
      if (!ReadG1(f, i)) throw std::runtime_error("Read u1 failed");
    }

    for (auto& i : u1_wm_) {
      if (!ReadG1wm(f, header.u1wm_len, i)) {
        throw std::runtime_error("Read u1_wm failed");
      }
    }

    for (auto& i : u2_) {
      if (!ReadG2(f, i)) throw std::runtime_error("Read u2 failed");
    }

    for (auto& i : u2_wm_) {
      if (!ReadG2wm(f, header.u2wm_len, i)) {
        throw std::runtime_error("Read u2_wm failed");
      }
    }
  }

 private:
  enum {
    kFpBinSize = 32,
    kG1BufSize = 1 + kFpBinSize * 2,
    kFp2BinSize = 64,
    kG2BufSize = 1 + kFp2BinSize * 2
  };

  static void ComputePdsBase(uint64_t index, G1* g) {
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
    if (!WriteUint(f, v.u1_size)) return false;
    if (!WriteUint(f, v.u2_size)) return false;
    if (!WriteUint(f, v.u1wm_size)) return false;
    if (!WriteUint(f, v.u2wm_size)) return false;
    if (!WriteUint(f, v.g1wm_len)) return false;
    if (!WriteUint(f, v.g2wm_len)) return false;
    if (!WriteUint(f, v.u1wm_len)) return false;
    if (!WriteUint(f, v.u2wm_len)) return false;
    return true;
  }

  static bool ReadHeader(FILE* f, Header& v) {
    if (!ReadUint(f, v.header_size)) return false;
    if (v.header_size != sizeof(Header)) return false;
    if (!ReadUint(f, v.u1_size)) return false;
    if (!ReadUint(f, v.u2_size)) return false;
    if (!ReadUint(f, v.u1wm_size)) return false;
    if (!ReadUint(f, v.u2wm_size)) return false;
    if (!ReadUint(f, v.g1wm_len)) return false;
    if (!ReadUint(f, v.g2wm_len)) return false;
    if (!ReadUint(f, v.u1wm_len)) return false;
    if (!ReadUint(f, v.u2wm_len)) return false;
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

  static G1Ptr ReadG1(FILE* f) {
    G1Ptr ret = std::make_shared<G1>();
    uint8_t buf[kG1BufSize];
    if (fread(buf, kG1BufSize, 1, f) != 1) return G1Ptr();
    if (!FlatBinToG1(buf, ret.get())) return G1Ptr();
    return ret;
  }

  static bool WriteG2(FILE* f, G2 const& v) {
    uint8_t buf[kG2BufSize];
    G2ToFlatBin(buf, &v);
    return fwrite(buf, kG2BufSize, 1, f) == 1;
  }

  static bool ReadG2(FILE* f, G2& v) {
    uint8_t buf[kG2BufSize];
    if (fread(buf, kG2BufSize, 1, f) != 1) return false;
    if (!FlatBinToG2(buf, &v)) return false;
    return true;
  }

  static bool WriteG1wm(FILE* f, G1WM const& v) {
    std::vector<uint8_t> bin;
    if (!G1wmToFlatBin(v, bin)) return false;
    return fwrite(bin.data(), bin.size(), 1, f) == 1;
  }

  static bool ReadG1wm(FILE* f, uint64_t size, G1WM& v) {
    std::unique_ptr<uint8_t[]> buf(new uint8_t[size]);
    if (fread(buf.get(), size, 1, f) != 1) return false;
    return FlatBinToG1wm(buf.get(), size, v);
  }

  static bool WriteG2wm(FILE* f, G2WM const& v) {
    std::vector<uint8_t> bin;
    if (!G2wmToFlatBin(v, bin)) return false;
    return fwrite(bin.data(), bin.size(), 1, f) == 1;
  }

  static bool ReadG2wm(FILE* f, uint64_t size, G2WM& v) {
    std::unique_ptr<uint8_t[]> buf(new uint8_t[size]);
    if (fread(buf.get(), size, 1, f) != 1) return false;
    return FlatBinToG2wm(buf.get(), size, v);
  }

 private:
  G1WM g1_wm_;
  G2WM g2_wm_;
  std::array<G1, kU1Size> u1_;
  std::array<G1WM, kU1WmSize> u1_wm_;
  std::array<G2, kU2Size> u2_;
  std::array<G2WM, kU2WmSize> u2_wm_;
};

inline EccPub& GetEccPub(std::string const& file = "") {
  static std::unique_ptr<EccPub> _instance_(new EccPub(file));
  return *_instance_;
}

inline bool operator==(EccPub const& a, EccPub const& b) {
  if (a.g1_wm() != b.g1_wm()) return false;
  if (a.g2_wm() != b.g2_wm()) return false;

  if (a.u1() != b.u1()) return false;
  auto const& a_u1_wm = a.u1_wm();
  auto const& b_u1_wm = b.u1_wm();
  if (a_u1_wm.size() != b_u1_wm.size()) return false;
  for (size_t i = 0; i < a_u1_wm.size(); ++i) {
    if (a_u1_wm[i] != b_u1_wm[i]) return false;
  }

  if (a.u2() != b.u2()) return false;
  auto const& a_u2_wm = a.u2_wm();
  auto const& b_u2_wm = b.u2_wm();
  if (a_u2_wm.size() != b_u2_wm.size()) return false;
  for (size_t i = 0; i < a_u2_wm.size(); ++i) {
    if (a_u2_wm[i] != b_u2_wm[i]) return false;
  }
  return true;
}

inline bool operator!=(EccPub const& a, EccPub const& b) { return !(a == b); }

inline bool OpenOrCreateEccPub(std::string const& file) {
  auto LoadEccPub = [](std::string const& file) {
    try {
      GetEccPub(file);
      return true;
    } catch (std::exception&) {
      return false;
    }
  };

  if (LoadEccPub(file)) return true;

  try {
    boost::system::error_code ec;
    boost::filesystem::remove(file, ec);

    std::unique_ptr<EccPub> ecc_pub(new EccPub);
    if (!ecc_pub->Save(file)) {
      std::cerr << "Save ecc pub file" << file << " failed.\n";
      return false;
    }

    std::cout << "Create ecc pub file success.\n";
    if (!LoadEccPub(file)) return false;
    assert(*ecc_pub == GetEccPub());
    return true;
  } catch (std::exception& e) {
    std::cerr << "Create ecc pub file exception: " << e.what() << "\n";
    return false;
  }
}

template <typename GET_F>
G1 MultiExpU1(uint64_t count, GET_F const& get_f) {
  auto const& ecc_pub = GetEccPub();
  if (count > ecc_pub.kU1Size) {
    throw std::invalid_argument("count too large");
  }

  if (count <= 32) {
    G1 ret = G1Zero();
    for (uint64_t i = 0; i < count; ++i) {
      ret += ecc_pub.PowerU1(i, get_f(i));
    }
    return ret;
  } else {
    // if count > 32, MultiExpBdlo12Inner is faster
    auto const& u1 = ecc_pub.u1();
    auto get_g = [&u1](uint64_t i) -> G1 const& { return u1[i]; };
    return MultiExpBdlo12<G1>(get_g, get_f, (size_t)count);
  }
}
