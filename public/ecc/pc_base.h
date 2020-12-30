#pragma once

#include "./funcs.h"
#include "./multiexp.h"
#include "./types.h"
#include "log/tick.h"
#include "public.h"

extern bool BIG_MODE;

namespace pc {

// pedersen commitment base H&G
class Base : boost::noncopyable {
 public:
  static int64_t GSize() {
    constexpr int64_t kBig = 4096 * 1025 * 10;
    constexpr int64_t kNor = 4096 * 1025;  // 16384 * 1024 + 100;
    return BIG_MODE ? kBig : kNor;
  }

  Base(std::string const& file) {
    g_ = new G1[GSize()];
    LoadInternal(file);
  }

  Base() {
    g_ = new G1[GSize()];
    Create();
  }

  ~Base() { delete[] g_; }

  G1 const& u() const& { return u_; }
  G1 const& h() const& { return h_; }
  G1 const* g() const { return g_; }

  bool Save(std::string const& file) {
    try {
      SaveInternal(file);
      return true;
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
  }

 private:
  struct Header {
    int64_t header_size;
    int64_t g_size;
  };

  void Create() {
    Tick tick(__FN__);

    GenerateG1(0xffffffffffffffffULL, &h_);

    GenerateG1(0xfffffffffffffffeULL, &u_);

    auto parallel_f = [this](int64_t i) { GenerateG1(i, &g_[i]); };
    parallel::For(GSize(), parallel_f);
  }

  void SaveInternal(std::string const& file) {
    Tick tick(__FN__);
    FILE* f = fopen(file.c_str(), "wb+");
    CHECK(f, file);

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    header.header_size = sizeof(Header);
    header.g_size = GSize();

    CHECK(WriteHeader(f, header), "");
    CHECK(WriteG1(f, h_), "");
    CHECK(WriteG1(f, u_), "");
    for (auto i = 0; i < GSize(); ++i) {
      auto& g = g_[i];
      CHECK(WriteG1(f, g), "");
    }
  }

  void LoadInternal(std::string const& file) {
    Tick tick(__FN__);
    FILE* f = fopen(file.c_str(), "rb");
    CHECK(f, file);

    std::unique_ptr<FILE, decltype(&fclose)> auto_close(f, fclose);

    Header header;
    CHECK(ReadHeader(f, header), file);

    CHECK(header.g_size == GSize(), file);

    CHECK(ReadG1(f, h_), file);

    CHECK(ReadG1(f, u_), file);

    for (auto i = 0; i < GSize(); ++i) {
      auto& g = g_[i];
      CHECK(ReadG1(f, g), file);
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
  G1 u_;
  G1 h_;
  G1* g_;
};

inline Base& GetPcBase(std::string const& file = "") {
  static std::unique_ptr<Base> _instance_(new Base(file));
  return *_instance_;
}

inline bool operator==(Base const& a, Base const& b) {
  if (a.u() != b.u()) return false;
  if (a.h() != b.h()) return false;
  for (int64_t i = 0; i < Base::GSize(); ++i) {
    if (a.g()[i] != b.g()[i]) return false;
  }
  return true;
}

inline bool operator!=(Base const& a, Base const& b) { return !(a == b); }

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

    std::unique_ptr<Base> base(new Base);
    if (!base->Save(file)) {
      std::cerr << "Save ecc pds file" << file << " failed.\n";
      return false;
    }

    std::cout << "Create pc base file success.\n";
    if (!Load(file)) return false;
    DCHECK(*base == GetPcBase(), "");
    return true;
  } catch (std::exception& e) {
    std::cerr << "Create pc base file exception: " << e.what() << "\n";
    return false;
  }
}

inline G1 const* PcG() {
  auto const& base = GetPcBase();
  return base.g();
}

inline G1 const& PcG(int64_t i) {
  auto const& base = GetPcBase();
  if (i == -1) {
    return base.u();
  } else {
    CHECK(i < Base::GSize(), std::to_string(i));
    return base.g()[i];
  }
}

inline G1 const& PcU() {
  auto const& base = GetPcBase();
  return base.u();
}

inline G1 const& PcH() {
  auto const& base = GetPcBase();
  return base.h();
}

inline G1 const& PcHG(int64_t i) {
  if (i == 0) return PcH();
  return PcG(i - 1);
}

inline GetRefG1 const kGetRefG1 = [](int64_t i) -> G1 const& {
  CHECK(i < Base::GSize(), std::to_string(i));
  return pc::PcG()[i];
};

inline G1 ComputeCom(int64_t n, GetRefG1 const& g, G1 const& h,
                     GetRefFr const& x, Fr const& r, bool check_01 = false) {
  auto get_g = [&g, &h](int64_t i) -> G1 const& { return i ? g(i - 1) : h; };
  auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x(i - 1) : r; };
  return MultiExpBdlo12<G1>(get_g, get_f, n + 1, check_01);
}

inline G1 ComputeCom(int64_t n, GetRefG1 const& g, GetRefFr const& x,
                     Fr const& r, bool check_01 = false) {
  return ComputeCom(n, g, PcH(), x, r, check_01);
}

inline G1 ComputeCom(int64_t n, GetRefG1 const& g, Fr const* x, Fr const& r,
                     bool check_01 = false) {
  auto get_f = [&x](int64_t i) -> Fr const& { return x[i]; };
  return ComputeCom(n, g, get_f, r, check_01);
}

inline G1 ComputeCom(GetRefG1 const& g, std::vector<Fr> const& x, Fr const& r,
                     bool check_01 = false) {
  return ComputeCom(x.size(), g, x.data(), r, check_01);
}

inline G1 ComputeCom(int64_t n, G1 const* g, G1 const& h, Fr const* x,
                     Fr const& r, bool check_01 = false) {
  auto get_g = [&g, &h](int64_t i) -> G1 const& { return i ? g[i - 1] : h; };
  auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x[i - 1] : r; };
  return MultiExpBdlo12<G1>(get_g, get_f, n + 1, check_01);
}

inline G1 ComputeCom(int64_t n, G1 const* g, Fr const* x, Fr const& r,
                     bool check_01 = false) {
  return ComputeCom(n, g, PcH(), x, r, check_01);
}

inline G1 ComputeCom(int64_t n, Fr const* x, Fr const& r,
                     bool check_01 = false) {
  return ComputeCom(n, PcG(), PcH(), x, r, check_01);
}

inline G1 ComputeCom(std::vector<Fr> const& x, Fr const& r,
                     bool check_01 = false) {
  return ComputeCom(x.size(), x.data(), r, check_01);
}

inline G1 ComputeCom(int64_t n, GetRefFr const& x, Fr const& r,
                     bool check_01 = false) {
  auto get_g = [](int64_t i) -> G1 const& { return PcG(i); };
  return ComputeCom(n, get_g, PcH(), x, r, check_01);
}

inline G1 ComputeCom(G1 const& g, G1 const& h, Fr const& x, Fr const& r) {
  return g * x + h * r;
}

inline G1 ComputeCom(G1 const& g, Fr const& x, Fr const& r) {
  return g * x + PcH() * r;
}

inline G1 ComputeCom(Fr const& x, Fr const& r) {
  return PcG(0) * x + PcH() * r;
}

inline std::vector<G1> CopyG(GetRefG1 const& get_g, int64_t n) {
  std::vector<G1> ret(n);
  for (int64_t i = 0; i < n; ++i) {
    ret[i] = get_g(i);
  }
  return ret;
}

inline G1 ComputeSigmaG(GetRefG1 const& get_g, int64_t n) {
  G1 ret = G1Zero();
  for (int64_t i = 0; i < n; ++i) {
    ret += get_g(i);
  }
  return ret;
}

inline G1 ComputeSigmaG(size_t offset, int64_t n) {
  G1 ret = G1Zero();
  for (int64_t i = 0; i < n; ++i) {
    ret += PcG(i + offset);
  }
  return ret;
}

inline bool TestPcCommitment(int64_t n) {
  if (n > Base::GSize()) return false;
  std::vector<Fr> x(n);
  FrRand(x);
  Fr r = FrRand();
  auto& base = GetPcBase();
  G1 left, right;
  {
    Tick tick("multiexp1");
    left = ComputeCom(n, base.g(), base.h(), x.data(), r);
  }
  {
    Tick tick("multiexp2");
    right = MultiExp(base.g(), x.data(), n) + base.h() * r;
  }

  bool success = left == right;
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace pc
