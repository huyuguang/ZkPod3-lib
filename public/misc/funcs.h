#pragma once

#include "./types.h"
#include "public.h"

namespace misc {

inline std::string HexToStr(void const* p, size_t len) {
  std::string ret;
  ret.reserve(len * 4);
  char buf[4];
  uint8_t const* q = (uint8_t const*)p;
  for (size_t i = 0; i < len; ++i) {
    sprintf(buf, "%02x", q[i]);
    ret.append(buf);
  }
  return ret;
}

inline std::string HexToStr(h256_t const& h) {
  return HexToStr(h.data(), h.size());
}

inline bool StrToHex(char const* p, size_t len, uint8_t* q) {
  if (len % 2) return false;
  for (size_t i = 0; i < len / 2; ++i) {
    char buf[3];
    buf[0] = p[0];
    buf[1] = p[1];
    buf[2] = 0;
    int v = strtol(buf, NULL, 16);
    if (v > 256 || v < 0) return false;
    *q = (uint8_t)(v & 0xff);
    ++p;
    ++p;
    ++q;
  }
  return true;
}

template <size_t T>
std::array<uint8_t, T> StrToH(std::string const& s) {
  CHECK(s.size() >= T * 2, "");
  std::array<uint8_t, T> ret;
  StrToHex(s.data(), T * 2, ret.data());
  return ret;
}

inline bool StartWith(char const* p, char const* q) {
  for (;;) {
    char a = *p;
    char b = *q;
    if (!b) return true;
    if (!a) return false;
    if (a != b) return false;
    ++p;
    ++q;
  }
}

inline void HexStrToH256(std::string const& str, h256_t& h) {
  CHECK(str.size() == 32 * 2, "");
  StrToHex(str.c_str(), str.size(), h.data());
}

inline uint64_t Log2UB(uint64_t n) {
  assert(n);
  if (n == 1) return 0;
  if (n % 2) ++n;
  return 1 + Log2UB(n / 2);
}

inline uint64_t Pow2UB(uint64_t v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;
  return v;
}

inline bool Str2UInt(std::string s, uint64_t* v) {
  try {
    *v = boost::lexical_cast<uint64_t>(s);
    return true;
  } catch (boost::exception&) {
    return false;
  }
}

inline bool LoadTinyFile(std::string const& filename, std::string* text) {
  typedef std::unique_ptr<FILE, decltype(&fclose)> FileUniquePtr;

  struct stat st;
  if (::stat(filename.c_str(), &st)) return false;
  std::size_t size = st.st_size;
  if (size > 1024 * 1024 * 100) return false;

  text->resize(size);
  if (size == 0) return true;

  FileUniquePtr file(std::fopen(filename.c_str(), "rb"), fclose);
  if (!file) return false;

  return fread(&(*text)[0], 1, size, file.get()) == size;
}

inline bool SaveTinyFile(std::string const& filename, void const* data,
                         size_t size) {
  typedef std::unique_ptr<FILE, decltype(&fclose)> FileUniquePtr;
  FileUniquePtr file(std::fopen(filename.c_str(), "wb"), fclose);
  if (!file) return false;
  return fwrite(data, 1, size, file.get()) == size;
}

inline bool IsSameFile(std::string const& file1, std::string const& file2) {
  try {
    io::mapped_file_params params1;
    params1.path = file1;
    params1.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view1(params1);

    io::mapped_file_params params2;
    params2.path = file2;
    params2.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view2(params2);

    if (view1.size() != view2.size()) return false;

    return memcmp(view1.data(), view2.data(), view1.size()) == 0;
  } catch (std::exception&) {
    return false;
  }
}

inline bool GetFileSha256(std::string const& file, h256_t& h) {
  try {
    io::mapped_file_params params;
    params.path = file;
    params.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view(params);

    CryptoPP::SHA256 hash;
    hash.Update((uint8_t const*)view.data(), view.size());
    hash.Final(h.data());
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline std::string RandString(size_t max_len) {
  size_t len = rand() % max_len;
  if (!len) return std::string();
  std::string ret;
  ret.reserve(len);
  for (size_t i = 0; i < len; ++i) {
    ret.push_back((char)rand());
  }
  return ret;
}

inline void PrintVector(std::vector<Fr> const& a) {
  std::cout << "\n";
  for (auto const& i : a) {
    std::cout << i << "\n";
  }
  std::cout << "\n";
}

template <size_t N>
void PrintArray(std::array<Fr, N> const& a) {
  std::cout << "\n";
  for (auto const& i : a) {
    std::cout << i << "\n";
  }
  std::cout << "\n";
}

}  // namespace misc
