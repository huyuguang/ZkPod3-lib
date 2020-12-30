#pragma once

#include <assert.h>

#include <chrono>
#include <iostream>
#include <string>

struct Tick {
  Tick(std::string const& desc) : desc_(desc) {
    Uniform();
    start_ = std::chrono::steady_clock::now();
    std::cout << GetIndentString(GetIndent());
    std::cout << "==> " << desc_ << "\n";
    IncIndent();
  }
  Tick(std::string const& desc, std::string const& desc2) : desc_(desc) {
    Uniform();
    desc_ += " ";
    desc_ += desc2;
    start_ = std::chrono::steady_clock::now();
    std::cout << GetIndentString(GetIndent());
    std::cout << "==> " << desc_ << "\n";
    IncIndent();
  }
  ~Tick() {
    auto t = std::chrono::steady_clock::now() - start_;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(t);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t);
    DecIndent();
    std::cout << GetIndentString(GetIndent());
    if (s.count() < 100) {
      std::cout << "<== " << desc_ << " tick: " << ms.count() << " ms\n";
    } else {
      std::cout << "<== " << desc_ << " tick: " << s.count() << " seconds\n";
    }
  }

  static void IncIndent() { IndentInner(1, false); }

  static void DecIndent() { IndentInner(-1, false); }

  static void SetIndent(int v) { IndentInner(v, true); }

  static int GetIndent() { return IndentInner(0, false); }

  static std::string GetIndentString(int indent) {
    return std::string(2 * indent, ' ');
  }

  static std::string GetIndentString() { return GetIndentString(GetIndent()); }

 private:
  // return original indent
  static int IndentInner(int value, bool set) {
    static thread_local int indent = 0;
    int bak = indent;
    if (set) {
      indent = value;
    } else {
      indent += value;
    }
    assert(indent >= 0);
    return bak;
  }

  void Uniform() {
#ifdef __GNUC__
    auto pos = desc_.find_first_of('(');
    if (pos == std::string::npos) return;
    desc_.resize(pos);
#endif
  }

  std::string desc_;
  std::chrono::steady_clock::time_point start_;
};

#ifdef __GNUC__
#define __FN__ __PRETTY_FUNCTION__
#else
#define __FN__ __FUNCTION__
#endif

struct AutoTickIndent {
  AutoTickIndent(int new_indent) {
    old_indent = Tick::GetIndent();
    Tick::SetIndent(new_indent);
  }
  ~AutoTickIndent() { Tick::SetIndent(old_indent); }

 private:
  int old_indent;
};