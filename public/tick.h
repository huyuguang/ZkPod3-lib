#pragma once

#include <chrono>
#include <iostream>
#include <string>

struct Tick {
  Tick(std::string const& desc) : desc_(desc) {
    start_ = std::chrono::steady_clock::now();
    std::cout << GetIndentString(UpdateIndent(1) - 1);
    std::cout << "==> " << desc_ << "\n";
  }
  Tick(std::string const& desc, std::string const& desc2)
      : desc_(desc + " " + desc2) {
    start_ = std::chrono::steady_clock::now();
    std::cout << GetIndentString(UpdateIndent(1) - 1);
    std::cout << "==> " << desc_ << "\n";
  }
  ~Tick() {
    auto t = std::chrono::steady_clock::now() - start_;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(t);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t);
    std::cout << GetIndentString(UpdateIndent(-1));
    if (s.count() < 10) {
      std::cout << "<== " << desc_ << " tick: " << ms.count() << " ms\n";
    } else {
      std::cout << "<== " << desc_ << " tick: " << s.count() << " seconds\n";
    }
  }

  static int UpdateIndent(int inc) {
    static thread_local int indent = 0;
    indent += inc;
    assert(indent >= 0);
    return indent;
  }

  static std::string GetIndentString(int indent) {
    return std::string(2 * indent, ' ');
  }

  static std::string GetIndentString() {
    return std::string(2 * UpdateIndent(0), ' ');
  }

 private:
  std::string desc_;
  std::chrono::steady_clock::time_point start_;
};