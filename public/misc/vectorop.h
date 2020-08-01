#pragma once

#include <functional>
#include <vector>

#include "parallel/parallel.h"

template <typename T>
void VectorMul(std::vector<T>& c, std::vector<T> const& a, T const& b) {
  c.resize(a.size());
  auto parallel_f = [&c, &a, &b](size_t i) { c[i] = a[i] * b; };
  parallel::For(a.size(), parallel_f, a.size() < 16 * 1024);
}

template <typename T>
void VectorMul(std::vector<T>& c, int64_t n,
               std::function<T const&(int64_t)>& get_a, T const& b) {
  c.resize(n);
  auto parallel_f = [&c, &get_a, &b](size_t i) { c[i] = get_a(i) * b; };
  parallel::For(n, parallel_f, n < 16 * 1024);
}

template <typename T>
void VectorAdd(std::vector<T>& c, std::vector<T> const& a, T const& b) {
  c.resize(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    c[i] = a[i] + b;
  }
}

template <typename T>
void VectorAdd(std::vector<T>& c, int64_t n,
               std::function<T const&(int64_t)>& get_a, T const& b) {
  c.resize(n);
  for (int64_t i = 0; i < n; ++i) {
    c[i] = get_a(i) + b;
  }
}

template <typename T>
void VectorAdd(std::vector<T>& c, std::vector<T> const& a,
               std::vector<T> const& b) {
  auto const& aa = a.size() >= b.size() ? a : b;
  auto const& bb = a.size() >= b.size() ? b : a;

  c.resize(aa.size());
  for (size_t i = 0; i < aa.size(); ++i) {
    if (i < bb.size()) {
      c[i] = a[i] + b[i];
    } else {
      c[i] = aa[i];
    }
  }
}

template <typename T>
void VectorAdd(std::vector<T>& c, int64_t n,
               std::function<T const&(int64_t)>& get_a,
               std::function<T const&(int64_t)>& get_b) {
  c.resize(n);
  for (int64_t i = 0; i < n; ++i) {
    c[i] = get_a(i) + get_b(i);
  }
}

template <typename T>
void VectorInc(std::vector<T>& a, std::vector<T> const& b) {
  VectorAdd(a, a, b);
}

template <typename T>
std::vector<T> operator*(std::vector<T> const& a, T const& b) {
  std::vector<T> c;
  VectorMul(c, a, b);
  return c;
}

template <typename T>
std::vector<T>& operator*=(std::vector<T>& a, T const& b) {
  VectorMul(a, a, b);
  return a;
}

template <typename T>
std::vector<T> operator+(std::vector<T> const& a, T const& b) {
  std::vector<T> c;
  VectorAdd(c, a, b);
  return c;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, T const& b) {
  VectorAdd(a, a, b);
  return a;
}

template <typename T>
std::vector<T> operator+(std::vector<T> const& a, std::vector<T> const& b) {
  std::vector<T> c;
  VectorAdd(c, a, b);
  return c;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, std::vector<T> const& b) {
  VectorAdd(a, a, b);
  return a;
}

template <typename T>
std::vector<T> operator-(std::vector<T> const& a) {
  std::vector<T> c(a.size());
  for (size_t i = 0; i < c.size(); ++i) {
    c[i] = -a[i];
  }
  return c;
}
