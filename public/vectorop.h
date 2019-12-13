#pragma once

#include <vector>

#include "ecc.h"

template <typename T>
void VectorMul(std::vector<T>& c, std::vector<T> const& a, T const& b) {
  c.resize(a.size());
  auto parallel_f = [&c, &a, &b](size_t i) { c[i] = a[i] * b; };
  parallel::For(a.size(), parallel_f, a.size() < 16 * 1024);
}

template <typename T>
void VectorAdd(std::vector<T>& c, std::vector<T> const& a, T const& b) {
  c.resize(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    c[i] = a[i] + b;
  }
}

template <typename T>
void VectorAdd(std::vector<T>& c, std::vector<T> const& a,
               std::vector<T> const& b) {
  assert(a.size() == b.size());
  c.resize(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    c[i] = a[i] + b[i];
  }
}

template <typename T>
void VectorInc(std::vector<T>& a, std::vector<T> const& b) {
  assert(a.size() == b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] += b[i];
  }
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