#pragma once

#include "./multiexp.h"
#include "parallel/parallel.h"

template <typename G, typename GET_G>
G ParallelMultiExpBdlo12(GET_G const& get_g, std::vector<Fr> const& f, size_t n,
                         bool check_01 = false) {
  auto get_f = [&f](int64_t i) -> Fr const& { return f[i]; };
  return ParallelMultiExpBdlo12Inner<G>(get_g, get_f, n, check_01);
}
