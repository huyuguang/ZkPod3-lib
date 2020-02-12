#pragma once

#include "./multiexp.h"
#include "parallel/parallel.h"

template <typename G, typename GET_G, typename GET_F>
G ParallelMultiExpBdlo12(GET_G const& get_g, GET_F const& get_f, size_t n,
                         bool check_01 = false) {
  auto thread_num = parallel::tbb_thread_num;
  if (thread_num <= 1 || (int64_t)n < thread_num) {
    return MultiExpBdlo12<G, GET_G, GET_F>(get_g, get_f, n, check_01);
  }

  struct Item {
    size_t begin;
    size_t end;
    G1 ret;
  };

  auto size = n / thread_num;
  std::vector<Item> items(thread_num);
  for (int i = 0; i < thread_num; ++i) {
    auto& item = items[i];
    item.begin = i * size;
    item.end = item.begin + size;
  }
  items.back().end = n;

  auto f = [check_01, &items, &get_g, &get_f](int64_t i) {
    auto& item = items[i];
    auto item_get_g = [&get_g, &item](int64_t j) -> decltype(get_g(0)) {
      return get_g(j + item.begin);
    };
    auto item_get_f = [&get_f, &item](int64_t j) -> decltype(get_f(0)) {
      return get_f(j + item.begin);
    };

    item.ret = MultiExpBdlo12<G>(item_get_g, item_get_f, item.end - item.begin,
                                 check_01);
  };
  parallel::For(items.size(), f);

  G ret;
  ret.clear();

  return parallel::Accumulate(
      items.begin(), items.end(), ret,
      [](G const& a, Item const& b) { return a + b.ret; });
}

template <typename G, typename GET_G>
G ParallelMultiExpBdlo12(GET_G const& get_g, std::vector<Fr> const& f, size_t n,
                         bool check_01 = false) {
  auto get_f = [&f](int64_t i) -> Fr const& { return f[i]; };
  return ParallelMultiExpBdlo12<G>(get_g, get_f, n, check_01);
}
