#pragma once

#include <utility>

#include "./funcs.h"
#include "./types.h"
#include "log/tick.h"

extern bool DISABLE_TBB;

template <typename G, typename GET_G, typename GET_F>
G MultiExpBdlo12Inner(GET_G const& get_g, GET_F const& get_f, size_t n) {
  G result = GZero<G>();
  if (n == 0) return result;
  if (n == 1) return get_g(0) * get_f(0);

  //if (n < 32) {
  //  for (size_t i = 0; i < n; ++i) {
  //    result += get_g(i) * get_f(i);
  //  }
  //  return result;
  //}

  auto local_log2 = [](size_t n) {
    // returns ceil(log2(n)), so ((size_t)1)<<log2(n) is the smallest power of
    // 2, that is not less than n.
    size_t r = ((n & (n - 1)) == 0 ? 0 : 1);  // add 1 if n is not power of 2
    while (n > 1) {
      n >>= 1;
      r++;
    }
    return r;
  };

  size_t log2_length = local_log2(n);
  size_t c = log2_length - (log2_length / 3 - 2);
  std::vector<mpz_class> bn_exponents(n);
  size_t num_bits = 0;

  for (size_t i = 0; i < n; i++) {
    bn_exponents[i] = get_f(i).getMpz();
    auto bits = mpz_sizeinbase(bn_exponents[i].get_mpz_t(), 2);
    num_bits = std::max(num_bits, bits);
  }

  size_t num_groups = (num_bits + c - 1) / c;

  bool result_nonzero = false;

  for (size_t k = num_groups - 1; k <= num_groups; k--) {
    if (result_nonzero) {
      for (size_t i = 0; i < c; i++) {
        G::dbl(result, result);
      }
    }

    std::vector<G> buckets((size_t)1 << c);
    std::vector<bool> bucket_nonzero((size_t)1 << c);

    for (size_t i = 0; i < n; i++) {
      size_t id = 0;
      for (size_t j = 0; j < c; j++) {
        if (mpz_tstbit(bn_exponents[i].get_mpz_t(), k * c + j)) {
          id |= (size_t)1 << j;
        }
      }

      if (id == 0) {
        continue;
      }

      if (bucket_nonzero[id]) {
        buckets[id] = buckets[id] + get_g(i);
      } else {
        buckets[id] = get_g(i);
        bucket_nonzero[id] = true;
      }
    }

    G running_sum;
    running_sum.clear();
    bool running_sum_nonzero = false;

    for (size_t i = ((size_t)1 << c) - 1; i > 0; i--) {
      if (bucket_nonzero[i]) {
        if (running_sum_nonzero) {
          running_sum = running_sum + buckets[i];
        } else {
          running_sum = buckets[i];
          running_sum_nonzero = true;
        }
      }

      if (running_sum_nonzero) {
        if (result_nonzero) {
          result = result + running_sum;
        } else {
          result = running_sum;
          result_nonzero = true;
        }
      }
    }
  }

  return result;
}

template <typename G, typename GET_G, typename GET_F>
G MultiExpBdlo12Inner(GET_G const& get_g, GET_F const& get_f, size_t n,
                      bool check_01) {
  if (!check_01) {
    return MultiExpBdlo12Inner<G>(get_g, get_f, n);
  } else {
    int count_0 = 0;
    int count_1 = 0;
    G ret = GZero<G>();
    std::vector<size_t> pos;
    pos.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      if (get_f(i).isZero()) {
        ++count_0;
        continue;
      } else if (get_f(i).isOne()) {
        ret += get_g(i);
        ++count_1;
      } else {
        pos.push_back(i);
      }
    }

    if (pos.size() == n) {
      return MultiExpBdlo12Inner<G, GET_G, GET_F>(get_g, get_f, n);
    }
    // std::cout << "MultiExpBdlo12: " << n << ", "
    //          << "0: " << count_0 << ", 1: " << count_1
    //          << ", left: " << pos.size() << "\n";

    if (!pos.empty()) {
      auto new_get_g = [&get_g, &pos](int64_t i) { return get_g(pos[i]); };
      auto new_get_f = [&get_f, &pos](int64_t i) { return get_f(pos[i]); };
      ret += MultiExpBdlo12Inner<G>(new_get_g, new_get_f, pos.size());
    }
    return ret;
  }
}

template <typename G, typename GET_G, typename GET_F>
G ParallelMultiExpBdlo12Inner(GET_G const& get_g, GET_F const& get_f, size_t n,
                              bool check_01 = false) {
  auto thread_num = parallel::tbb_thread_num;
  if (DISABLE_TBB || thread_num <= 1 || (int64_t)n < thread_num) {
    return MultiExpBdlo12Inner<G, GET_G, GET_F>(get_g, get_f, n, check_01);
  }

  struct Item {
    size_t begin;
    size_t end;
    G1 ret;
  };

  if (thread_num > 4) thread_num = 4;

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

    item.ret = MultiExpBdlo12Inner<G>(item_get_g, item_get_f,
                                      item.end - item.begin, check_01);
  };
  parallel::For(items.size(), f);

  G ret;
  ret.clear();

  return parallel::Accumulate(
      items.begin(), items.end(), ret,
      [](G const& a, Item const& b) { return a + b.ret; });
}

template <typename G, typename GET_G, typename GET_F>
G MultiExpBdlo12(GET_G const& get_g, GET_F const& get_f, size_t n,
                 bool check_01 = false) {
  if (n > 2000000 && !DISABLE_TBB)
    return ParallelMultiExpBdlo12Inner<G>(get_g, get_f, n, check_01);

  return MultiExpBdlo12Inner<G>(get_g, get_f, n, check_01);
}

template <typename G, typename GET_G>
G MultiExpBdlo12(GET_G const& get_g, std::vector<Fr> const& f, size_t n,
                 bool check_01 = false) {
  auto get_f = [&f](int64_t i) -> Fr const& { return f[i]; };
  return MultiExpBdlo12<G>(get_g, get_f, n, check_01);
}

template <typename G, typename GET_G>
G MultiExpBdlo12(GET_G const& get_g, Fr const* f, size_t n,
                 bool check_01 = false) {
  auto get_f = [f](int64_t i) -> Fr const& { return f[i]; };
  return MultiExpBdlo12<G>(get_g, get_f, n, check_01);
}

template <typename G, typename GET_F>
G MultiExpBdlo12(std::vector<G> const& g, GET_F const& get_f, size_t n,
                 bool check_01 = false) {
  auto get_g = [&g](int64_t i) -> Fr const& { return g[i]; };
  return MultiExpBdlo12<G>(get_g, get_f, n, check_01);
}

template <typename G>
G MultiExpBdlo12(G const* pg, Fr const* pf, size_t n, bool check_01 = false) {
  auto get_g = [pg](size_t i) -> G const& { return pg[i]; };
  auto get_f = [pf](size_t i) -> Fr const& { return pf[i]; };
  return MultiExpBdlo12<G>(get_g, get_f, n, check_01);
}

template <typename G>
G MultiExpBdlo12(std::vector<G> const& g, std::vector<Fr> const& f,
                 bool check_01 = false) {
  auto size = std::min(g.size(), f.size());
  return MultiExpBdlo12<G>(g.data(), f.data(), size, check_01);
}

template <typename G>
G MultiExpBdlo12(std::vector<G const*> const& g,
                 std::vector<Fr const*> const& f, bool check_01 = false) {
  auto size = std::min(g.size(), f.size());
  auto get_g = [&g](size_t i) -> G const& { return *g[i]; };
  auto get_f = [&f](size_t i) -> Fr const& { return *f[i]; };
  return MultiExpBdlo12<G>(get_g, get_f, size, check_01);
}

template <typename G>
G1 MultiExpBdlo12(std::vector<G const*> const& g, std::vector<Fr> const& f,
                  uint64_t offset, uint64_t count, bool check_01 = false) {
  auto get_g = [&g, offset](size_t i) -> G const& { return *g[i + offset]; };
  auto get_f = [&f, offset](size_t i) -> Fr const& { return f[i + offset]; };
  return MultiExpBdlo12<G>(get_g, get_f, count, check_01);
}

inline bool TestEccDbl(int64_t n) {
  Tick tick(__FN__);
  G1 g = G1Rand();
  for (int64_t i = 0; i < n; ++i) {
    G1::dbl(g, g);
  }
  return g != G1Zero();
}

inline bool TestEccAdd(int64_t n) {
  Tick tick(__FN__);
  G1 g = G1Rand();
  G1 h = G1Rand();
  for (int64_t i = 0; i < n; ++i) {
    g = g + h;
  }
  return g != G1Zero();
}

inline bool TestMultiexp(int64_t n) {
  std::vector<Fr> f(n);
  FrRand(f);
  std::vector<G1> g(n);
  G1Rand(g);
  auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < 10; ++i) {
    Tick tick(__FN__, std::to_string(i));
    MultiExpBdlo12<G1>(g, f);
  }
  auto end = std::chrono::steady_clock::now();
  auto t = (end - start)/10;
  auto s = std::chrono::duration_cast<std::chrono::seconds>(t);
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t);
  std::cout << Tick::GetIndentString();
  if (s.count() < 100) {
    std::cout << "average: " << ms.count() << " ms\n";
  } else {
    std::cout << "average: " << s.count() << " seconds\n";
  }
  return true;
}