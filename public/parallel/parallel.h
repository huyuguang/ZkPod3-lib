#pragma once

#include <tbb/scalable_allocator.h>
#include <tbb/tbb.h>
#include <tbb/tbb_allocator.h>
#include <tbb/tbbmalloc_proxy.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>
namespace parallel {

inline void CheckAllocationHook() {
#ifdef _WIN32
#ifndef _DEBUG
  char** func_replacement_log;
  int func_replacement_status =
      TBB_malloc_replacement_log(&func_replacement_log);

  if (func_replacement_status != 0) {
    std::cout << "tbbmalloc_proxy cannot replace memory allocation routines\n";
    for (char** log_string = func_replacement_log; *log_string != 0;
         log_string++) {
      std::cout << *log_string << "\n";
    }
    abort();
  }
#endif
#endif
}

typedef std::function<void()> Task;

template <typename T, typename F>
void For(T count, F& f, bool direct = false) {
  if (!count) return;
  if (count == 1) {
    f(0);
    return;
  }

  if (direct) {
    for (T i = 0; i < count; ++i) f(i);
    return;
  }

  auto f2 = [&f](const tbb::blocked_range<T>& range) {
    for (T i = range.begin(); i != range.end(); ++i) {
      f(i);
    }
  };
  tbb::parallel_for(tbb::blocked_range<T>(0, count), f2);
}

template <typename TaskContainer>
void Invoke(TaskContainer& tasks, bool direct = false) {
  if (tasks.empty()) return;

  if (tasks.size() == 1) return tasks[0]();

  if (direct) {
    for (auto& task : tasks) {
      task();
    }
    return;
  }

  if (tasks.size() == 2) {
    tbb::parallel_invoke(tasks[0], tasks[1]);
  } else if (tasks.size() == 3) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2]);
  } else if (tasks.size() == 4) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3]);
  } else if (tasks.size() == 5) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3], tasks[4]);
  } else if (tasks.size() == 6) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3], tasks[4],
                         tasks[5]);
  } else if (tasks.size() == 7) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3], tasks[4],
                         tasks[5], tasks[6]);
  } else if (tasks.size() == 8) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3], tasks[4],
                         tasks[5], tasks[6], tasks[7]);
  } else if (tasks.size() == 9) {
    tbb::parallel_invoke(tasks[0], tasks[1], tasks[2], tasks[3], tasks[4],
                         tasks[5], tasks[6], tasks[7], tasks[8]);
  } else {
    throw std::runtime_error("");
  }
}

template <class InputIt, class T>
T Accumulate(InputIt first, InputIt last, T init) {
  auto count = std::distance(first, last);
  if (count < 16 * 1024) {
    return std::accumulate(first, last, init);
  }

  return tbb::parallel_reduce(
      tbb::blocked_range<InputIt>(first, last), init,
      [](tbb::blocked_range<InputIt> const& range, T init) {
        return std::accumulate(range.begin(), range.end(), init);
      },
      std::plus<T>());
}

template <class InputIt, class T, class BinaryOperation>
T Accumulate(InputIt first, InputIt last, T init, BinaryOperation op) {
  auto count = std::distance(first, last);
  if (count < 16 * 1024) {
    return std::accumulate(first, last, init, std::move(op));
  }

  return tbb::parallel_reduce(
      tbb::blocked_range<InputIt>(first, last), init,
      [&op](tbb::blocked_range<InputIt> const& range, T init) {
        return std::accumulate(range.begin(), range.end(), init, op);
      },
      std::plus<T>());
}

template <typename T, typename F>
void For(bool* all_success, T count, F& f, bool direct = false) {
  std::vector<int64_t> rets(count);
  auto f2 = [&f, &rets](T i) { rets[i] = f(i); };
  For(count, f2, direct);
  auto is_true = [](int64_t const& r) { return !!r; };
  *all_success = std::all_of(rets.begin(), rets.end(), is_true);
  assert(*all_success);
}

template <typename T, typename F>
void For(T begin, T end, F& f, bool direct = false) {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, direct);
}

template <typename T, typename F>
void For(bool* all_success, T begin, T end, F& f, bool direct = false) {
  auto count = end - begin;
  auto f2 = [all_success, begin, &f](T i) { return f(begin + i); };
  For(all_success, count, f2, direct);
}
}  // namespace parallel
