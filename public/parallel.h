#pragma once

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include "tick.h"

#ifdef USE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
// #include <tbb/tbb.h>
namespace parallel {

typedef std::function<void()> Task;

template <typename T, typename F>
void For(T count, F& f, bool direct = false) {
  if (!count) return;
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

template <typename T, typename F>
void For(T begin, T end, F& f, bool direct = false) {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, direct);
}

inline void Invoke(std::vector<Task>& tasks, bool direct = false) {
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

}  // namespace parallel

#else
#include "recursive_taskpool.h"
namespace parallel {

namespace details {

inline SimpleTaskPool& GetTaskPool() {
  static SimpleTaskPool instance;
  return instance;
}

}  // namespace details

typedef std::function<void()> Task;

template <typename T, typename F>
void For(T count, F& f, bool direct = false) {
  if (!count) return;
  if (direct) {
    for (T i = 0; i < count; ++i) f(i);
    return;
  }

  auto& task_pool = details::GetTaskPool();
  size_t thread_sum = task_pool.thread_sum();
  if (!thread_sum) {
    for (T i = 0; i < count; ++i) {
      f(i);
    }
    return;
  }

  std::vector<Task> tasks;
  if ((size_t)count > thread_sum) {
    size_t bundle_count = (size_t)count / thread_sum;
    size_t left_count = count % thread_sum;
    tasks.reserve(thread_sum + left_count);
    for (size_t i = 0; i < thread_sum; ++i) {
      size_t begin = i * bundle_count;
      size_t end = begin + bundle_count;
      tasks.emplace_back([&f, begin, end]() mutable {
        for (size_t i = begin; i < end; ++i) f(i);
      });
    }
    for (size_t i = 0; i < left_count; ++i) {
      size_t offset = thread_sum * bundle_count + i;
      tasks.emplace_back([&f, offset]() mutable { f(offset); });
    }
  } else {
    tasks.reserve(count);
    for (T i = 0; i < count; ++i) {
      tasks.emplace_back([&f, i]() { f(i); });
    }
  }
  task_pool.Execute(tasks);
}

template <typename T, typename F>
void For(T begin, T end, F& f, bool direct = false) {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, direct);
}

inline void Invoke(std::vector<Task>& tasks, bool direct = false) {
  if (tasks.empty()) return;

  if (tasks.size() == 1) return tasks[0]();

  if (direct) {
    for (auto& task : tasks) {
      task();
    }
    return;
  }

  details::GetTaskPool().Execute(tasks);
}

template <class InputIt, class T>
T Accumulate(InputIt first, InputIt last, T init) {
  return std::accumulate(std::forward<InputIt>(first),
                         std::forward<InputIt>(last), std::forward<T>(init));
}

template <class InputIt, class T, class BinaryOperation>
T Accumulate(InputIt first, InputIt last, T init, BinaryOperation op) {
  return std::accumulate(std::forward<InputIt>(first),
                         std::forward<InputIt>(last), std::forward<T>(init),
                         std::forward<BinaryOperation>(op));
}

}  // namespace parallel
#endif
