#pragma once

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "tick.h"
#include "recursive_taskpool.h"

#if 0
#if defined(_MSC_VER)
#pragma warning(push, 0)
#endif
#include <boost/asio.hpp>
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

namespace parallel {

namespace details {
inline thread_local bool running_in_this_thread = false;

inline boost::asio::thread_pool basic_thread_pool;
}  // namespace details

template <typename T, typename F>
void For(T count, F& f, std::string const& desc = "") {
  using details::basic_thread_pool;
  if (!count) return;
  Tick _tick_(__FUNCTION__, desc + " " + std::to_string(count));
  if (basic_thread_pool.get_executor().running_in_this_thread()) {
    for (T i = 0; i < count; ++i) {
      f(i);
    }
  } else {
    std::mutex mutex;
    std::condition_variable cv;
    auto left_count = count;

    auto notify = [&left_count, &mutex, &cv]() {
      std::unique_lock<std::mutex> lock(mutex);
      if (0 == --left_count) {
        cv.notify_one();
      }
    };

    auto wait = [&left_count, &mutex, &cv]() {
      std::unique_lock<std::mutex> lock(mutex);
      while (left_count) cv.wait(lock);
    };

    for (T i = 0; i < count; ++i) {
      boost::asio::post(basic_thread_pool,
                        [i, &left_count, &f, &notify]() mutable {
                          f(i);
                          notify();
                        });
    }
    wait();
  }
}

template <typename T, typename F>
void For(T begin, T end, F& f, std::string const& desc = "") {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, desc);
}

}  // namespace parallel
#endif

#if 0

#include <boost/thread/executors/basic_thread_pool.hpp>

namespace parallel {

namespace details {
inline thread_local bool running_in_this_thread = false;

inline void AtThreadEntry(boost::basic_thread_pool& /*pool*/) {
  running_in_this_thread = true;
}

inline boost::basic_thread_pool basic_thread_pool(
    boost::thread::hardware_concurrency(), &AtThreadEntry);
}  // namespace details

template <typename T, typename F>
void For(T count, F& f, std::string const& desc = "") {
  using details::basic_thread_pool;
  using details::running_in_this_thread;
  if (!count) return;
  std::cout << "parallel::For " << desc << " " << count << "\n";

  std::mutex mutex;
  std::condition_variable cv;
  auto left_count = count;
  auto notify = [&left_count, &mutex, &cv]() {
    std::unique_lock<std::mutex> lock(mutex);
    if (0 == --left_count) {
      cv.notify_one();
    }
  };

  auto wait = [&left_count, &mutex, &cv]() {
    std::unique_lock<std::mutex> lock(mutex);
    while (left_count) cv.wait(lock);
  };

  auto is_ready = [&left_count, &mutex]() {
    std::unique_lock<std::mutex> lock(mutex);
    return left_count == 0;
  };

  for (T i = 0; i < count; ++i) {
    basic_thread_pool.submit([i, &f, &notify]() mutable {
      f(i);
      notify();
    });
  }

  if (running_in_this_thread) {
    if (!basic_thread_pool.reschedule_until(std::move(is_ready))) wait();
  } else {
    wait();
  }
}

template <typename T, typename F>
void For(T begin, T end, F& f, std::string const& desc = "") {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, desc);
}

inline void Exec(std::vector<std::function<void()>>& executors,
                 std::string const& desc = "") {
  using details::basic_thread_pool;
  using details::running_in_this_thread;
  if (executors.empty()) return;
  std::cout << "parallel::Exec " << desc << " " << executors.size() << "\n";

  std::mutex mutex;
  std::condition_variable cv;
  auto left_count = executors.size();
  
  auto notify = [&left_count, &mutex, &cv]() {
    std::unique_lock<std::mutex> lock(mutex);
    if (0 == --left_count) {
      cv.notify_one();
    }
  };

  auto wait = [&left_count, &mutex, &cv]() {
    std::unique_lock<std::mutex> lock(mutex);
    while (left_count) cv.wait(lock);
  };

  auto is_ready = [&left_count, &mutex]() {
    std::unique_lock<std::mutex> lock(mutex);
    return left_count == 0;
  };

  for (auto& executor : executors) {
    basic_thread_pool.submit([&executor, &notify]() mutable {
      executor();
      notify();
    });
  }

  if (running_in_this_thread) {
    if (!basic_thread_pool.reschedule_until(std::move(is_ready))) wait();
  } else {
    wait();
  }
}
}  // namespace parallel
#endif


namespace parallel {

namespace details {

inline RecursiveTaskPool& GetTaskPool() {
  static RecursiveTaskPool recursive_task_pool;
  return recursive_task_pool;
}

}  // namespace details

template <typename T, typename F>
void For(T count, F& f, std::string const& desc = "") {
  if (!count) return;
  (void)desc;

  // Tick _tick_(__FUNCTION__, desc + " " + std::to_string(count));
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
  task_pool.PostAndWait(tasks);
}

template <typename T, typename F>
void For(T begin, T end, F& f, std::string const& desc = "") {
  auto count = end - begin;
  auto f2 = [begin, &f](T i) { return f(begin + i); };
  For(count, f2, desc);
}

inline void SyncExec(std::vector<Task>& tasks) {
  auto& task_pool = details::GetTaskPool();
  static const size_t thread_sum = task_pool.thread_sum();
  if (!thread_sum) {
    for (auto& task : tasks) {
      task();
    }
    return;
  }

  task_pool.PostAndWait(tasks);
}
}  // namespace parallel