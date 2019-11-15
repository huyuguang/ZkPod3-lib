#pragma once

#include <condition_variable>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <queue>
#include <deque>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <assert.h>

namespace parallel {

namespace details {
inline thread_local bool in_pool_thread = false;
}  // namespace details

typedef std::function<void()> Task;

class SimpleTaskPool {
 public:
  SimpleTaskPool() {
    char const* thread_num_env = std::getenv("options:thread_num");  
    size_t thread_num =
        thread_num_env ? std::strtoul(thread_num_env, nullptr, 10) : 0;
    if (!thread_num) thread_num = std::thread::hardware_concurrency();   

    if (thread_num > 1) {
      threads_.resize(thread_num);
      for (auto& t : threads_) {
        t = std::thread([this]() { ThreadFunc(); });
      }
    }
  }

  ~SimpleTaskPool() {
    exit_ = true;
    cv_.notify_all();
    for (auto& t : threads_) {
      t.join();
    }
  }

  size_t thread_sum() const { return threads_.size(); }

  void Execute(std::vector<Task>& tasks) {
    if (threads_.empty()) {
      ExecuteSerial(tasks);
    } else if (details::in_pool_thread) {
      ExecuteSerial(tasks);
    } else {
      ExecuteParallel(tasks);
    }
  }

 private:
  void ThreadFunc() {
    details::in_pool_thread = true;
#ifdef _WIN32
    SetThreadDescription(GetCurrentThread(), L"task_pool");
#endif
    Task task;
    for (;;) {      
      if (PopTask(task)) {
        task();
      } else {
        if (exit_) break;
        WaitNewTask();
      }
    }
  }

  void WaitNewTask() {
    if (exit_) return;
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this]() { return !tasks_.empty() || exit_; });
  }

  bool PopTask(Task& task) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (tasks_.empty()) return false;
    task = std::move(tasks_.front());
    tasks_.pop();
    return true;
  }

  void PushTasks(std::vector<Task>&& tasks) {
    std::unique_lock<std::mutex> lock(mutex_);
    for (auto&& i : tasks) {
      tasks_.emplace(std::move(i));
    }
    cv_.notify_all();
  }

  void ExecuteParallel(std::vector<Task>& tasks) {
    assert(!details::in_pool_thread);

    struct Context {
      Context(size_t left_count) : left_count(left_count) {}
      std::atomic<size_t> left_count;
      std::mutex mutex;
      std::condition_variable cv;
    };
    auto context = std::make_shared<Context>(tasks.size());

    std::vector<Task> wrapped_tasks(tasks.size());
    for (size_t i = 0; i < tasks.size(); ++i) {
      auto& task = tasks[i];
      wrapped_tasks[i] = [&task, context]() {
        task();
        std::unique_lock<std::mutex> lock(context->mutex);
        if (0 == --context->left_count) {
          context->cv.notify_one(); // NOTE: here is context, not cv_ and mutex_
        }
      };
    }

    PushTasks(std::move(wrapped_tasks));

    std::unique_lock<std::mutex> lock(context->mutex);
    context->cv.wait(lock, [context]() { return !context->left_count; });
  }

  void ExecuteSerial(std::vector<Task>& tasks) {
    for (auto& t : tasks) {
      t();
    }
  }

 private:
  std::vector<std::thread> threads_;
  std::mutex mutex_;
  std::condition_variable cv_;
  std::queue<Task> tasks_;
  std::atomic<bool> exit_{false};
};

}  // namespace parallel