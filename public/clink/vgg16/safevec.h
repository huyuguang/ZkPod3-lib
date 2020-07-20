#pragma once

#include <thread>
#include <vector>

namespace clink::vgg16{

template<typename T>
class SafeVec {
 public:
  void emplace(T&& item) {
    std::lock_guard<std::mutex> lock(mutex_);
    items_.emplace_back(std::move(item));
  }
  void take(std::vector<T>& items) {
    std::lock_guard<std::mutex> lock(mutex_);
    items = std::move(items_);
  }
 private:
  std::vector<T> items_;
  std::mutex mutex_;
};

}