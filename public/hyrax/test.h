#include "hyrax.h"

namespace hyrax {
inline void Test() {
  std::vector<bool> ret;

  ret.push_back(a1::Test());
  ret.push_back(a2::Test(14));
  ret.push_back(a3::Test(13));
  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}

inline void TestPerformance() {
  std::vector<bool> ret;

  ret.push_back(a1::Test());
  ret.push_back(a2::Test(1024));
  ret.push_back(a3::Test(1024));
  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}

}  // namespace hyrax