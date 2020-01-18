#include "cmd.h"

namespace cmd {
inline void Test() {
  std::vector<bool> ret;

  ret.push_back(match_query::Test());
  ret.push_back(substr_query::Test());

  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}
}