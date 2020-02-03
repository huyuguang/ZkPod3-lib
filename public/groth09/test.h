#include "groth09.h"

namespace groth09 {
inline void Test() {
  std::vector<bool> ret;

  ret.push_back(hyrax::a2::Test(10));

  ret.push_back(sec51a::Test(17));
  ret.push_back(sec51b::Test(13));  

  ret.push_back(sec52a::Test(17, 20));
  ret.push_back(sec52b::Test(17, 20));

  ret.push_back(sec53a::Test(13, 10));
  ret.push_back(sec53b::Test(13, 10));

  ret.push_back(sec43b::Test(2, 14));

  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}

inline void TestPerformance() {
  std::vector<bool> ret;

  ret.push_back(groth09::sec51a::Test(123));
  ret.push_back(groth09::sec51b::Test(123));  

  ret.push_back(groth09::sec52a::Test(17, 20));
  ret.push_back(groth09::sec52b::Test(17, 20));

  ret.push_back(groth09::sec53a::Test(128, 1024 * 16));
  ret.push_back(groth09::sec53b::Test(128, 1024 * 16));

  ret.push_back(groth09::sec43b::Test(256, 1024 * 16));
  ret.push_back(groth09::sec43b::Test(1, 1024 * 16));

  ret.push_back(hyrax::a2::Test(1000));

  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}
}  // namespace groth09