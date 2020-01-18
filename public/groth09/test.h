#include "groth09.h"

namespace groth09 {
inline void Test() {
  std::vector<bool> ret;

  ret.push_back(hyrax::a2::TestRom(10));

  ret.push_back(sec51a::TestRom(17));
  ret.push_back(sec51b::TestRom(13));  

  ret.push_back(sec52a::TestRom(17, 20));
  ret.push_back(sec52b::TestRom(17, 20));

  ret.push_back(sec53a::TestRom(13, 10));
  ret.push_back(sec53b::TestRom(13, 10));

  ret.push_back(sec43b::TestRom(2, 14));

  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}

inline void TestPerformance() {
  std::vector<bool> ret;

  ret.push_back(groth09::sec51a::TestRom(123));
  ret.push_back(groth09::sec51b::TestRom(123));  

  ret.push_back(groth09::sec52a::TestRom(17, 20));
  ret.push_back(groth09::sec52b::TestRom(17, 20));

  ret.push_back(groth09::sec53a::TestRom(128, 1024 * 16));
  ret.push_back(groth09::sec53b::TestRom(128, 1024 * 16));

  ret.push_back(groth09::sec43b::TestRom(256, 1024 * 16));
  ret.push_back(groth09::sec43b::TestRom(1, 1024 * 16));

  ret.push_back(hyrax::a2::TestRom(1000));

  std::cout << __FUNCTION__ << " summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}
}  // namespace groth09