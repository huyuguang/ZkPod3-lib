#include "groth09.h"

namespace groth09 {
inline void Test() {
  std::vector<bool> ret;

#ifdef _DEBUG
  //ret.push_back(sec52::TestRom(17, 20));
  ret.push_back(sec53::TestRom(13, 104));
  ret.push_back(sec43::TestRom(12, 14));
  ret.push_back(sec51::TestRom(123));
  ret.push_back(hyrax::a2::TestRom(10));
#else
  // ret.push_back(groth09::sec52::TestRom(17, 20));
  // ret.push_back(groth09::sec53::TestRom(128, 1024*1024));
  // ret.push_back(groth09::sec43::TestRom(256, 1024 * 1024));
  ret.push_back(groth09::sec43::TestRom(120, 1024 * 128));
  // ret.push_back(groth09::sec51::TestRom(123));
  // ret.push_back(hyrax::a2::TestRom(10));
#endif
  std::cout << "summary:\n";
  for (auto i : ret) {
    std::cout << (i ? "success" : "failed") << "\n";
  }
}
}  // namespace groth09