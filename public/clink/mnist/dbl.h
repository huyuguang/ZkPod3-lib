#pragma once

#include "../details.h"

namespace clink {
namespace mnist {
namespace dbl {
struct ConvPara {
  std::array<double, 9> coef;
  double bias;
};

struct Dense1Para {
  std::array<double, 13 * 13 * 5> coef;
  double bias;
};

struct Dense2Para {
  std::array<double, 10> coef;
  double bias;
};

struct Para {
  std::array<ConvPara, 5> conv;
  std::array<Dense1Para, 10> dense1;
  std::array<Dense2Para, 10> dense2;
};

typedef std::array<uint8_t, 28 * 28> RawData;
typedef std::array<double, 28 * 28> UniData;

typedef std::array<double, 26 * 26> ConvData;
typedef std::array<ConvData, 5> ConvsData;

typedef std::array<double, 13 * 13> PoolingData;
typedef std::array<PoolingData, 5> PoolingsData;

typedef std::array<double, 13 * 13 * 5> FlatData;

typedef std::array<double, 10> Dense1Data;
typedef std::array<double, 10> Dense2Data;

}  // namespace dbl
}  // namespace mnist
}  // namespace clink