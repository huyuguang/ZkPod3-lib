#pragma once

#include <memory>
#include <vector>

#include "clink/parallel_r1cs.h"

namespace clink::vgg16 {
struct BaseR1csSec {
  std::unique_ptr<R1csInfo> r1cs_info;
  std::vector<Fr> com_w_r;
  std::vector<std::vector<Fr>> mutable w;
};

}  // namespace clink::vgg16