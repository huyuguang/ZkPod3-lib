#pragma once

#include "./sec51b.h"
#include "./sec51c.h"
#include "./sec43b.h"
#include "hyrax/a2.h"
#include "hyrax/a3.h"

namespace groth09 {

struct OrdinaryPolicy {
  using HyraxA = hyrax::A2;
  using Sec51 = Sec51b;
  using Sec53 = Sec53b<Sec51>;
  using Sec43 = Sec43b<Sec53, HyraxA>;
};

struct SuccinctPolicy {
  using HyraxA = hyrax::A3;
  using Sec51 = Sec51c;
  using Sec53 = Sec53b<Sec51>;
  using Sec43 = Sec43b<Sec53, HyraxA>;
};
}