#pragma once

#include "log/tick.h"

extern bool DEBUG_CHECK;

#define DCHECK(X, S)                                                 \
  if (DEBUG_CHECK && !(X)) {                                         \
    std::cout << __FN__ << ":" << __LINE__ << " oops " << S << "\n"; \
    throw std::runtime_error("oops");                                \
  }

#define CHECK(X, S)                                                  \
  if (!(X)) {                                                        \
    std::cout << __FN__ << ":" << __LINE__ << " oops " << S << "\n"; \
    throw std::runtime_error("oops");                                \
  }
