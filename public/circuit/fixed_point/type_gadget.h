#pragma once

#include "./sign_gadget.h"

namespace circuit::fixed_point {

template <size_t W> 
using TypeGadget = SignGadget<W>;
}  // namespace circuit::fixed_point