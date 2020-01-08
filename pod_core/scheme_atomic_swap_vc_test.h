#pragma once

#include <string>

#include "misc/misc.h"

namespace scheme::plain::atomic_swap_vc {
bool Test(std::string const& publish_path, std::string const& output_path,
          std::vector<Range> const& demands, bool test_evil);
}  // namespace scheme::plain::atomic_swap_vc

namespace scheme::table::atomic_swap_vc {
bool Test(std::string const& publish_path, std::string const& output_path,
          std::vector<Range> const& demands, bool test_evil);
}  // namespace scheme::table::atomic_swap_vc
