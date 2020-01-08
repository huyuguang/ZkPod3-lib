#pragma once

#include <string>

#include "misc/misc.h"

namespace scheme::plain::complaint {
bool Test(std::string const& publish_path, std::string const& output_path,
          std::vector<Range> const& demands, bool test_evil);
}  // namespace scheme::plain::complaint

namespace scheme::table::complaint {
bool Test(std::string const& publish_path, std::string const& output_path,
          std::vector<Range> const& demands, bool test_evil);
}  // namespace scheme::table::complaint
