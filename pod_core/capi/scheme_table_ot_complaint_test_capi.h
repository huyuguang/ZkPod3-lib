#pragma once

#include <string>

#include "misc/misc.h"

namespace scheme::table::ot_complaint::capi {
bool Test(std::string const& publish_path, std::string const& output_path,
          std::vector<Range> const& demands, std::vector<Range> const& phantoms,
          bool test_evil);
}  // namespace scheme::table::ot_complaint::capi