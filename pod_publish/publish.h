#pragma once

#include <stdint.h>

#include <string>
#include <vector>

#include "scheme/public_misc.h"
#include "scheme/plain_misc.h"
#include "scheme/table_misc.h"

bool PublishTable(std::string publish_file, std::string output_path,
                  scheme::table::Type table_type,
                  std::vector<uint64_t> vrf_colnums_index,
                  std::vector<bool> unique_key);

bool PublishPlain(std::string publish_file, std::string output_path,
                  uint64_t column_num);