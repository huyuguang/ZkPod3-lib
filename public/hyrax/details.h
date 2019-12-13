#pragma once

#include <stdint.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "../ecc.h"
#include "../ecc_pub.h"
#include "../fst.h"
#include "../misc.h"
#include "../multiexp.h"
#include "../parallel.h"
#include "../pds_pub.h"
#include "../tick.h"
#include "../vectorop.h"

namespace hyrax::details {

inline G1 const& GetPdsG() {
  auto const& pds_pub = GetPdsPub();
  return pds_pub.g()[0];
}

inline G1 const& GetPdsH() {
  auto const& pds_pub = GetPdsPub();
  return pds_pub.h();
}

inline G1 ComputeCommitment(std::vector<Fr> const& x, Fr const& r) {
  auto const& pds_pub = GetPdsPub();
  // Tick tick(__FUNCTION__, std::to_string(x.size()));
  assert(PdsPub::kGSize >= x.size());
  auto get_g = [&pds_pub](int64_t i) -> G1 const& {
    return i ? pds_pub.g()[i - 1] : pds_pub.h();
  };
  auto get_f = [&x, &r](int64_t i) -> Fr const& { return i ? x[i - 1] : r; };
  return MultiExpBdlo12Inner<G1>(get_g, get_f, x.size() + 1);
}

inline G1 ComputeCommitment(Fr const& x, Fr const& r) {
  return ComputeCommitment(std::vector<Fr>{x}, r);
}

}  // namespace hyrax::details