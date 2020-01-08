#pragma once

#include <stdint.h>

#include <memory>
#include <string>

#include "misc/misc.h"
#include "bp/bp.h"
#include "scheme/bulletin_table.h"
#include "utils/mkl_tree.h"
#include "scheme/public_misc.h"
#include "vrf/vrf.h"
#include "scheme/vrf_meta.h"

namespace scheme::table {

class AliceData {
 public:
  AliceData(std::string const& publish_path);
  Bulletin const& bulletin() const { return bulletin_; }
  VrfMeta const& vrf_meta() const { return vrf_meta_; }
  vrf::Pk<> const& vrf_pk() const { return vrf_pk_; }
  vrf::Sk<> const& vrf_sk() const { return vrf_sk_; }
  std::vector<G1> const& sigmas() const { return sigmas_; }
  std::vector<Fr> const& m() const { return m_; }

 public:
  VrfKeyMeta const* GetKeyMetaByName(std::string const& name);

 private:
  std::string const publish_path_;
  scheme::table::Bulletin bulletin_;
  vrf::Pk<> vrf_pk_;
  vrf::Sk<> vrf_sk_;  // secret
  VrfMeta vrf_meta_;
  std::vector<bp::P1Proof> vrf_key_bp_proofs_;
  std::vector<G1> sigmas_;
  mkl::Tree sigma_mkl_tree_;
  std::vector<Fr> m_;  // secret
  std::vector<std::vector<Fr>> key_m_;
};

typedef std::shared_ptr<AliceData> AliceDataPtr;
}  // namespace scheme::table