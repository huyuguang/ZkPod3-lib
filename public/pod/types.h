#pragma once

#include <functional>

#include "ecc/ecc.h"
#include "vrs/vrs.h"

namespace pod {

typedef std::function<Fr const&(int64_t i, int64_t j)> GetMatrix;
typedef std::function<Fr const&(int64_t i)> GetR;
typedef std::function<G1 const&(int64_t i)> GetCom;

// com is base on g_offset = 0
struct CommitedData {
  int64_t n = 0;
  int64_t s = 0;
  GetMatrix get_m;
  GetR get_r;
  GetCom get_com;
};

template <typename Policy>
struct ProvedData {
  std::vector<G1> k;   // size = n+1, k=com(v)
  std::vector<Fr> em;  // size = n*(s+1), encrypted r&m
  std::vector<Fr> vw;  // s + 1
  h256_t vrs_plain_seed;
  Fr vw_com_r;
  std::vector<vrs::Proof<Policy>> vrs_proofs;
  bool operator==(ProvedData const& b) const {
    return k == b.k && em == b.em && vw == b.vw &&
           vrs_plain_seed == b.vrs_plain_seed && vw_com_r == b.vw_com_r &&
           vrs_proofs == b.vrs_proofs;
  }

  bool operator!=(ProvedData const& b) const { return !(*this == b); }
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, typename ProvedData<groth09::OrdinaryPolicy> const& t) {
  ar& YAS_OBJECT_NVP("pod.pd", ("k", t.k), ("em", t.em), ("vw", t.vw),
                     ("seed", t.vrs_plain_seed), ("r", t.vw_com_r),
                     ("p", t.vrs_proofs));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, typename ProvedData<groth09::OrdinaryPolicy>& t) {
  ar& YAS_OBJECT_NVP("pod.pd", ("k", t.k), ("em", t.em), ("vw", t.vw),
                     ("seed", t.vrs_plain_seed), ("r", t.vw_com_r),
                     ("p", t.vrs_proofs));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, typename ProvedData<groth09::SuccinctPolicy> const& t) {
  ar& YAS_OBJECT_NVP("pod.pd", ("k", t.k), ("em", t.em), ("vw", t.vw),
                     ("seed", t.vrs_plain_seed), ("r", t.vw_com_r),
                     ("p", t.vrs_proofs));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, typename ProvedData<groth09::SuccinctPolicy>& t) {
  ar& YAS_OBJECT_NVP("pod.pd", ("k", t.k), ("em", t.em), ("vw", t.vw),
                     ("seed", t.vrs_plain_seed), ("r", t.vw_com_r),
                     ("p", t.vrs_proofs));
}

struct Receipt {
  G1 h;
  G1 g;
  G1 seed0_com;
};

inline bool operator==(Receipt const& a, Receipt const& b) {
  return a.h == b.h && a.g == b.g && a.seed0_com == b.seed0_com;
}
inline bool operator!=(Receipt const& a, Receipt const& b) { return !(a == b); }

struct Secret {
  Fr seed0;
  Fr seed0_com_r;
};

template <typename VrsScheme, typename Policy>
struct ProveOutput {
  ProvedData<Policy> proved_data;
  Receipt receipt;
  Secret secret;
  vrs::AutoCacheFileUPtr<VrsScheme> auto_cache;
};

struct VerifyOutput {
  Receipt receipt;
  std::vector<Fr> plain;
  std::vector<Fr> w;
  Fr sigma_vw;
};
}  // namespace pod