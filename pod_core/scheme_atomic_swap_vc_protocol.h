#pragma once

#include <string>
#include <vector>
#include "ecc.h"
#include "vrs/vrs.h"

namespace scheme::atomic_swap_vc {
struct Request {
  h256_t bob_nonce;
  std::vector<Range> demands;
};

struct Response {
  h256_t alice_nonce;
  std::vector<G1> k;   // n+1
  std::vector<Fr> m;   // n*s
  std::vector<Fr> vw;  // s
  h256_t vrs_plain_seed;
  Fr vw_com_r;
  std::vector<vrs::Proof> vrs_proofs;  
};

inline bool operator==(Response const& left, Response const& right) {
  return left.alice_nonce == right.alice_nonce && left.k == right.k &&
         left.m == right.m && left.vw == right.vw &&
         left.vw_com_r == right.vw_com_r && left.vrs_proofs == right.vrs_proofs;
}

inline bool operator!=(Response const& left,
                       Response const& right) {
  return !(left == right);
}

struct Receipt {
  G1 h;
  G1 g;
  G1 seed0_com;
};

inline bool operator==(Receipt const& left, Receipt const& right) {
  return left.h == right.h && left.g == right.g &&
         left.seed0_com == right.seed0_com;
}

inline bool operator!=(Receipt const& left, Receipt const& right) {
  return !(left == right); 
}

struct Secret {
  Fr seed0;
  Fr seed0_r;
};

inline bool operator==(Secret const& left, Secret const& right) {
  return left.seed0 == right.seed0 && left.seed0_r == right.seed0_r;
}

inline bool operator!=(Secret const& left, Secret const& right) {
  return !(left == right); 
}

}  // namespace scheme::atomic_swap_vc
