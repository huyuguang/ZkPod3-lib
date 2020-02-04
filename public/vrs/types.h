#pragma once

#include <functional>
#include <memory>

#include "ecc/ecc.h"
#include "groth09/groth09.h"
#include "hyrax/a2.h"
#include "misc/misc.h"

namespace vrs {
struct Cache {
  std::string type;
  int64_t max_unit_per_zkp;
  int64_t count = 0;
  h256_t seed;  // plain seed
  Fr key;
  Fr key_com_r;
  std::vector<std::vector<G1>> var_coms;
  std::vector<std::vector<Fr>> var_coms_r;
};
typedef std::unique_ptr<Cache> CacheUPtr;

inline bool operator==(Cache const& left, Cache const& right) {
  return left.type == right.type &&
         left.max_unit_per_zkp == right.max_unit_per_zkp &&
         left.count == right.count && left.seed == right.seed &&
         left.key == right.key && left.key_com_r == right.key_com_r &&
         left.var_coms == right.var_coms && left.var_coms_r == right.var_coms_r;
}

inline bool operator!=(Cache const& left, Cache const& right) {
  return !(left == right);
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Cache const& t) {
  ar& YAS_OBJECT_NVP("cache", ("type", t.type), ("max", t.max_unit_per_zkp),
                     ("c", t.count), ("s", t.seed), ("k", t.key),
                     ("kcr", t.key_com_r), ("vc", t.var_coms),
                     ("vcr", t.var_coms_r));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Cache& t) {
  ar& YAS_OBJECT_NVP("cache", ("type", t.type), ("max", t.max_unit_per_zkp),
                     ("c", t.count), ("s", t.seed), ("k", t.key),
                     ("kcr", t.key_com_r), ("vc", t.var_coms),
                     ("vcr", t.var_coms_r));
}

struct PublicInput {  // change to function
  PublicInput(int64_t n, std::function<Fr(int64_t i)> f)
      : count(n), get_p(std::move(f)) {}
  int64_t const count;
  std::function<Fr(int64_t i)> get_p;
};

struct SecretInput {
  SecretInput() {}
  SecretInput(Fr const& key, Fr const& key_com_r, Fr const& vw_com_r)
      : key(key), key_com_r(key_com_r), vw_com_r(vw_com_r) {}
  Fr key;
  Fr key_com_r;
  Fr vw_com_r;
};

struct CheckInput {
  CheckInput(std::vector<Fr> const& v, Fr const& vw) : v(v), vw(vw) {}
  std::vector<Fr> const& v;
  Fr const& vw;
};

template <typename Policy>
struct Proof {
  std::vector<G1> var_coms;
  G1 com_vw;
  typename Policy::Sec43::Proof proof_hp;
  typename Policy::HyraxA::Proof proof_ip;
  bool operator==(Proof const& right) const {
    if (var_coms != right.var_coms) std::cout << __LINE__ << "\n";
    if (com_vw != right.com_vw) std::cout << __LINE__ << "\n";
    if (proof_hp != right.proof_hp) std::cout << __LINE__ << "\n";
    if (proof_ip != right.proof_ip) std::cout << __LINE__ << "\n";
    return var_coms == right.var_coms && com_vw == right.com_vw &&
           proof_hp == right.proof_hp && proof_ip == right.proof_ip;
  }
  bool operator!=(Proof const& right) const { return !(*this == right); }
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof<groth09::OrdinaryPolicy> const& t) {
  ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof<groth09::OrdinaryPolicy>& t) {
  ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Proof<groth09::SuccinctPolicy> const& t) {
  ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Proof<groth09::SuccinctPolicy>& t) {
  ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", t.var_coms), ("com_vw", t.com_vw),
                     ("php", t.proof_hp), ("pip", t.proof_ip));
}

struct ProveOutput {
  G1 g;
  G1 h;
  G1 key_com;
};

inline bool operator==(ProveOutput const& left, ProveOutput const& right) {
  return left.g == right.g && left.h == right.h &&
         left.key_com == right.key_com;
}

inline bool operator!=(ProveOutput const& left, ProveOutput const& right) {
  return !(left == right);
}

struct VerifyOutput {
  G1 g;
  G1 h;
  G1 key_com;
};

inline bool operator==(VerifyOutput const& left, VerifyOutput const& right) {
  return left.g == right.g && left.h == right.h &&
         left.key_com == right.key_com;
}

inline bool operator!=(VerifyOutput const& left, VerifyOutput const& right) {
  return !(left == right);
}

}  // namespace vrs