#pragma once

#include "./details.h"
#include "circuit/mimc5_gadget.h"
#include "circuit/sha256c_gadget.h"

namespace clink {

struct VrsSha256cScheme {
  VrsSha256cScheme() : gadget(pb, "Sha256cGadget") {
    pb.set_input_sizes(kPrimaryInputSize);
    constraint_system = pb.get_constraint_system();
  }
  static Fr Generate(Fr const& plain, Fr const& key) {
    return circuit::Sha256Enc(plain, key);
  }

  int64_t num_variables() const { return (int64_t)pb.num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb.num_constraints(); }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 8;
#else
  static const int64_t kMaxUnitPerZkp = 512;
#endif
  static_assert(kMaxUnitPerZkp <= (int64_t)PcBase::kGSize / 3,
                "kMaxUnitPerZkp too large");

  static std::string const& type() {
    static const std::string a = "sha256c";
    return a;
  };

  libsnark::protoboard<Fr> pb;
  circuit::Sha256cGadget gadget;
  libsnark::r1cs_constraint_system<Fr> constraint_system;
};

struct VrsMimc5Scheme {
  VrsMimc5Scheme() : gadget(pb, "Mimc5Gadget") {
    pb.set_input_sizes(kPrimaryInputSize);
    constraint_system = pb.get_constraint_system();
  }
  static Fr Generate(Fr const& plain, Fr const& key) {
    return circuit::Mimc5Enc(plain, key);
  }

  int64_t num_variables() const { return (int64_t)pb.num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb.num_constraints(); }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 32;
#else
  static const int64_t kMaxUnitPerZkp = 1024 * 32;
#endif
  static_assert(kMaxUnitPerZkp <= (int64_t)PcBase::kGSize / 3,
                "kMaxUnitPerZkp too large");

  static std::string const& type() {
    static const std::string a = "mimc5";
    return a;
  };

  libsnark::protoboard<Fr> pb;
  circuit::Mimc5Gadget gadget;
  libsnark::r1cs_constraint_system<Fr> constraint_system;
};
}  // namespace clink