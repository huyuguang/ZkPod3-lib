#pragma once

#include "../details.h"
#include "circuit/mimc5_gadget.h"
#include "circuit/poseidon_gadget.h"
#include "circuit/sha256c_gadget.h"

namespace clink {

struct VrsSha256cScheme {
  VrsSha256cScheme() {
    libsnark::protoboard<Fr> pb;
    auto gadget = CreateGadget(pb);
    r1cs_info.reset(new R1csInfo(pb));
  }

  std::unique_ptr<circuit::Sha256cGadget> CreateGadget(
      libsnark::protoboard<Fr>& pb) {
    std::unique_ptr<circuit::Sha256cGadget> ret;
    ret.reset(new circuit::Sha256cGadget(pb, "Sha256cGadget"));
    pb.set_input_sizes(kPrimaryInputSize);
    return ret;
  }

  static Fr Generate(Fr const& plain, Fr const& key) {
    return circuit::Sha256Enc(plain, key);
  }

  int64_t num_variables() const { return r1cs_info->num_variables; }

  int64_t num_constraints() const { return r1cs_info->num_constraints; }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 8;
#else
  static const int64_t kMaxUnitPerZkp = 128;
#endif

  static std::string const& type() {
    static const std::string a = "sha256c";
    return a;
  };

  std::unique_ptr<R1csInfo> r1cs_info;
};

struct VrsMimc5Scheme {
  VrsMimc5Scheme() {
    libsnark::protoboard<Fr> pb;
    auto gadget = CreateGadget(pb);
    r1cs_info.reset(new R1csInfo(pb));
  }

  std::unique_ptr<circuit::Mimc5Gadget> CreateGadget(
      libsnark::protoboard<Fr>& pb) {
    std::unique_ptr<circuit::Mimc5Gadget> ret;
    ret.reset(new circuit::Mimc5Gadget(pb, "Mimc5Gadget"));
    pb.set_input_sizes(kPrimaryInputSize);
    return ret;
  }

  static Fr Generate(Fr const& plain, Fr const& key) {
    return circuit::Mimc5Enc(plain, key);
  }

  int64_t num_variables() const { return r1cs_info->num_variables; }

  int64_t num_constraints() const { return r1cs_info->num_constraints; }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 32;
#else
  static const int64_t kMaxUnitPerZkp = 1024 * 16;
#endif

  static std::string const& type() {
    static const std::string a = "mimc5";
    return a;
  };

  std::unique_ptr<R1csInfo> r1cs_info;
};

struct VrsPoseidonScheme {
  VrsPoseidonScheme() {
    libsnark::protoboard<Fr> pb;
    auto gadget = CreateGadget(pb);
    r1cs_info.reset(new R1csInfo(pb));
  }

  std::unique_ptr<circuit::VrsPoseidon> CreateGadget(
      libsnark::protoboard<Fr>& pb) {
    std::unique_ptr<circuit::VrsPoseidon> ret;
    ret.reset(new circuit::VrsPoseidon(pb, "PoseidonGadget"));
    pb.set_input_sizes(kPrimaryInputSize);
    return ret;
  }

  static Fr Generate(Fr const& plain, Fr const& key) {
    std::array<Fr, 2> inputs{{plain, key}};
    auto ret = circuit::Poseidon<5, 1, 6, 52, 2, 1>(inputs);
    assert(ret.size() == 1);
    assert(ret[0] == circuit::VrsPoseidon::permute(plain, key));
    return ret[0];
  }

  int64_t num_variables() const { return r1cs_info->num_variables; }

  int64_t num_constraints() const { return r1cs_info->num_constraints; }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 32;
#else
  static const int64_t kMaxUnitPerZkp = 1024 * 16;
#endif
  // static_assert(kMaxUnitPerZkp <= pc::Base::GSize() / 3,
  //              "kMaxUnitPerZkp too large");

  static std::string const& type() {
    static const std::string a = "poseidon";
    return a;
  };

  std::unique_ptr<R1csInfo> r1cs_info;
};
}  // namespace clink