#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./mimc5.h"

namespace vrs {
class Mimc5Gadget : public libsnark::gadget<Fr> {
 public:
  Mimc5Gadget(libsnark::protoboard<Fr>& pb,
              const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix), constants_(Mimc5Const()) {
    plain_.allocate(pb, "plain");
    key_.allocate(pb, "key");
    rounds_x2_.allocate(pb, constants_.size(),
                        FMT(annotation_prefix, " rounds_x2"));
    rounds_x4_.allocate(pb, constants_.size(),
                        FMT(annotation_prefix, " rounds_x4"));
    rounds_x5_.allocate(pb, constants_.size(),
                        FMT(annotation_prefix, " rounds_x5"));

    generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> plain() { return plain_; }

  libsnark::pb_variable<Fr> key() { return key_; }

  libsnark::pb_variable<Fr> result() {
    return rounds_x5_[constants_.size() - 1];
  }

  void Assign(Fr const& plain, Fr const& key) {
    this->pb.val(plain_) = plain;
    this->pb.val(key_) = key;
    generate_r1cs_witness();
  }

 private:
  void generate_r1cs_constraints() {
    auto data = plain_;
    for (size_t i = 0; i < constants_.size(); ++i) {
      auto x1 = data + key_ + constants_[i];
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(x1, x1, rounds_x2_[i]), "x1*x1");
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(rounds_x2_[i], rounds_x2_[i],
                                        rounds_x4_[i]),
          "x2*x2");

      if (i < constants_.size() - 1) {
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(rounds_x4_[i], x1, rounds_x5_[i]),
            "x4*x1");
        data = rounds_x5_[i];
      } else {
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(rounds_x4_[i], x1,
                                          rounds_x5_[i] - key_),
            "x5");
      }
    }
  }

  void generate_r1cs_witness() {
    auto const& plain = this->pb.val(plain_);
    auto const& key = this->pb.val(key_);
    auto data = plain;
    for (size_t i = 0; i < constants_.size(); ++i) {
      auto x1 = data + key + constants_[i];
      auto x2 = x1 * x1;
      auto x4 = x2 * x2;
      auto x5 = x4 * x1;
      this->pb.val(rounds_x2_[i]) = x2;
      this->pb.val(rounds_x4_[i]) = x4;
      if (i < constants_.size() - 1) {
        this->pb.val(rounds_x5_[i]) = x5;
        // std::cout << x5 << "\n";
        data = x5;
      } else {
        this->pb.val(rounds_x5_[i]) = x5 + key;
      }
    }
    assert(Mimc5Enc(plain, key) == pb.val(result()));
  }

 private:
  std::vector<Fr> const& constants_;
  libsnark::pb_variable<Fr> plain_;
  libsnark::pb_variable<Fr> key_;
  libsnark::pb_variable_array<Fr> rounds_x2_;
  libsnark::pb_variable_array<Fr> rounds_x4_;
  libsnark::pb_variable_array<Fr> rounds_x5_;
};

struct Mimc5Scheme {
  Mimc5Scheme() : gadget(pb, "Mimc5Gadget") {    
    pb.set_input_sizes(kPrimaryInputSize);
    constraint_system = pb.get_constraint_system();
  }
  Fr Generate(Fr const& plain, Fr const& key) { return Mimc5Enc(plain, key); }

  int64_t num_variables() const { return (int64_t)pb.num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb.num_constraints(); }

  static const int64_t kPrimaryInputSize = 1;
#ifdef _DEBUG
  static const int64_t kMaxUnitPerZkp = 32;
#else
  static const int64_t kMaxUnitPerZkp = 1024 * 32;
#endif
  static_assert(kMaxUnitPerZkp <= (int64_t)PcBase::kGSize/3,
                "kMaxUnitPerZkp too large");
  
  static std::string const& type() {
    static const std::string a = "mimc5";
    return a;
  };

  libsnark::protoboard<Fr> pb;
  Mimc5Gadget gadget;
  libsnark::r1cs_constraint_system<Fr> constraint_system;
};
}  // namespace vrs
