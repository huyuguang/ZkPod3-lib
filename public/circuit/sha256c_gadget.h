#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./sha256c.h"

namespace circuit {
class Sha256cGadget : public libsnark::gadget<Fr> {
 public:
  Sha256cGadget(libsnark::protoboard<Fr>& pb,
                const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    plain_.allocate(pb, "plain");
    key_.allocate(pb, "key");

    dual_plain_.reset(new libsnark::dual_variable_gadget<Fr>(
        this->pb, plain_, 256, "dual_plain"));

    dual_key_.reset(new libsnark::dual_variable_gadget<Fr>(this->pb, key_, 256,
                                                           "dual_key"));

    dual_output_.reset(
        new libsnark::dual_variable_gadget<Fr>(this->pb, 256, "dual_output"));

    sha_left_.reset(new libsnark::digest_variable<Fr>(
        this->pb, 256, dual_plain_->bits, 0, "sha_left"));

    sha_right_.reset(new libsnark::digest_variable<Fr>(
        this->pb, 256, dual_key_->bits, 0, "sha_right"));

    sha_output_.reset(new libsnark::digest_variable<Fr>(
        this->pb, 256, dual_output_->bits, 0, "sha_output"));

    sha_.reset(new libsnark::sha256_two_to_one_hash_gadget<Fr>(
        this->pb, *sha_left_, *sha_right_, *sha_output_, "sha"));

    output_.allocate(pb, "output");  // must in the last
    generate_r1cs_constraints();
  }

  libsnark::pb_variable<Fr> plain() { return plain_; }

  libsnark::pb_variable<Fr> key() { return key_; }

  libsnark::pb_variable<Fr> result() { return output_; }

  void Assign(Fr const& plain, Fr const& key) {
    this->pb.val(plain_) = plain;
    this->pb.val(key_) = key;
    generate_r1cs_witness();
    //assert(this->pb.is_satisfied());
  }

 private:
  void generate_r1cs_constraints() {
    // auto data = plain_;
    dual_plain_->generate_r1cs_constraints(true);
    dual_key_->generate_r1cs_constraints(true);
    dual_output_->generate_r1cs_constraints(true);

    sha_->generate_r1cs_constraints();

    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(output_, 1, dual_output_->packed),
        "output==dual_output.packed");
  }

  void generate_r1cs_witness() {
    dual_plain_->generate_r1cs_witness_from_packed();
    dual_key_->generate_r1cs_witness_from_packed();

    sha_->generate_r1cs_witness();

    dual_output_->generate_r1cs_witness_from_bits();

    this->pb.val(output_) = this->pb.val(dual_output_->packed);
  }

 private:
  libsnark::pb_variable<Fr> plain_;
  libsnark::pb_variable<Fr> key_;
  libsnark::pb_variable<Fr> output_;
  std::unique_ptr<libsnark::dual_variable_gadget<Fr>> dual_plain_;
  std::unique_ptr<libsnark::dual_variable_gadget<Fr>> dual_key_;
  std::unique_ptr<libsnark::dual_variable_gadget<Fr>> dual_output_;

  std::unique_ptr<libsnark::digest_variable<Fr>> sha_left_;
  std::unique_ptr<libsnark::digest_variable<Fr>> sha_right_;
  std::unique_ptr<libsnark::digest_variable<Fr>> sha_output_;
  std::unique_ptr<libsnark::sha256_two_to_one_hash_gadget<Fr>> sha_;
};

inline bool TestSha256cGadget() {
  libsnark::protoboard<Fr> pb;
  Sha256cGadget gadget(pb, "Sha256cGadget");

  Fr plain = FrRand();
  Fr key = FrRand();

  Fr ret = Sha256Enc(plain, key);
  std::cout << ret << "\n";

  gadget.Assign(plain, key);

  auto output = pb.val(gadget.result());
  std::cout << "output: " << output << "\n";

  assert(ret == output);
  return ret == output;
}

}  // namespace circuit