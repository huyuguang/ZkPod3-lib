#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./misc.h"

namespace vrs {
class BaseOnewayGadget : public libsnark::gadget<Fr> {
 public:
  BaseOnewayGadget(libsnark::protoboard<Fr>& pb,
             const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix) {}

  virtual libsnark::pb_variable<Fr> plain() = 0;

  virtual libsnark::pb_variable<Fr> key() = 0;

  virtual libsnark::pb_variable<Fr> result() = 0;

  virtual void Assign(Fr const& plain, Fr const& key) = 0;
};

class BaseOnewayScheme {
 public:
  virtual int64_t MaxUnitPerZkp() = 0;
  virtual Fr Generate(Fr const& plain, Fr const& key) = 0;

 public:
  libsnark::protoboard<Fr>& pb() { return pb_; }
  libsnark::pb_variable<Fr> plain() { return gadget_->plain(); }
  libsnark::pb_variable<Fr> key() { return gadget_->key(); }
  libsnark::pb_variable<Fr> result() { return gadget_->result(); }
  void Assign(Fr const& plain, Fr const& key) {
    return gadget_->Assign(plain, key);
  }

 protected:
  libsnark::protoboard<Fr> pb_;
  std::unique_ptr<BaseOnewayGadget> gadget_;
};
}  // namespace vrs