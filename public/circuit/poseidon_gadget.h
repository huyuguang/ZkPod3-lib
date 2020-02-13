#pragma once

#include <stdlib.h>

#include <iostream>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./poseidon.h"

// NOTE: Fr must initialize as 0!!!
namespace circuit {

namespace utils {
inline std::vector<libsnark::linear_combination<Fr>> VariableArrayT_to_lc(
    const libsnark::pb_variable_array<Fr>& in_vars) {
  std::vector<libsnark::linear_combination<Fr>> ret;
  ret.reserve(in_vars.size());
  for (const auto& var : in_vars) {
    ret.emplace_back(var);
  }
  return ret;
}

inline Fr lc_val(const libsnark::protoboard<Fr>& pb,
                 const libsnark::linear_combination<Fr>& in_lc) {
  Fr sum = 0;
  for (const auto& term : in_lc.terms) {
    sum += term.coeff * pb.val(libsnark::pb_variable<Fr>(term.index));
  }
  return sum;
}
inline std::vector<Fr> vals(
    const libsnark::protoboard<Fr>& pb,
    const std::vector<libsnark::linear_combination<Fr>>& in_lcs) {
  std::vector<Fr> ret;
  ret.reserve(in_lcs.size());
  for (const auto& lc : in_lcs) {
    ret.emplace_back(lc_val(pb, lc));
  }
  return ret;
}

inline std::vector<Fr> vals(const libsnark::protoboard<Fr>& pb,
                            const libsnark::pb_variable_array<Fr>& in_vars) {
  return in_vars.get_vals(pb);
}

inline const libsnark::pb_variable_array<Fr> make_var_array(
    libsnark::protoboard<Fr>& in_pb, size_t n, const std::string& annotation) {
  libsnark::pb_variable_array<Fr> x;
  x.allocate(in_pb, n, annotation);
  return x;
}

inline const libsnark::pb_variable_array<Fr> make_var_array(
    libsnark::protoboard<Fr>& in_pb, const std::string& annotation,
    std::vector<Fr> values) {
  auto vars = make_var_array(in_pb, values.size(), annotation);
  for (unsigned i = 0; i < values.size(); i++) {
    in_pb.val(vars[i]) = values[i];
  }
  return vars;
}
}  // namespace utils

class FifthPower_gadget : public libsnark::gadget<Fr> {
 public:
  libsnark::pb_variable<Fr> x2;
  libsnark::pb_variable<Fr> x4;
  libsnark::pb_variable<Fr> x5;

  FifthPower_gadget(libsnark::protoboard<Fr>& pb,
                    const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    x2.allocate(pb, FMT(annotation_prefix, ".x2"));
    x4.allocate(pb, FMT(annotation_prefix, ".x4"));
    x5.allocate(pb, FMT(annotation_prefix, ".x5"));
  }

  void generate_r1cs_constraints(
      const libsnark::linear_combination<Fr>& x) const {
    pb.add_r1cs_constraint(libsnark::r1cs_constraint<Fr>(x, x, x2),
                           ".x^2 = x * x");
    pb.add_r1cs_constraint(libsnark::r1cs_constraint<Fr>(x2, x2, x4),
                           ".x^4 = x2 * x2");
    pb.add_r1cs_constraint(libsnark::r1cs_constraint<Fr>(x, x4, x5),
                           ".x^5 = x * x4");
  }

  void generate_r1cs_witness(const Fr& val_x) const {
    const auto val_x2 = val_x * val_x;
    const auto val_x4 = val_x2 * val_x2;
    const auto val_x5 = val_x4 * val_x;
    this->pb.val(x2) = val_x2;
    this->pb.val(x4) = val_x4;
    this->pb.val(x5) = val_x5;
  }

  const libsnark::pb_variable<Fr>& result() const { return x5; }
};

/**
 * One round of the Base permutation:
 *
 *    - takes a state of `t` elements
 *    - adds the round constant to each element in the state
 *    - performs exponentiation on the first `n` elements of the state
 *    - creates `o` outputs, mixed using a matrix vector transform
 *
 * This generic version can be used as either a 'full', 'partial' or 'last'
 * round. It avoids computing as many constraints as is possible, given all the
 * information.
 */
template <unsigned param_t, unsigned nSBox, unsigned nInputs, unsigned nOutputs>
class Poseidon_Round : public libsnark::gadget<Fr> {
 public:
  const Fr& C_i;
  const std::vector<Fr>& M;
  const std::vector<libsnark::linear_combination<Fr>> state;
  const std::vector<FifthPower_gadget> sboxes;
  const std::vector<libsnark::linear_combination<Fr>> outputs;

  static std::vector<FifthPower_gadget> make_sboxes(
      libsnark::protoboard<Fr>& in_pb, const std::string& annotation_prefix) {
    std::vector<FifthPower_gadget> ret;

    ret.reserve(nSBox);
    for (unsigned h = 0; h < nSBox; h++) {
      ret.emplace_back(in_pb, FMT(annotation_prefix, ".sbox[%u]", h));
    }

    return ret;
  }

  static std::vector<libsnark::linear_combination<Fr>> make_outputs(
      libsnark::protoboard<Fr>& in_pb, const Fr& in_C_i,
      const std::vector<Fr>& in_M,
      const std::vector<libsnark::linear_combination<Fr>>& in_state,
      const std::vector<FifthPower_gadget>& in_sboxes) {
    (void)in_pb;
    std::vector<libsnark::linear_combination<Fr>> ret;

    for (unsigned i = 0; i < nOutputs; i++) {
      const unsigned M_offset = i * param_t;

      // Any element which isn't passed through an sbox
      // Can be accumulated separately as part of the constant term
      Fr constant_term(0);
      for (unsigned j = nSBox; j < param_t; j++) {
        constant_term += in_C_i * in_M[M_offset + j];
      }

      libsnark::linear_combination<Fr> lc;
      lc.terms.reserve(param_t);
      if (nSBox < param_t) {
        lc.add_term(libsnark::pb_variable<Fr>(0), constant_term);
      }

      // Add S-Boxes to the output row
      for (unsigned s = 0; s < nSBox; s++) {
        lc.add_term(in_sboxes[s].result(), in_M[M_offset + s]);
      }

      // Then add inputs (from the state) multiplied by the matrix element
      for (unsigned k = nSBox; k < nInputs; k++) {
        lc = lc + (in_state[k] * in_M[M_offset + k]);
      }

      ret.emplace_back(lc);
    }
    return ret;
  }

  Poseidon_Round(libsnark::protoboard<Fr>& in_pb, const Fr& in_C_i,
                 const std::vector<Fr>& in_M,
                 const libsnark::pb_variable_array<Fr>& in_state,
                 const std::string& annotation_prefix)
      : Poseidon_Round(in_pb, in_C_i, in_M,
                       utils::VariableArrayT_to_lc(in_state),
                       annotation_prefix) {}

  Poseidon_Round(libsnark::protoboard<Fr>& in_pb, const Fr& in_C_i,
                 const std::vector<Fr>& in_M,
                 const std::vector<libsnark::linear_combination<Fr>>& in_state,
                 const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(in_pb, annotation_prefix),
        C_i(in_C_i),
        M(in_M),
        state(in_state),
        sboxes(make_sboxes(in_pb, annotation_prefix)),
        outputs(make_outputs(in_pb, in_C_i, in_M, in_state, sboxes)) {
    static_assert(nInputs <= param_t, "nInputs <= param_t");
    static_assert(nOutputs <= param_t, "nOutputs <= param_t");
  }

  void generate_r1cs_witness() const {
    for (unsigned h = 0; h < nSBox; h++) {
      auto value = C_i;
      if (h < nInputs) {
        value += utils::lc_val(this->pb, state[h]);  // this->pb.val(state[h]);
      }
      sboxes[h].generate_r1cs_witness(value);
    }
  }

  void generate_r1cs_constraints() const {
    for (unsigned h = 0; h < nSBox; h++) {
      if (h < nInputs) {
        sboxes[h].generate_r1cs_constraints(state[h] + C_i);
      } else {
        sboxes[h].generate_r1cs_constraints(C_i);
      }
    }
  }
};

template <unsigned param_t, unsigned param_c, unsigned param_F,
          unsigned param_P, unsigned nInputs, unsigned nOutputs,
          bool constrainOutputs = true>
class Poseidon_gadget_T : public libsnark::gadget<Fr> {
 protected:
  typedef Poseidon_Round<param_t, param_t, nInputs, param_t>
      FirstRoundT;  // ingests `nInput` elements, expands to `t` elements using
                    // round constants
  typedef Poseidon_Round<param_t, param_c, param_t, param_t>
      PartialRoundT;  // partial round only runs sbox on `c` elements (capacity)
  typedef Poseidon_Round<param_t, param_t, param_t, param_t>
      FullRoundT;  // full bandwidth
  typedef Poseidon_Round<param_t, param_t, param_t, nOutputs>
      LastRoundT;  // squeezes state into `nOutputs`

  typedef const std::vector<libsnark::linear_combination<Fr>>& lc_outputs_t;
  typedef const libsnark::linear_combination<Fr>& lc_output_t;
  typedef const libsnark::pb_variable<Fr>& var_output_t;
  typedef const libsnark::pb_variable_array<Fr>& var_outputs_t;

  static constexpr unsigned partial_begin = (param_F / 2);
  static constexpr unsigned partial_end = (partial_begin + param_P);
  static constexpr unsigned total_rounds = param_F + param_P;

 public:
  const libsnark::pb_variable_array<Fr>& inputs;
  const PoseidonConstants& constants;

  FirstRoundT first_round;
  std::vector<FullRoundT> prefix_full_rounds;
  std::vector<PartialRoundT> partial_rounds;
  std::vector<FullRoundT> suffix_full_rounds;
  LastRoundT last_round;

  // When `constrainOutputs==true`, need variables to store outputs
  const libsnark::pb_variable_array<Fr> _output_vars;

  template <typename T>
  static const std::vector<T> make_rounds(
      unsigned n_begin, unsigned n_end, libsnark::protoboard<Fr>& pb,
      const std::vector<libsnark::linear_combination<Fr>>& inputs,
      const PoseidonConstants& constants,
      const std::string& annotation_prefix) {
    std::vector<T> result;
    result.reserve(n_end - n_begin);

    for (unsigned i = n_begin; i < n_end; i++) {
      const auto& state = (i == n_begin) ? inputs : result.back().outputs;
      result.emplace_back(pb, constants.C[i], constants.M, state,
                          FMT(annotation_prefix, ".round[%u]", i));
    }

    return result;
  }

  static std::vector<Fr> permute(std::vector<Fr> const& inputs) {
    libsnark::protoboard<Fr> pb;

    assert(inputs.size() == nInputs);
    auto var_inputs = utils::make_var_array(pb, "input", inputs);

    Poseidon_gadget_T<param_t, param_c, param_F, param_P, nInputs, nOutputs>
        gadget(pb, var_inputs, "gadget");
    gadget.generate_r1cs_witness();

    // auto a1 = utils::vals(pb, gadget.first_round.outputs);
    // misc::PrintVector(a1);
    // auto a2 = utils::vals(pb, gadget.prefix_full_rounds.back().outputs);
    // misc::PrintVector(a2);
    // auto a3 = utils::vals(pb, gadget.partial_rounds.back().outputs);
    // misc::PrintVector(a3);
    // auto a4 = utils::vals(pb, gadget.suffix_full_rounds.back().outputs);
    // misc::PrintVector(a4);

    /*
    // Debugging statements
    gadget.generate_r1cs_constraints();

    unsigned i = 0;
    const auto first_outputs = gadget.first_round.outputs;
    for( unsigned j = 0; j < first_outputs.size(); j++ ) {
            std::cout << "o[" << i << "][" << j << "] = ";
            pb.val(first_outputs[j]).print();
    }
    std::cout << std::endl;

    for( const auto prefix_round : gadget.prefix_full_rounds )
    {
            i += 1;
            const auto outputs = prefix_round.outputs;
            for( unsigned j = 0; j < outputs.size(); j++ ) {
                    std::cout << "o[" << i << "][" << j << "] = ";
                    pb.val(outputs[j]).print();
            }
    }
    std::cout << std::endl;

    for( const auto partial_round : gadget.partial_rounds )
    {
            i += 1;
            const auto outputs = partial_round.outputs;
            for( unsigned j = 0; j < outputs.size(); j++ ) {
                    std::cout << "o[" << i << "][" << j << "] = ";
                    pb.val(outputs[j]).print();
            }
    }
    std::cout << std::endl;

    for( const auto suffix_round : gadget.suffix_full_rounds )
    {
            i += 1;
            const auto outputs = suffix_round.outputs;
            for( unsigned j = 0; j < outputs.size(); j++ ) {
                    std::cout << "o[" << i << "][" << j << "] = ";
                    pb.val(outputs[j]).print();
            }
    }
    std::cout << std::endl;

    const auto last_outputs = gadget.last_round.outputs;
    for( unsigned j = 0; j < last_outputs.size(); j++ ) {
            std::cout << "o[" << i << "][" << j << "] = ";
            pb.val(last_outputs[j]).print();
    }
    std::cout << std::endl;

    if( ! pb.is_satisfied() ) {
            std::cerr << "Not satisfied\n";
    }

    std::cout << pb.num_constraints() << " constraints" << std::endl;
    */

    return utils::vals(pb, gadget.results());
  }

  Poseidon_gadget_T(libsnark::protoboard<Fr>& pb,
                    const libsnark::pb_variable_array<Fr>& in_inputs,
                    const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        inputs(in_inputs),
        constants(poseidon_params<param_t, param_F, param_P>()),
        first_round(pb, constants.C[0], constants.M, in_inputs,
                    FMT(annotation_prefix, ".round[0]")),
        prefix_full_rounds(
            make_rounds<FullRoundT>(1, partial_begin, pb, first_round.outputs,
                                    constants, annotation_prefix)),
        partial_rounds(make_rounds<PartialRoundT>(
            partial_begin, partial_end, pb, prefix_full_rounds.back().outputs,
            constants, annotation_prefix)),
        suffix_full_rounds(make_rounds<FullRoundT>(
            partial_end, total_rounds - 1, pb, partial_rounds.back().outputs,
            constants, annotation_prefix)),
        last_round(pb, constants.C.back(), constants.M,
                   suffix_full_rounds.back().outputs,
                   FMT(annotation_prefix, ".round[%u]", total_rounds - 1)),
        _output_vars(constrainOutputs
                         ? utils::make_var_array(pb, nOutputs, ".output")
                         : libsnark::pb_variable_array<Fr>()) {}

  template <bool x = constrainOutputs>
  typename std::enable_if<!x, lc_outputs_t>::type results() const {
    return last_round.outputs;
  }

  template <bool x = constrainOutputs>
  typename std::enable_if<x, var_outputs_t>::type results() const {
    return _output_vars;
  }

  template <bool x = constrainOutputs, unsigned n = nOutputs>
  typename std::enable_if<!x && n == 1, lc_output_t>::type result() const {
    return last_round.outputs[0];
  }

  template <bool x = constrainOutputs, unsigned n = nOutputs>
  typename std::enable_if<x && n == 1, var_output_t>::type result() const {
    return _output_vars[0];
  }

  void generate_r1cs_constraints() const {
    first_round.generate_r1cs_constraints();

    for (auto& prefix_round : prefix_full_rounds) {
      prefix_round.generate_r1cs_constraints();
    }

    for (auto& partial_round : partial_rounds) {
      partial_round.generate_r1cs_constraints();
    }

    for (auto& suffix_round : suffix_full_rounds) {
      suffix_round.generate_r1cs_constraints();
    }

    last_round.generate_r1cs_constraints();

    if (constrainOutputs) {
      unsigned i = 0;
      for (const auto& lc : last_round.outputs) {
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(lc, libsnark::pb_variable<Fr>(0),
                                          _output_vars[i]),
            FMT(this->annotation_prefix, ".output[%u] = last_round.output[%u]",
                i, i));
        i += 1;
      }
    }
  }

  void generate_r1cs_witness() const {
    first_round.generate_r1cs_witness();

    for (auto& prefix_round : prefix_full_rounds) {
      prefix_round.generate_r1cs_witness();
    }

    for (auto& partial_round : partial_rounds) {
      partial_round.generate_r1cs_witness();
    }

    for (auto& suffix_round : suffix_full_rounds) {
      suffix_round.generate_r1cs_witness();
    }

    last_round.generate_r1cs_witness();

    // When outputs are constrained, fill in the variable
    if (constrainOutputs) {
      unsigned i = 0;
      for (const auto& value : utils::vals(pb, last_round.outputs)) {
        this->pb.val(_output_vars[i++]) = value;
      }
    }
  }
};

template <unsigned nInputs, unsigned nOutputs, bool constrainOutputs = true>
using Poseidon128 =
    Poseidon_gadget_T<6, 1, 8, 57, nInputs, nOutputs, constrainOutputs>;

//
class VrsPoseidon : public libsnark::gadget<Fr> {
 public:
  typedef Poseidon_gadget_T<5, 1, 6, 52, 2, 1, true> Base;
  VrsPoseidon(libsnark::protoboard<Fr>& pb,
              const std::string& annotation_prefix)
      : libsnark::gadget<Fr>(pb, annotation_prefix) {
    inputs.allocate(pb, 2, FMT(annotation_prefix, " inputs"));
    poseidon.reset(new Base(pb, inputs, annotation_prefix));
    poseidon->generate_r1cs_constraints();
  }

  static Fr permute(Fr const& plain, Fr const& key) {
    auto outputs = Base::permute(std::vector<Fr>{{plain, key}});
    assert(outputs.size() == 1);
    return outputs[0];
  }

  libsnark::pb_variable<Fr> plain() { return inputs[0]; }

  libsnark::pb_variable<Fr> key() { return inputs[1]; }

  libsnark::pb_variable<Fr> result() { return poseidon->result(); }

  void Assign(Fr const& plain, Fr const& key) {
    std::vector<Fr> vals{{plain, key}};
    this->inputs.fill_with_field_elements(this->pb, vals);
    poseidon->generate_r1cs_witness();
  }

  libsnark::pb_variable_array<Fr> inputs;
  std::unique_ptr<Base> poseidon;
};

}  // namespace circuit
