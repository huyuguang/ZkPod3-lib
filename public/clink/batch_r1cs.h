#pragma once

#include <libsnark/gadgetlib1/protoboard.hpp>

#include "./details.h"
#include "./parallel_r1cs.h"

// r1cs_info: r1cs of one circuit. m=r1cs_info.num_constraints()
// s=r1cs_info.num_variables()
// w: matrix<Fr, s, n>, all var(witness) of the circuits.
// com_w_r: array<Fr, s>, random fr.
// com_w: array<Fr, s>, com_w[i] = com(w[i], com_w_r[i]).
// for open n instance of circuit, open com_w, prove consistency of the com_w

namespace clink {

template <typename Policy>
struct BatchR1cs {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using Sec43 = typename Policy::Sec43;
  using Proof = typename Sec43::Proof;
  using BaseR1cs = ParallelR1cs<Policy>;
  using ProveInput = typename BaseR1cs::ProveInput;
  using VerifyInput = typename BaseR1cs::VerifyInput;

  static void UpdateSeed(h256_t& seed,
                         std::vector<ProveInput*> const& sorted_inputs) {
    // TODO
    (void)seed;
    (void)sorted_inputs;
  }

  static void UpdateSeed(h256_t& seed,
                         std::vector<VerifyInput*> const& sorted_inputs) {
    // TODO
    (void)seed;
    (void)sorted_inputs;
  }

  // w: s*n
  static void Prove(Proof& proof, h256_t seed,
                    std::vector<ProveInput*>&& inputs) {
    Tick tick(__FN__);
    if (inputs.empty()) return;

    auto const& get_g = inputs[0]->get_g;
    for (auto const& i : inputs) {
      if (i->get_g(0) != get_g(0)) throw std::runtime_error("oops");
    }

    std::sort(inputs.begin(), inputs.end(), [](auto const& a, auto const& b) {
      return a->unique_tag < b->unique_tag;
    });

    UpdateSeed(seed, inputs);
    // std::cout << Tick::GetIndentString() << " " << misc::HexToStr(seed) <<
    // "\n";

    std::vector<typename Sec43::CommitmentPub> com_pubs(inputs.size());
    std::vector<typename Sec43::CommitmentSec> com_secs(inputs.size());

    auto pf = [&inputs, &com_pubs, &com_secs](int64_t i) {
      auto const& input = *inputs[i];
      auto& com_pub = com_pubs[i];
      auto& com_sec = com_secs[i];
      BaseR1cs::BuildHpCom(input.m, input.n, input.com_w, input.com_w_r,
                           input.constraints(), input.get_g, com_pub, com_sec);
    };
    parallel::For(inputs.size(), pf);

    size_t combined_m = 0;
    for (auto const& input : inputs) {
      combined_m += input->m;
    }
    std::vector<std::vector<Fr>> combined_x(combined_m);
    std::vector<std::vector<Fr>> combined_y(combined_m);
    std::vector<std::vector<Fr>> combined_z(combined_m);
    std::vector<G1> combined_a(combined_m);
    std::vector<G1> combined_b(combined_m);
    std::vector<G1> combined_c(combined_m);
    std::vector<Fr> combined_r(combined_m);
    std::vector<Fr> combined_s(combined_m);
    std::vector<Fr> combined_t(combined_m);

    size_t cursor = 0;
    for (size_t i = 0; i < inputs.size(); ++i) {
      auto& input = *inputs[i];
      auto& com_pub = com_pubs[i];
      auto& com_sec = com_secs[i];
      for (size_t j = 0; j < (size_t)input.m; ++j) {
        combined_x[cursor] = std::move(input.x[j]);
        combined_y[cursor] = std::move(input.y[j]);
        combined_z[cursor] = std::move(input.z[j]);
        combined_a[cursor] = com_pub.a[j];
        combined_b[cursor] = com_pub.b[j];
        combined_c[cursor] = com_pub.c[j];
        combined_r[cursor] = com_sec.r[j];
        combined_s[cursor] = com_sec.s[j];
        combined_t[cursor] = com_sec.t[j];
        ++cursor;
      }
    }

    typename Sec43::ProveInput input_43(
        std::move(combined_x), std::move(combined_y), std::move(combined_z),
        get_g, get_g, get_g);

    typename Sec43::CommitmentPub com_pub_43;
    typename Sec43::CommitmentSec com_sec_43;
    com_pub_43.a = std::move(combined_a);
    com_pub_43.b = std::move(combined_b);
    com_pub_43.c = std::move(combined_c);
    com_sec_43.r = std::move(combined_r);
    com_sec_43.s = std::move(combined_s);
    com_sec_43.t = std::move(combined_t);
    BaseR1cs::DebugCheckHpCom(combined_m, input_43, com_pub_43, com_sec_43);

    Sec43::Prove(proof, seed, std::move(input_43), std::move(com_pub_43),
                 std::move(com_sec_43));
  }

  static bool Verify(Proof const& proof, h256_t seed,
                     std::vector<VerifyInput*>&& inputs) {
    Tick tick(__FN__);

    for (auto const& input : inputs) {
      CHECK(input->Check(), input->unique_tag);
    }

    auto const& get_g = inputs[0]->get_g;
    for (auto const& input : inputs) {
      CHECK(input->get_g(0) == get_g(0), input->unique_tag);
    }

    std::sort(inputs.begin(), inputs.end(), [](auto const& a, auto const& b) {
      return a->unique_tag < b->unique_tag;
    });

    UpdateSeed(seed, inputs);
    std::cout << __FN__ << " " << misc::HexToStr(seed) << "\n";
    std::vector<typename Sec43::CommitmentPub> com_pubs(inputs.size());

    auto pf = [&inputs, &com_pubs](int64_t i) {
      auto const& input = *inputs[i];
      auto& com_pub = com_pubs[i];
      BaseR1cs::BuildHpCom(input.m, input.n, input.com_w, input.constraints(),
                           input.get_g, com_pub);
    };
    parallel::For(inputs.size(), pf);

    std::vector<size_t> mn;
    for (auto const& i : inputs) {
      mn.resize(mn.size() + i->m, i->n);
    }
    std::vector<G1> combined_a(mn.size());
    std::vector<G1> combined_b(mn.size());
    std::vector<G1> combined_c(mn.size());

    size_t cursor = 0;
    for (size_t i = 0; i < inputs.size(); ++i) {
      auto& input = *inputs[i];
      auto& com_pub = com_pubs[i];
      for (size_t j = 0; j < (size_t)input.m; ++j) {
        combined_a[cursor] = com_pub.a[j];
        combined_b[cursor] = com_pub.b[j];
        combined_c[cursor] = com_pub.c[j];
        ++cursor;
      }
    }

    typename Sec43::CommitmentPub com_pub_43;
    com_pub_43.a = std::move(combined_a);
    com_pub_43.b = std::move(combined_b);
    com_pub_43.c = std::move(combined_c);

    typename Sec43::VerifyInput input_43(mn, com_pub_43, get_g, get_g, get_g);
    return Sec43::Verify(proof, seed, input_43);
  }

  // see match.h
  static bool Test() { return true; }

 private:
};

}  // namespace clink