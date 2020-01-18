#pragma once

#include "./mimc5.h"
#include "./mimc5_gadget.h"
#include "./misc.h"
#include "./types.h"
#include "./sha256c.h"
#include "./sha256c_gadget.h"
#include "utils/fst.h"

// Verifiable Random Sequence Prover

namespace vrs {

class Prover {
 public:
  Prover(PublicInput const& public_input, SecretInput const& secret_input,
         std::vector<G1> cached_var_coms, std::vector<Fr> cached_var_coms_r)
      : public_input_(public_input),
        secret_input_(secret_input),
        cached_var_coms_(std::move(cached_var_coms)),
        cached_var_coms_r_(std::move(cached_var_coms_r)) {
    auto count = public_input_.count;
    pds_sigma_g_ = PcComputeSigmaG(g_offset_, count);
    v_.resize(count);
    values_.resize(count);
    pb_.reset(new libsnark::protoboard<Fr>);
    gadget_.reset(new Mimc5Gadget(*pb_, "Mimc5Gadget"));
    pb_->set_input_sizes(primary_input_size_);  // var_plain is public statement

    assert(cached_var_coms_.size() == cached_var_coms_r_.size());
  }

  void DebugCheckInput(CheckInput const& check_input) {
#ifdef _DEBUG
    assert(check_input.v == v_);
    assert(vw_ == check_input.vw);
#else
    (void)check_input;
#endif
  }

  PublicInput const& public_input() const { return public_input_; }

  SecretInput const& secret_input() const { return secret_input_; }

  std::vector<Fr> const& v() const { return v_; }

  Fr const& vw() const { return vw_; }

  std::vector<G1> const& var_coms() const { return var_coms_; }

  // We need to compute the plain and v even have cache because we can not cache
  // the v (too large)
  void Evaluate() {
    // Tick tick(__FUNCTION__);
    auto count = public_input_.count;
    for (int64_t i = 0; i < count; ++i) {
      auto plain = public_input_.get_p(i);
      gadget_->Assign(plain, secret_input_.key);
      assert(pb_->is_satisfied());
      v_[i] = pb_->val(gadget_->result());
      values_[i] = pb_->full_variable_assignment();
    }
  }

  void Prove(h256_t seed, std::function<Fr(int64_t)> get_w, Proof& proof,
             ProveOutput& output) {
    // Tick tick(__FUNCTION__);
    BuildVarComs();

    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, var_coms_);
    hash.Final(seed.data());

    HpProve(seed, proof.proof_hp);
    IpProve(seed, std::move(get_w), proof.proof_ip);
    proof.var_coms = var_coms_;
    proof.com_vw = com_vw_;

    output.g = pds_sigma_g_;
    output.h = GetPcBase().h();
    output.key_com = var_coms_[1];  // var_coms_[0]:plain, var_coms_[1]:key
  }

 private:
  int64_t num_variables() const { return (int64_t)pb_->num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb_->num_constraints(); }

  void BuildVarComs() {
    if ((int64_t)cached_var_coms_.size() == num_variables()) {
      var_coms_ = std::move(cached_var_coms_);
      var_coms_r_ = std::move(cached_var_coms_r_);
      for (int64_t i = 0; i < primary_input_size_; ++i) {
        assert(var_coms_r_[i] == FrZero());
      }
      assert(var_coms_r_[primary_input_size_] == secret_input_.key_com_r);
    } else {
      assert(cached_var_coms_.empty());
      ComputeVarComs();
    }
  }

  void ComputeVarComs() {
    // Tick tick(__FUNCTION__);
    auto count = public_input_.count;
    var_coms_.resize(num_variables());
    var_coms_r_.resize(var_coms_.size());

    auto parallel_f = [this, count](int64_t i) {
      std::vector<Fr> data(count);
      auto& var_com = var_coms_[i];
      auto& var_com_r = var_coms_r_[i];
      for (int64_t j = 0; j < count; ++j) {
        data[j] = values_[j][i];
      }
      // var_coms_[0]:plain, var_coms_[1]:key
      if (i < primary_input_size_) {
        var_com_r = FrZero();  // public val
      } else if (i == primary_input_size_) {
        var_com_r = secret_input_.key_com_r;
      } else {
        var_com_r = FrRand();
      }
      var_com = PcComputeCommitmentG(g_offset_, data, var_com_r);
    };
    parallel::For((int64_t)var_coms_.size(), parallel_f);
  }

  // <A,X>=constraint.a, <B,X>=constraint.b, <C,X>=constraint.c
  groth09::sec43b::ProverInput BuildHpInput() {
    // Tick tick(__FUNCTION__);
    auto m = num_constraints();
    auto n = public_input_.count;
    std::vector<std::vector<Fr>> x(m);
    std::vector<std::vector<Fr>> y(m);
    std::vector<std::vector<Fr>> z(m);
    for (auto& i : x) i.resize(n);
    for (auto& i : y) i.resize(n);
    for (auto& i : z) i.resize(n);
    // std::cout << "sec43b: " << m << "*" << n << "\n";

    auto constraint_system = pb_->get_constraint_system();
    auto const& constraints = constraint_system.constraints;

    auto parallel_f = [this, &constraints, &x, &y, &z, m](int64_t j) {
      auto const& values = values_[j];
      for (int64_t i = 0; i < m; ++i) {
        auto const& constraint = constraints[i];
        x[i][j] = constraint.a.evaluate(values);
        y[i][j] = constraint.b.evaluate(values);
        z[i][j] = constraint.c.evaluate(values);
        assert(z[i][j] == x[i][j] * y[i][j]);
      }
    };
    parallel::For(n, parallel_f);

    // now we do not need values_
    values_.clear();
    values_.shrink_to_fit();

    return groth09::sec43b::ProverInput(std::move(x), std::move(y),
                                        std::move(z), g_offset_, g_offset_,
                                        g_offset_);
  }

  void BuildHpCom(groth09::sec43b::CommitmentPub& com_pub,
                  groth09::sec43b::CommitmentSec& com_sec) {
    // Tick tick(__FUNCTION__);
    auto m = num_constraints();
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);
    com_sec.r.resize(m);
    FrZero(com_sec.r);
    com_sec.s.resize(m);
    FrZero(com_sec.s);
    com_sec.t.resize(m);
    FrZero(com_sec.t);

    auto constraint_system = pb_->get_constraint_system();
    auto const& constraints = constraint_system.constraints;

    auto parallel_f = [this, &com_pub, &com_sec,
                       &constraints](int64_t i) mutable {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto& com_sec_r = com_sec.r[i];
      auto& com_sec_s = com_sec.s[i];
      auto& com_sec_t = com_sec.t[i];
      auto const& constraint = constraints[i];
      BuildHpCom(constraint.a, com_pub_a, com_sec_r, pds_sigma_g_);
      BuildHpCom(constraint.b, com_pub_b, com_sec_s, pds_sigma_g_);
      BuildHpCom(constraint.c, com_pub_c, com_sec_t, pds_sigma_g_);
    };
    parallel::For(m, parallel_f);
  }

  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  void BuildHpCom(libsnark::linear_combination<Fr> const& lc, G1& com_pub,
                  Fr& com_sec, G1 const& sigma_g) {
    // Tick tick(__FUNCTION__);
    for (auto const& term : lc.terms) {
      if (term.index == 0) {
        com_pub += sigma_g * term.coeff;
      } else {
        com_pub += var_coms_[term.index - 1] * term.coeff;
        com_sec += var_coms_r_[term.index - 1] * term.coeff;
      }
    }
  }

  void DebugCheckHpCom(groth09::sec43b::ProverInput const& input,
                       groth09::sec43b::CommitmentPub const& com_pub,
                       groth09::sec43b::CommitmentSec const& com_sec) {
#ifdef _DEBUG
    Tick tick(__FUNCTION__);
    auto m = num_constraints();
    for (int64_t i = 0; i < m; ++i) {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto& com_sec_r = com_sec.r[i];
      auto& com_sec_s = com_sec.s[i];
      auto& com_sec_t = com_sec.t[i];

      auto const& xi = input.x[i];
      G1 check_com_pub_a = PcComputeCommitmentG(g_offset_, xi, com_sec_r);
      assert(check_com_pub_a == com_pub_a);

      auto const& yi = input.y[i];
      G1 check_com_pub_b = PcComputeCommitmentG(g_offset_, yi, com_sec_s);
      assert(check_com_pub_b == com_pub_b);

      auto const& zi = input.z[i];
      G1 check_com_pub_c = PcComputeCommitmentG(g_offset_, zi, com_sec_t);
      assert(check_com_pub_c == com_pub_c);
    }
#else
    (void)input;
    (void)com_pub;
    (void)com_sec;
#endif
  }

  void HpProve(h256_t const& seed, groth09::sec43b::RomProof& rom_proof) {
    // Tick tick(__FUNCTION__);
    auto input = BuildHpInput();
    groth09::sec43b::CommitmentPub com_pub;
    groth09::sec43b::CommitmentSec com_sec;
    BuildHpCom(com_pub, com_sec);
    DebugCheckHpCom(input, com_pub, com_sec);

    gadget_.reset();
    pb_.reset();

    groth09::sec43b::AlignData(input, com_pub, com_sec);
    groth09::sec43b::RomProve(rom_proof, seed, std::move(input),
                              std::move(com_pub), std::move(com_sec));
  }

  void IpProve(h256_t const& seed, std::function<Fr(int64_t)> get_w,
               hyrax::a2::RomProof& rom_proof) {
    // Tick tick(__FUNCTION__);
    std::vector<Fr> input_w(public_input_.count);
    for (int64_t i = 0; i < public_input_.count; ++i) {
      input_w[i] = get_w(i);
    }
    vw_ = InnerProduct(v_, input_w);
    int64_t x_g_offset = g_offset_;
    int64_t y_g_offset = -1;
    hyrax::a2::ProverInput input(v_, input_w, vw_, x_g_offset, y_g_offset);
    hyrax::a2::CommitmentPub com_pub;
    hyrax::a2::CommitmentSec com_sec;
    com_sec.r_xi = var_coms_r_.back();
    com_sec.r_tau = secret_input_.vw_com_r;
    com_pub.xi = var_coms_.back();  // the last var is the result of the mimc
    com_pub.tau = PcComputeCommitmentG(y_g_offset, input.y, com_sec.r_tau);
    com_vw_ = com_pub.tau;
    hyrax::a2::RomProve(rom_proof, seed, std::move(input), std::move(com_pub),
                        std::move(com_sec));
  }

 private:
  int64_t const primary_input_size_ = 1;
  PublicInput public_input_;
  SecretInput secret_input_;
  std::vector<G1> cached_var_coms_;
  std::vector<Fr> cached_var_coms_r_;
  std::unique_ptr<libsnark::protoboard<Fr>> pb_;
  std::unique_ptr<Mimc5Gadget> gadget_;
  std::vector<std::vector<Fr>> values_;
  std::vector<G1> var_coms_;
  std::vector<Fr> var_coms_r_;
  std::vector<Fr> v_;
  Fr vw_;
  G1 com_vw_;
  G1 pds_sigma_g_;
  int64_t const g_offset_ = 0;  // hardcode 0 because the cache use 0
};

}  // namespace vrs
