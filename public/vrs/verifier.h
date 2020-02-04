#pragma once

#include "./mimc5_gadget.h"
#include "./misc.h"
#include "./prover.h"
#include "./types.h"

// Verifiable Random Sequence Verifier

namespace vrs {

template <typename Scheme, typename Policy>
class Verifier {
 public:
  using HyraxA = typename Policy::HyraxA;
  using Sec43 = typename Policy::Sec43;

  Verifier(PublicInput const& public_input) : public_input_(public_input) {
    pds_sigma_g_ = PcComputeSigmaG(g_offset_, public_input_.count);
    scheme_.reset(new Scheme);
  }

  bool Verify(h256_t seed, std::function<Fr(int64_t)> get_w,
              Proof<Policy> const& proof, VerifyOutput& output) {
    // Tick tick(__FUNCTION__);
    auto n = public_input_.count;
    auto m = num_constraints();
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, proof.var_coms);
    hash.Final(seed.data());

    std::array<parallel::Task, 3> tasks;
    // check com_plain
    bool ret_com_plain = false;
    tasks[0] = [this, &proof, &ret_com_plain, n]() {
      G1 const& com_plain = proof.var_coms[0];
      Fr com_plain_r = FrZero();  // public val
      std::vector<Fr> data(n);
      for (int64_t i = 0; i < n; ++i) {
        data[i] = public_input_.get_p(i);
      }
      int64_t x_g_offset = 0;  // hardcode 0, because prover always use 0
      auto check_value = PcComputeCommitmentG(x_g_offset, data, com_plain_r);
      ret_com_plain = check_value == com_plain;
    };

    // check hadamard product
    bool ret_hp = false;
    tasks[1] = [this, &proof, &ret_hp, &seed, m, n]() {
      Sec43::CommitmentPub com_pub_hp;
      BuildHpCom(proof, com_pub_hp);
      com_pub_hp.Align();
      int64_t x_g_offset = 0;
      int64_t y_g_offset = 0;
      int64_t z_g_offset = 0;
      Sec43::VerifierInput input_hp(m, n, com_pub_hp, x_g_offset, y_g_offset,
                                    z_g_offset);
      ret_hp = Sec43::Verify(proof.proof_hp, seed, input_hp);
    };

    // check inner product
    bool ret_ip = false;
    tasks[2] = [this, &proof, &ret_ip, n, &get_w, &seed]() {
      std::vector<Fr> input_w(n);
      for (int64_t i = 0; i < (int64_t)input_w.size(); ++i) {
        input_w[i] = get_w(i);
      }
      HyraxA::CommitmentPub com_pub_ip;
      BuildIpCom(proof, com_pub_ip);
      int64_t x_g_offset = 0;
      int64_t y_g_offset = -1;
      HyraxA::VerifierInput input_ip(input_w, com_pub_ip, x_g_offset,
                                     y_g_offset);
      ret_ip = HyraxA::Verify(proof.proof_ip, seed, input_ip);
    };

    parallel::Invoke(tasks);

    if (!ret_com_plain || !ret_hp || !ret_ip) {
      std::cout << "ret_com_plain: " << ret_com_plain << ", ret_hp: " << ret_hp
                << ", ret_ip: " << ret_ip << "\n";
      assert(false);
      return false;
    }

    output.g = pds_sigma_g_;
    output.h = GetPcBase().h();
    output.key_com = proof.var_coms[1];
    return true;
  }

 private:
  int64_t num_variables() const { return scheme_->num_variables(); }

  int64_t num_constraints() const { return scheme_->num_constraints(); }

  void BuildIpCom(Proof<Policy> const& proof,
                  typename HyraxA::CommitmentPub& com_pub) {
    com_pub.xi = proof.var_coms.back();
    com_pub.tau = proof.com_vw;
  }

  void BuildHpCom(Proof<Policy> const& proof,
                  typename Sec43::CommitmentPub& com_pub) {
    // Tick tick(__FUNCTION__);
    auto m = num_constraints();
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);

    auto const& constraint_system = scheme_->constraint_system;
    auto const& constraints = constraint_system.constraints;

    auto parallel_f = [this, &com_pub, &constraints,
                       &proof](int64_t i) mutable {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto const& constraint = constraints[i];
      BuildHpCom(constraint.a, com_pub_a, pds_sigma_g_, proof);
      BuildHpCom(constraint.b, com_pub_b, pds_sigma_g_, proof);
      BuildHpCom(constraint.c, com_pub_c, pds_sigma_g_, proof);
    };
    parallel::For(m, parallel_f);
  }

  // com(<A,X>) or com(<B,X>) or com(<C,X>)
  void BuildHpCom(libsnark::linear_combination<Fr> const& lc, G1& com_pub,
                  G1 const& sigma_g, Proof<Policy> const& proof) {
    // Tick tick(__FUNCTION__);
    for (auto const& term : lc.terms) {
      if (term.index == 0) {
        com_pub += sigma_g * term.coeff;
      } else {
        com_pub += proof.var_coms[term.index - 1] * term.coeff;
      }
    }
  }

 private:
  PublicInput public_input_;
  std::unique_ptr<Scheme> scheme_;
  G1 pds_sigma_g_;
  int64_t const g_offset_ = 0;
};

}  // namespace vrs