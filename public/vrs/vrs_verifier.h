#pragma once

#include "vrs_mimc5_gadget.h"
#include "vrs_misc.h"
#include "vrs_prover.h"
#include "vrs_types.h"

// Verifiable Random Sequence Verifier

namespace vrs {

class Verifier {
 public:
  Verifier(PublicInput const& public_input)
      : public_input_(public_input) {
    pds_sigma_g_ = ComputePdsSigmaG(public_input_.count);    
    Evaluate();
  }

  bool Verify(h256_t const& rom_seed, std::function<Fr(int64_t)> get_w,
              Proof const& proof, VerifyOutput& output) {
    Tick tick(__FUNCTION__);
    using groth09::details::ComputeCommitment;
    auto count = public_input_.count;

    // check com_plain
    G1 const& com_plain = proof.var_coms[0];
    Fr com_plain_r = FrZero();
    std::vector<Fr> data(count);
    for (int64_t i = 0; i < count; ++i) {
      data[i] = public_input_.get_p(i);      
    }
    auto check_value = ComputeCommitment(data, com_plain_r);
    if (check_value != com_plain) {
      std::cout << "check_value != com_plain\n";
      assert(false);
      return false;
    }

    // check hadamard product
    groth09::sec43::CommitmentPub com_pub_hp;
    BuildHpCom(proof, com_pub_hp);
    com_pub_hp.Align();
    groth09::sec43::VerifierInput input_hp(com_pub_hp);
    if (!groth09::sec43::RomVerify(proof.proof_hp, rom_seed,
                                   input_hp)) {
      std::cout << "groth09::sec43::RomVerify failed\n";
      assert(false);
      return false;
    }

    // check inner product
    std::vector<Fr> input_w(public_input_.count);
    for (int64_t i = 0; i < (int64_t)input_w.size(); ++i) {
      input_w[i] = get_w(i);
    }
    hyrax::a2::CommitmentPub com_pub_ip;
    BuildIpCom(proof, com_pub_ip);
    hyrax::a2::VerifierInput input_ip(input_w, com_pub_ip);
    if (!hyrax::a2::RomVerify(proof.proof_ip, rom_seed, input_ip)) {
      std::cout << "hyrax::a2::RomVerify failed\n";
      assert(false);
      return false;
    }

    output.g = pds_sigma_g_;
    output.h = GetPdsPub().h();
    output.key_com = proof.var_coms[1];
    return true;
  }

 private:
  int64_t num_variables() const { return (int64_t)pb_.num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb_.num_constraints(); }

  void Evaluate() {
    //Tick tick(__FUNCTION__);
    gadget_.reset(new Mimc5Gadget(pb_));
    pb_.set_input_sizes(primary_input_size_);  // var_plain is public statement
  }

  void BuildIpCom(Proof const& proof, hyrax::a2::CommitmentPub& com_pub) {
    com_pub.xi = proof.var_coms.back();
    com_pub.tau = proof.com_vw;
  }

  void BuildHpCom(Proof const& proof, groth09::sec43::CommitmentPub& com_pub) {
    Tick tick(__FUNCTION__);
    auto m = num_constraints();
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);

    auto constraint_system = pb_.get_constraint_system();
    auto const& constraints = constraint_system.constraints;    

//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
    for (int64_t i = 0; i < m; ++i) {
      auto& com_pub_a = com_pub.a[i];
      auto& com_pub_b = com_pub.b[i];
      auto& com_pub_c = com_pub.c[i];
      auto const& constraint = constraints[i];
      BuildHpCom(constraint.a, com_pub_a, pds_sigma_g_, proof);
      BuildHpCom(constraint.b, com_pub_b, pds_sigma_g_, proof);
      BuildHpCom(constraint.c, com_pub_c, pds_sigma_g_, proof);
    }
  }

  void BuildHpCom(libsnark::linear_combination<Fr> const& lc, G1& com_pub,
                  G1 const& sigma_g, Proof const& proof) {
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
  int64_t const primary_input_size_ = 1;
  PublicInput public_input_;
  libsnark::protoboard<Fr> pb_;
  std::shared_ptr<Mimc5Gadget> gadget_;
  G1 pds_sigma_g_;
};

}  // namespace vrs