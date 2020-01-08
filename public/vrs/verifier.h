#pragma once

#include "./mimc5_gadget.h"
#include "./misc.h"
#include "./prover.h"
#include "./types.h"

// Verifiable Random Sequence Verifier

namespace vrs {

class Verifier {
 public:
  Verifier(PublicInput const& public_input) : public_input_(public_input) {
    pds_sigma_g_ = PcComputeSigmaG(public_input_.count);
    Evaluate();
  }

  bool Verify(h256_t seed, std::function<Fr(int64_t)> get_w,
              Proof const& proof, VerifyOutput& output) {
    // Tick tick(__FUNCTION__);
    auto count = public_input_.count;

    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, proof.var_coms);
    hash.Final(seed.data());

    std::array<parallel::Task, 3> tasks;
    // check com_plain
    bool ret_com_plain = false;
    tasks[0] = [this, &proof, &ret_com_plain, count]() {
      G1 const& com_plain = proof.var_coms[0];
      Fr com_plain_r = FrZero(); // public val
      std::vector<Fr> data(count);
      for (int64_t i = 0; i < count; ++i) {
        data[i] = public_input_.get_p(i);
      }
      auto check_value = PcComputeCommitment(data, com_plain_r);
      ret_com_plain = check_value == com_plain;
    };

    // check hadamard product
    bool ret_hp = false;
    tasks[1] = [this, &proof, &ret_hp, &seed]() {
      groth09::sec43::CommitmentPub com_pub_hp;
      BuildHpCom(proof, com_pub_hp);
      com_pub_hp.Align();
      groth09::sec43::VerifierInput input_hp(com_pub_hp);
      ret_hp = groth09::sec43::RomVerify(proof.proof_hp, seed, input_hp);
    };

    // check inner product
    bool ret_ip = false;
    tasks[2] = [this, &proof, &ret_ip, count, &get_w, &seed]() {
      std::vector<Fr> input_w(count);
      for (int64_t i = 0; i < (int64_t)input_w.size(); ++i) {
        input_w[i] = get_w(i);
      }
      hyrax::a2::CommitmentPub com_pub_ip;
      BuildIpCom(proof, com_pub_ip);
      hyrax::a2::VerifierInput input_ip(input_w, com_pub_ip);
      ret_ip = hyrax::a2::RomVerify(proof.proof_ip, seed, input_ip);
    };

    parallel::Invoke(tasks);

    if (!ret_com_plain || !ret_hp || !ret_ip) {
      std::cout << "ret_com_plain: " << ret_com_plain << ", ret_hp: " << ret_hp
                << ", ret_ip" << ret_ip << "\n";
      assert(false);
      return false;
    }

    output.g = pds_sigma_g_;
    output.h = GetPcBase().h();
    output.key_com = proof.var_coms[1];
    return true;
  }

 private:
  int64_t num_variables() const { return (int64_t)pb_.num_variables(); }

  int64_t num_constraints() const { return (int64_t)pb_.num_constraints(); }

  void Evaluate() {
    // Tick tick(__FUNCTION__);
    gadget_.reset(new Mimc5Gadget(pb_));
    pb_.set_input_sizes(primary_input_size_);  // var_plain is public statement
  }

  void BuildIpCom(Proof const& proof, hyrax::a2::CommitmentPub& com_pub) {
    com_pub.xi = proof.var_coms.back();
    com_pub.tau = proof.com_vw;
  }

  void BuildHpCom(Proof const& proof, groth09::sec43::CommitmentPub& com_pub) {
    // Tick tick(__FUNCTION__);
    auto m = num_constraints();
    com_pub.a.resize(m);
    G1Zero(com_pub.a);
    com_pub.b.resize(m);
    G1Zero(com_pub.b);
    com_pub.c.resize(m);
    G1Zero(com_pub.c);

    auto constraint_system = pb_.get_constraint_system();
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