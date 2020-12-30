#pragma once

#include "../details.h"
#include "../parallel_r1cs.h"
#include "./vrs_pub.h"
#include "circuit/mimc5_gadget.h"
#include "circuit/sha256c_gadget.h"

namespace clink {

template <typename Scheme, typename Policy>
struct VrsBasic {
  using R1cs = typename clink::ParallelR1cs<Policy>;
  using HyraxA = typename Policy::HyraxA;

  typedef std::function<Fr const&(int64_t)> GetP;
  typedef std::function<Fr const&(int64_t)> GetW;

  // since cache always use 0 as g_offset, here also use 0 as g_offset

  struct ProveInput {
    ProveInput(int64_t n, Fr const& k, Fr const& k_com_r, GetP&& iget_p,
               GetW&& iget_w, Fr const& vw_com_r, G1 const& gvw)
        : k(k),
          k_com_r(k_com_r),
          get_p(std::move(iget_p)),
          get_w(std::move(iget_w)),
          vw_com_r(vw_com_r),
          gvw(gvw),
          scheme(new Scheme),
          n(n),
          m(scheme->num_variables()) {}

    std::vector<std::vector<Fr>> EvaluateAndCommit(
        std::vector<G1>&& icached_var_coms,
        std::vector<Fr>&& icached_var_coms_r) {
      std::vector<std::vector<Fr>> vars(m);
      for (auto& i : vars) {
        i.resize(n);
      }

      libsnark::protoboard<Fr> pb;
      auto gadget = scheme->CreateGadget(pb);
      for (int64_t j = 0; j < n; ++j) {
        Fr const& p = get_p(j);
        gadget->Assign(p, k);
        assert(pb.is_satisfied());
        auto var = pb.full_variable_assignment();
        for (int64_t i = 0; i < m; ++i) {
          vars[i][j] = var[i];
        }
        assert(vars[0][j] == p);
        assert(vars[1][j] == k);
      }

      BuildVarComs(vars, std::move(icached_var_coms),
                   std::move(icached_var_coms_r));

      return vars;
    }

    G1 ComputeSigmaG() const { return pc::ComputeSigmaG(0, n); }

    R1csInfo const& r1cs_info() const { return *scheme->r1cs_info; }

    // hardcode as 0 because cache always use 0
    Fr const& k;
    Fr const& k_com_r;
    GetP get_p;
    GetW get_w;
    Fr const& vw_com_r;
    G1 const& gvw;
    std::unique_ptr<Scheme> scheme;
    int64_t const n;
    int64_t const m;
    std::vector<G1> var_coms;
    std::vector<Fr> var_coms_r;

   private:
    void BuildVarComs(std::vector<std::vector<Fr>> const& vars,
                      std::vector<G1>&& icached_var_coms,
                      std::vector<Fr>&& icached_var_coms_r) {
      auto constexpr kPrimaryInputSize = Scheme::kPrimaryInputSize;

      assert(icached_var_coms.size() == icached_var_coms_r.size());

      if ((int64_t)icached_var_coms.size() == m) {
        assert(CheckCache(vars, icached_var_coms, icached_var_coms_r));
        var_coms = std::move(icached_var_coms);
        var_coms_r = std::move(icached_var_coms_r);
      } else {
        Tick _tick_(__FN__);
        assert(icached_var_coms.empty());
        var_coms.resize(m);
        var_coms_r.resize(m);
        for (int64_t i = 0; i < m; ++i) {
          // var_coms[0]:plain, var_coms[1]:key
          if (i < kPrimaryInputSize) {
            var_coms_r[i] = FrZero();  // public val
          } else if (i == kPrimaryInputSize) {
            var_coms_r[i] = k_com_r;
          } else {
            var_coms_r[i] = FrRand();
          }
        }
        auto parallel_f = [this, &vars](int64_t i) {
          var_coms[i] = pc::ComputeCom(vars[i], var_coms_r[i]);
        };
        parallel::For(m, parallel_f);
      }
    }

    bool CheckCache(std::vector<std::vector<Fr>> const& vars,
                    std::vector<G1> const& cached_var_coms,
                    std::vector<Fr> const& cached_var_coms_r) {
      auto constexpr kPrimaryInputSize = Scheme::kPrimaryInputSize;
      for (int64_t i = 0; i < kPrimaryInputSize; ++i) {
        if (cached_var_coms_r[i] != FrZero()) {
          assert(false);
          return false;
        }
      }
      if (cached_var_coms_r[kPrimaryInputSize] != k_com_r) {
        assert(false);
        return false;
      }

      std::vector<G1> check_cached_var_coms(m);
      auto parallel_f = [&vars, &check_cached_var_coms,
                         &cached_var_coms_r](int64_t i) {
        check_cached_var_coms[i] =
            pc::ComputeCom(vars[i], cached_var_coms_r[i]);
      };
      parallel::For(m, parallel_f);

      if (cached_var_coms != check_cached_var_coms) {
        assert(false);
        return false;
      }
      return true;
    }
  };

  struct Proof {
    std::vector<G1> var_coms;
    G1 vw_com;
    typename R1cs::Proof r1cs_proof;
    typename HyraxA::Proof ip_proof;
    bool operator==(Proof const& right) const {
      return var_coms == right.var_coms && vw_com == right.vw_com &&
             r1cs_proof == right.r1cs_proof && ip_proof == right.ip_proof;
    }
    bool operator!=(Proof const& right) const { return !(*this == right); }
    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", var_coms), ("vw_com", vw_com),
                         ("r1cs_proof", r1cs_proof), ("ip_proof", ip_proof));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("vrs.pf", ("var_coms", var_coms), ("vw_com", vw_com),
                         ("r1cs_proof", r1cs_proof), ("ip_proof", ip_proof));
    }
  };

  struct ProveOutput {
    G1 g;
    G1 h;
    G1 key_com;
  };

  static Fr Prove(Proof& proof, ProveOutput& output, h256_t seed,
                  ProveInput&& input, std::vector<G1>&& icached_var_coms,
                  std::vector<Fr>&& icached_var_coms_r) {
    Tick _tick_(__FN__);
    auto constexpr kPrimaryInputSize = Scheme::kPrimaryInputSize;
    auto vars = input.EvaluateAndCommit(std::move(icached_var_coms),
                                        std::move(icached_var_coms_r));

    std::vector<Fr> v(input.n);
    for (int64_t i = 0; i < input.n; ++i) {
      v[i] = vars.back()[i];
      assert(Scheme::Generate(input.get_p(i), input.k) == v[i]);
    }

    std::vector<Fr> w(input.n);
    for (int64_t i = 0; i < input.n; ++i) {
      w[i] = input.get_w(i);
    }

    Fr vw = InnerProduct(v, w);
    G1 vw_com = pc::ComputeCom(input.gvw, vw, input.vw_com_r);

    // proof v[i] = hash(p[i],key)
    typename R1cs::ProveInput r1cs_input(input.r1cs_info(), "vrs",
                                         std::move(vars), input.var_coms,
                                         input.var_coms_r, pc::kGetRefG1);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));

    // proof vw = inner_product(v, w)
    typename HyraxA::ProveInput ip_input("vrs", v, w, vw, pc::kGetRefG1,
                                         input.gvw);
    typename HyraxA::CommitmentPub ip_com_pub;
    typename HyraxA::CommitmentSec ip_com_sec;
    ip_com_sec.r_xi = input.var_coms_r.back();
    ip_com_sec.r_tau = input.vw_com_r;
    ip_com_pub.xi = input.var_coms.back();  // com(v)
    ip_com_pub.tau = vw_com;                // com(vw)
    HyraxA::Prove(proof.ip_proof, seed, std::move(ip_input),
                  std::move(ip_com_pub), std::move(ip_com_sec));

    proof.var_coms = input.var_coms;
    proof.vw_com = vw_com;

    output.g = input.ComputeSigmaG();
    output.h = pc::PcH();
    output.key_com = input.var_coms[kPrimaryInputSize];
    return vw;
  }

  struct VerifyInput {
    VerifyInput(int64_t n, GetP&& iget_p, GetW&& iget_w, G1 const& gvw)
        : get_p(std::move(iget_p)),
          get_w(std::move(iget_w)),
          gvw(gvw),
          scheme(new Scheme),
          n(n),
          m(scheme->num_variables()) {}

    G1 ComputeSigmaG() const { return pc::ComputeSigmaG(0, n); }

    R1csInfo const& r1cs_info() const { return *scheme->r1cs_info; }

    GetP get_p;
    GetW get_w;
    G1 const& gvw;
    std::unique_ptr<Scheme> scheme;
    int64_t const n;
    int64_t const m;
  };

  struct VerifyOutput {
    G1 g;
    G1 h;
    G1 key_com;
  };

  static bool Verify(VerifyOutput& output, Proof const& proof, h256_t seed,
                     VerifyInput&& input) {
    std::vector<Fr> p(input.n);
    for (int64_t i = 0; i < input.n; ++i) {
      p[i] = input.get_p(i);
    }
    std::vector<Fr> w(input.n);
    for (int64_t i = 0; i < input.n; ++i) {
      w[i] = input.get_w(i);
    }

    // check com_plain
    if (proof.var_coms[0] != pc::ComputeCom(p, FrZero())) {
      assert(false);
      return false;
    }

    // check ip product
    typename HyraxA::CommitmentPub ip_com_pub;
    // the last var is the v
    ip_com_pub.xi = proof.var_coms.back();  // com(v)
    ip_com_pub.tau = proof.vw_com;          // com(vw)
    typename HyraxA::VerifyInput ip_input("vrs", w, ip_com_pub, pc::kGetRefG1,
                                          input.gvw);
    if (!HyraxA::Verify(proof.ip_proof, seed, ip_input)) {
      assert(false);
      return false;
    }

    // check r1cs (hp product)
    std::vector<std::vector<Fr>> public_w{std::move(p)};
    typename R1cs::VerifyInput r1cs_input(input.n, input.r1cs_info(), "vrs",
                                          proof.var_coms, public_w,
                                          pc::kGetRefG1);
    if (!R1cs::Verify(proof.r1cs_proof, seed, r1cs_input)) {
      assert(false);
      return false;
    }

    output.g = input.ComputeSigmaG();
    output.h = pc::PcH();
    output.key_com = proof.var_coms[1];
    return true;
  }

  static bool Test(int64_t n);
};

template <typename Scheme, typename Policy>
bool VrsBasic<Scheme, Policy>::Test(int64_t n) {
  Tick tick(__FN__);
  auto seed = misc::RandH256();
  std::vector<Fr> p(n);
  FrRand(p);
  Fr k = FrRand();
  Fr k_com_r = FrRand();
  std::vector<Fr> w(n);
  FrRand(w);
  Fr vw_com_r = FrRand();
  G1 gvw = pc::PcU();
  std::vector<G1> icached_var_coms;
  std::vector<Fr> icached_var_coms_r;
  ProveInput prove_input(
      n, k, k_com_r, [&p](int64_t i) -> Fr const& { return p[i]; },
      [&w](int64_t i) -> Fr const& { return w[i]; }, vw_com_r, gvw);

  Proof proof;
  ProveOutput prove_output;
  Prove(proof, prove_output, seed, std::move(prove_input),
        std::move(icached_var_coms), std::move(icached_var_coms_r));

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  if (proof != proof2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  VerifyInput verify_input(
      n, [&p](int64_t i) -> Fr const& { return p[i]; },
      [&w](int64_t i) -> Fr const& { return w[i]; }, gvw);
  VerifyOutput verify_output;
  bool success = Verify(verify_output, proof, seed, std::move(verify_input));
  if (success) {
    assert(prove_output.g == verify_output.g);
    assert(prove_output.h == verify_output.h);
    assert(prove_output.key_com == verify_output.key_com);
    success = VrsPub<Scheme>::VerifySecret(prove_output.h, prove_output.g,
                                           prove_output.key_com, k_com_r, k);
  }
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink