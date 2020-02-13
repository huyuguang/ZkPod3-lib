#pragma once

#include "./match.h"
#include "./pack.h"

namespace clink {

template <typename Policy>
struct MatchPack {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using R1cs = typename clink::ParallelR1cs<Policy>;

  struct Proof {
    typename Match<Policy>::Proof match_proof;
    typename Pack<HyraxA>::Proof pack_proof;
    G1 com_pack_y;
    bool operator==(Proof const& b) const {
      return match_proof == b.match_proof && pack_proof == b.pack_proof &&
             com_pack_y == b.com_pack_y;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("mp.p", ("s", match_proof), ("p", pack_proof),
                         ("y", com_pack_y));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("mp.p", ("s", match_proof), ("p", pack_proof),
                         ("y", com_pack_y));
    }
  };

  struct ProveOutput {
    typename Match<Policy>::Proof match_proof;
    typename Pack<HyraxA>::Proof pack_proof;
    std::vector<Fr> y;
    Fr com_y_r;
    G1 com_pack_y;
    std::vector<Fr> pack_y;
    Fr com_pack_y_r;
    Proof BuildProof() const {
      Proof ret;
      ret.match_proof = match_proof;
      ret.pack_proof = pack_proof;
      ret.com_pack_y = com_pack_y;
      return ret;
    }
  };

  struct ProveInput {
    ProveInput(Fr const& k, std::vector<Fr> const& x, G1 const& com_x,
               Fr const& com_x_r, int64_t x_g_offset, int64_t py_g_offset)
        : k(k),
          x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          n((int64_t)x.size()),
          x_g_offset(x_g_offset),
          py_g_offset(py_g_offset) {
#ifdef _DEBUG
      assert(com_x == PcComputeCommitmentG(x_g_offset, x, com_x_r, true));
#endif
    }
    Fr const& k;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    int64_t const n;
    int64_t const x_g_offset;
    int64_t const py_g_offset;
  };

  static void Prove(ProveOutput& output, h256_t seed, ProveInput const& input) {
    // prove y[i]=x[i] == k, i=[0,s)
    int64_t n = input.n;
    output.y.resize(n);
    for (int64_t i = 0; i < n; ++i) {
      output.y[i] = input.x[i] == input.k ? 1 : 0;
    }

    output.com_y_r = FrRand();
    auto com_y =
        PcComputeCommitmentG(input.x_g_offset, output.y, output.com_y_r, true);

    typename Match<Policy>::ProveInput m_input(
        input.k, input.x, input.com_x, input.com_x_r, output.y, com_y,
        output.com_y_r, input.x_g_offset);
    auto& match_proof = output.match_proof;
    Match<Policy>::Prove(match_proof, seed, m_input);
    assert(match_proof.com_w.back() == com_y);

    // Pack<HyraxA> y to pack_y
    output.pack_y = FrBitsToFrs(output.y);
    output.com_pack_y_r = FrRand();
    output.com_pack_y = PcComputeCommitmentG(input.py_g_offset, output.pack_y,
                                             output.com_pack_y_r);
    auto& pack_proof = output.pack_proof;
    typename Pack<HyraxA>::ProveInput p_input(
        output.y, com_y, output.com_y_r, input.x_g_offset, output.pack_y,
        output.com_pack_y, output.com_pack_y_r, input.py_g_offset);
    Pack<HyraxA>::Prove(pack_proof, seed, p_input);
  }

  struct VerifyInput {
    VerifyInput(int64_t n, Fr const& k, int64_t const x_g_offset,
                int64_t const py_g_offset)
        : n(n), k(k), x_g_offset(x_g_offset), py_g_offset(py_g_offset) {}
    int64_t n;
    Fr const& k;
    int64_t const x_g_offset;
    int64_t const py_g_offset;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    typename Match<Policy>::VerifyInput m_input(input.n, input.k,
                                                input.x_g_offset);
    if (!Match<Policy>::Verify(proof.match_proof, seed, m_input)) return false;

    G1 const& com_y = proof.match_proof.com_w.back();
    typename Pack<HyraxA>::VerifyInput p_input(
        input.n, com_y, input.x_g_offset, proof.com_pack_y, input.py_g_offset);
    return Pack<HyraxA>::Verify(proof.pack_proof, seed, p_input);
  }

  static bool Test(int64_t n, std::string const& k);
};

template <typename Policy>
bool MatchPack<Policy>::Test(int64_t n, std::string const& k) {
  if (k.size() > 31) {
    std::cout << "invalid parameter: k.size() must <= 31.\n";
    return false;
  }
  if (n >= PcBase::kGSize / 2) {
    std::cout << "invalid parameter: n must < " << PcBase::kGSize / 2 << "\n";
    return false;
  }

  auto seed = misc::RandH256();
  int64_t x_g_offset = 0;
  int64_t py_g_offset = 20;
  Fr fr_k = PackStrToFr(k.c_str());
  std::vector<Fr> x(n);
  FrRand(x);
  x[rand() % n] = fr_k;
  x[rand() % n] = fr_k;
  auto com_x_r = FrRand();
  auto com_x = PcComputeCommitmentG(x_g_offset, x, com_x_r);

  ProveInput prove_input(fr_k, x, com_x, com_x_r, x_g_offset, py_g_offset);
  ProveOutput output;
  Prove(output, seed, prove_input);

  auto proof = output.BuildProof();
  if (proof.match_proof.com_x() != com_x) {
    assert(false);
    return false;
  }

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

  VerifyInput verify_input(n, fr_k, x_g_offset, py_g_offset);
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink