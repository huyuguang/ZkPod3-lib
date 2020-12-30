#pragma once

#include "./pack.h"
#include "./substr.h"

namespace clink {

template <typename Policy>
struct SubstrPack {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using R1cs = typename clink::ParallelR1cs<Policy>;

  struct Proof {
    typename Substr<Policy>::Proof substr_proof;
    typename Pack<HyraxA>::Proof pack_proof;
    G1 com_pack_y;
    bool operator==(Proof const& b) const {
      return substr_proof == b.substr_proof && pack_proof == b.pack_proof &&
             com_pack_y == b.com_pack_y;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("sp.p", ("s", substr_proof), ("p", pack_proof),
                         ("y", com_pack_y));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("sp.p", ("s", substr_proof), ("p", pack_proof),
                         ("y", com_pack_y));
    }
  };

  struct ProveOutput {
    typename Substr<Policy>::Proof substr_proof;
    typename Pack<HyraxA>::Proof pack_proof;
    std::vector<Fr> y;
    Fr com_y_r;
    G1 com_pack_y;
    std::vector<Fr> pack_y;
    Fr com_pack_y_r;
    Proof BuildProof() const {
      Proof ret;
      ret.substr_proof = substr_proof;
      ret.pack_proof = pack_proof;
      ret.com_pack_y = com_pack_y;
      return ret;
    }
  };

  struct ProveInput {
    ProveInput(std::string const& k, std::vector<Fr> const& x, G1 const& com_x,
               Fr const& com_x_r, GetRefG1 const& get_gx,
               GetRefG1 const& get_gpy)
        : k(k),
          x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          n((int64_t)x.size()),
          get_gx(get_gx),
          get_gpy(get_gpy) {
#ifdef _DEBUG
      assert(com_x == pc::ComputeCom(get_gx, x, com_x_r, true));
#endif
    }
    std::string const& k;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    int64_t const n;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gpy;
  };

  static void Prove(ProveOutput& output, h256_t seed, ProveInput const& input) {
    // prove y[i]=substr(x[i],k), i=[0,s)
    int64_t n = input.n;
    output.y.resize(n);
    std::vector<std::string> x_str(n);
    for (int64_t i = 0; i < n; ++i) {
      auto str = UnPackStrFromFr(input.x[i]);
      output.y[i] = str.find(input.k) != std::string::npos ? 1 : 0;
    }

    output.com_y_r = FrRand();
    auto com_y = pc::ComputeCom(input.get_gx, output.y, output.com_y_r, true);

    typename Substr<Policy>::ProveInput s_input(input.k, input.x, input.com_x,
                                                input.com_x_r, output.y, com_y,
                                                output.com_y_r, input.get_gx);
    auto& substr_proof = output.substr_proof;
    Substr<Policy>::Prove(substr_proof, seed, s_input);
    assert(substr_proof.com_w.back() == com_y);

    // pack y to pack_y
    output.pack_y = FrBitsToFrs(output.y);
    output.com_pack_y_r = FrRand();
    output.com_pack_y =
        pc::ComputeCom(input.get_gpy, output.pack_y, output.com_pack_y_r);
    auto& pack_proof = output.pack_proof;
    typename Pack<HyraxA>::ProveInput p_input(
        output.y, com_y, output.com_y_r, input.get_gx, output.pack_y,
        output.com_pack_y, output.com_pack_y_r, input.get_gpy);

    clink::Pack<HyraxA>::Prove(pack_proof, seed, p_input);
  }

  struct VerifyInput {
    VerifyInput(int64_t n, std::string const& k, GetRefG1 const& get_gx,
                GetRefG1 const& get_gpy)
        : n(n), k(k), get_gx(get_gx), get_gpy(get_gpy) {}
    int64_t n;
    std::string const& k;
    GetRefG1 const& get_gx;
    GetRefG1 const& get_gpy;
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    typename Substr<Policy>::VerifyInput s_input(input.n, input.k,
                                                 input.get_gx);
    if (!Substr<Policy>::Verify(proof.substr_proof, seed, s_input))
      return false;

    G1 const& com_y = proof.substr_proof.com_w.back();
    typename Pack<HyraxA>::VerifyInput p_input(input.n, com_y, input.get_gx,
                                               proof.com_pack_y, input.get_gpy);
    return Pack<HyraxA>::Verify(proof.pack_proof, seed, p_input);
  }

  static bool Test(int64_t n, std::string const& k);
};

template <typename Policy>
bool SubstrPack<Policy>::Test(int64_t n, std::string const& k) {
  if (k.size() > 31) {
    std::cout << "invalid parameter: k.size() must <= 31.\n";
    return false;
  }
  if (n >= pc::Base::GSize() / 2) {
    std::cout << "invalid parameter: n must < " << pc::Base::GSize() / 2
              << "\n";
    return false;
  }

  auto seed = misc::RandH256();
  int64_t x_g_offset = 0;    // any value
  int64_t py_g_offset = 20;  // any value
  GetRefG1 get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  GetRefG1 get_gpy = [py_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[py_g_offset + i];
  };

  std::vector<Fr> x(n);
  for (int64_t i = 0; i < n; ++i) {
    x[i] = PackStrToFr(misc::RandString(31).c_str());
  }

  if (k.size() == 31) {
    x[rand() % n] = PackStrToFr(k.c_str());
    x[rand() % n] = PackStrToFr(k.c_str());
  } else {
    x[rand() % n] = PackStrToFr((k + 'a').c_str());
    x[rand() % n] = PackStrToFr((std::string("b") + k).c_str());
  }

  auto com_x_r = FrRand();
  auto com_x = pc::ComputeCom(get_gx, x, com_x_r);

  ProveInput prove_input(k, x, com_x, com_x_r, get_gx, get_gpy);
  ProveOutput output;
  Prove(output, seed, prove_input);

  auto proof = output.BuildProof();
  if (proof.substr_proof.com_x() != com_x) {
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

  VerifyInput verify_input(n, k, get_gx, get_gpy);
  bool success = Verify(proof, seed, verify_input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink