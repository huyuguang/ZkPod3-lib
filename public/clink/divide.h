#pragma once

#include "./details.h"
#include "./equal_ip.h"

// x: vector<Fr>, size = n, n%s = 0
// s[i]: vector<Fr>, size = sn, i = [0~n/sn), sn can be 1.
// open com(x), com(s[i])
// prove x = {s[i]}

namespace clink {

// HyraxA: hyrax::A2 or hyrax::A3
template <typename HyraxA>
struct Divide {
  using Proof = typename EqualIp<HyraxA>::Proof;

  struct ProveInput {
    ProveInput(std::vector<Fr> const& x, G1 const& com_x, Fr const& com_x_r,
               int64_t x_g_offset, int64_t sn, std::vector<G1> const& com_s,
               std::vector<Fr> const& com_s_r, int64_t s_g_offset)
        : x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          x_g_offset(x_g_offset),
          sn(sn),
          com_s(com_s),
          com_s_r(com_s_r),
          s_g_offset(s_g_offset) {
      assert(sn > 0);
      assert(x.size() % sn == 0);
      assert(x.size() / sn == com_s.size());
      assert(com_s.size() == com_s_r.size());
    }
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    int64_t const x_g_offset;
    int64_t const sn;
    std::vector<G1> const& com_s;
    std::vector<Fr> const& com_s_r;
    int64_t const s_g_offset;
  };

  static void UpdateSeed(h256_t& seed, G1 const& com_x, int64_t sn,
                         std::vector<G1> const& com_s) {
    // update seed
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, sn);
    HashUpdate(hash, com_x);
    HashUpdate(hash, com_s);
    hash.Final(seed.data());
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input) {
    UpdateSeed(seed, input.com_x, input.sn, input.com_s);

    // d, e
    std::vector<Fr> d(input.com_s.size());
    ComputeFst(seed, "divide:d", d);
    std::vector<Fr> e(input.sn);
    ComputeFst(seed, "divide:e", e);

    // f
    std::vector<Fr> f;
    f.reserve(input.x.size());
    for (size_t i = 0; i < input.com_s.size(); ++i) {
      auto de = e * d[i];
      f.insert(f.end(), de.begin(), de.end());
    }

    // y
    std::vector<Fr> y(input.sn, FrZero());
    auto get_s = [&input](int64_t i) {
      auto begin = input.x.begin() + i * input.sn;
      auto end = begin + input.sn;
      return std::vector<Fr>(begin, end);
    };
    for (size_t i = 0; i < input.com_s.size(); ++i) {
      auto ds = get_s(i) * d[i];
      y = y + ds;
    }

    // com_y_r, com_y
    Fr com_y_r = InnerProduct(input.com_s_r, d);
    G1 com_y =
        PcComputeCommitmentG(input.s_g_offset, y, com_y_r);  // multiexp(sn+1)
#ifdef _DEBUG
    G1 check_com_y = MultiExpBdlo12(input.com_s, d);  // multiexp(n/sn)
    assert(check_com_y == com_y);
#endif

    // prove <x, f> == <y, e>
    auto z = InnerProduct(input.x, f);
    assert(z == InnerProduct(y, e));

    typename EqualIp<HyraxA>::ProveInput eip_input(
        input.x, f, input.com_x, input.com_x_r, input.x_g_offset, y, e, com_y,
        com_y_r, input.s_g_offset, z);

    EqualIp<HyraxA>::Prove(proof, seed, eip_input);
  }

  struct VerifyInput {
    VerifyInput(int64_t sn, G1 const& com_x, int64_t x_g_offset,
                std::vector<G1> const& com_s, int64_t s_g_offset)
        : sn(sn),
          com_x(com_x),
          x_g_offset(x_g_offset),
          com_s(com_s),
          s_g_offset(s_g_offset) {}
    int64_t const sn;
    G1 const& com_x;
    int64_t const x_g_offset;
    std::vector<G1> const& com_s;
    int64_t const s_g_offset;
  };

  static bool Verify(h256_t seed, VerifyInput const& input,
                     Proof const& proof) {
    UpdateSeed(seed, input.com_x, input.sn, input.com_s);
    int64_t n = input.com_s.size() * input.sn;

    // d, e
    std::vector<Fr> d(input.com_s.size());
    ComputeFst(seed, "divide:d", d);
    std::vector<Fr> e(input.sn);
    ComputeFst(seed, "divide:e", e);

    // f
    std::vector<Fr> f;
    f.reserve(n);
    for (size_t i = 0; i < input.com_s.size(); ++i) {
      auto de = e * d[i];
      f.insert(f.end(), de.begin(), de.end());
    }

    // com_y
    G1 com_y = MultiExpBdlo12(input.com_s, d);  // multiexp(com_s.size())

    typename EqualIp<HyraxA>::VerifyInput eip_input(
        f, input.com_x, input.x_g_offset, e, com_y, input.s_g_offset);
    return EqualIp<HyraxA>::Verify(seed, proof, eip_input);
  }

  static bool Test();
};

template <typename HyraxA>
bool Divide<HyraxA>::Test() {
  auto seed = misc::RandH256();

  int64_t x_g_offset = 50;
  int64_t s_g_offset = 770;
  std::vector<Fr> x(15);
  FrRand(x);
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(x_g_offset, x, com_x_r);

  int64_t sn = 3;
  int64_t m = x.size() / sn;
  std::vector<G1> com_s(m);
  std::vector<Fr> com_s_r(m);
  for (auto i = 0; i < m; ++i) {
    auto begin = x.begin() + i * sn;
    auto end = begin + sn;
    std::vector<Fr> s(begin, end);
    com_s_r[i] = FrRand();
    com_s[i] = PcComputeCommitmentG(s_g_offset, s, com_s_r[i]);
  }

  ProveInput prove_input(x, com_x, com_x_r, x_g_offset, sn, com_s, com_s_r,
                         s_g_offset);
  Proof proof;
  Prove(proof, seed, prove_input);

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

  VerifyInput verify_input(sn, com_x, x_g_offset, com_s, s_g_offset);
  bool success = Verify(seed, verify_input, proof);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink