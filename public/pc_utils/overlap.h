#pragma once

#include "./details.h"
#include "./equal_ip.h"
#include "./types.h"

// for vector x and y, open com(gx, x) and com(gy, y), prove overlap

namespace pc_utils {

template <typename HyraxA>
struct Overlap {
  using Proof = typename EqualIp<HyraxA>::Proof;

  struct OverlapPosition {
    int64_t x;
    int64_t y;
  };

  static void UpdateSeed(h256_t& seed, G1 const& c1, G1 const& c2, int64_t xn,
                         int64_t yn,
                         std::vector<OverlapPosition> const& overlap) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, c1);
    HashUpdate(hash, c2);
    HashUpdate(hash, xn);
    HashUpdate(hash, yn);
    for (auto& i : overlap) {
      HashUpdate(hash, i.x);
      HashUpdate(hash, i.y);
    }
    hash.Final(seed.data());
  }

  struct ProverInput {
    ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& y,
                std::vector<OverlapPosition> const& overlap, G1 const& com_x,
                G1 const& com_y, Fr const& com_x_r, Fr const& com_y_r,
                int64_t x_g_offset, int64_t y_g_offset)
        : x(x),
          y(y),
          overlap(overlap),
          com_x(com_x),
          com_y(com_y),
          com_x_r(com_x_r),
          com_y_r(com_y_r),
          x_g_offset(x_g_offset),
          y_g_offset(y_g_offset) {}
    std::vector<Fr> const& x;
    std::vector<Fr> const& y;
    std::vector<OverlapPosition> const& overlap;
    G1 const& com_x;
    G1 const& com_y;
    Fr const& com_x_r;
    Fr const& com_y_r;
    int64_t const x_g_offset;
    int64_t const y_g_offset;
    int64_t xn() const { return (int64_t)x.size(); }
    int64_t yn() const { return (int64_t)y.size(); }
  };

  static void Prove(Proof& proof, h256_t seed, ProverInput const& input) {
    int64_t xn = input.xn();
    int64_t yn = input.yn();

    UpdateSeed(seed, input.com_x, input.com_y, xn, yn, input.overlap);
    std::vector<Fr> c(input.overlap.size());
    ComputeFst(seed, "consistency::overlap::c", c);

    std::vector<Fr> a(xn, FrZero());
    std::vector<Fr> b(yn, FrZero());
    for (size_t i = 0; i < input.overlap.size(); ++i) {
      auto const& o = input.overlap[i];
      assert(input.x[o.x] == input.y[o.y]);
      a[o.x] += c[i];  // NOTE: here is +=
      b[o.y] += c[i];
    }

    Fr z = InnerProduct(input.x, a);
    assert(z == InnerProduct(input.y, b));

    EqualIp<HyraxA>::ProverInput eip_input(
        input.x, a, input.com_x, input.com_x_r, input.x_g_offset, input.y, b,
        input.com_y, input.com_y_r, input.y_g_offset, z);
    EqualIp<HyraxA>::Prove(proof, seed, eip_input);
  }

  struct VerifierInput {
    VerifierInput(int64_t xn, G1 const& com_x, int64_t x_g_offset, int64_t yn,
                  G1 const& com_y, int64_t y_g_offset,
                  std::vector<OverlapPosition> const& overlap)
        : xn(xn),
          com_x(com_x),
          x_g_offset(x_g_offset),
          yn(yn),
          com_y(com_y),
          y_g_offset(y_g_offset),
          overlap(overlap) {}
    int64_t const xn;
    G1 const& com_x;
    int64_t const x_g_offset;
    int64_t const yn;
    G1 const& com_y;
    int64_t const y_g_offset;
    std::vector<OverlapPosition> const& overlap;
  };

  static bool Verify(h256_t seed, VerifierInput const& input,
                     Proof const& proof) {
    UpdateSeed(seed, input.com_x, input.com_y, input.xn, input.yn,
               input.overlap);
    std::vector<Fr> c(input.overlap.size());
    ComputeFst(seed, "consistency::overlap::c", c);

    std::vector<Fr> a(input.xn, FrZero());
    std::vector<Fr> b(input.yn, FrZero());
    for (size_t i = 0; i < input.overlap.size(); ++i) {
      auto const& o = input.overlap[i];
      a[o.x] += c[i];  // NOTE: here is +=
      b[o.y] += c[i];
    }

    EqualIp<HyraxA>::VerifierInput eip_input(a, input.com_x, input.x_g_offset,
                                             b, input.com_y, input.y_g_offset);
    return EqualIp<HyraxA>::Verify(seed, proof, eip_input);
  }

  static bool Test();
};

template <typename HyraxA>
bool Overlap<HyraxA>::Test() {
  auto seed = misc::RandH256();

  int64_t x_g_offset = 20;
  int64_t y_g_offset = 40;

  std::vector<Fr> x(10);
  FrRand(x);
  std::vector<Fr> y(7);
  FrRand(y);

  y[0] = x[1];
  y[1] = x[3];
  y[2] = x[9];
  y[6] = x[9];
  y[4] = x[7];

  std::vector<OverlapPosition> overlap(5);
  overlap[0].x = 1;
  overlap[0].y = 0;
  overlap[1].x = 3;
  overlap[1].y = 1;
  overlap[2].x = 9;
  overlap[2].y = 2;
  overlap[3].x = 9;
  overlap[3].y = 6;
  overlap[4].x = 7;
  overlap[4].y = 4;

  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(x_g_offset, x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitmentG(y_g_offset, y, com_y_r);

  ProverInput prover_input(x, y, overlap, com_x, com_y, com_x_r, com_y_r,
                           x_g_offset, y_g_offset);

  Proof proof;
  Prove(proof, seed, prover_input);

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

  auto xn = (int64_t)x.size();
  auto yn = (int64_t)y.size();
  VerifierInput verifier_input(xn, com_x, x_g_offset, yn, com_y, y_g_offset,
                               overlap);
  bool success = Verify(seed, verifier_input, proof);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace pc_utils