#pragma once

#include "./details.h"
#include "./equal_ip.h"

// for vector x and y, open com(gx, x) and com(gy, y), prove overlap

namespace clink {

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

  struct ProveInput {
    ProveInput(std::vector<Fr> const& x, std::vector<Fr> const& y,
               std::vector<OverlapPosition> const& overlap, G1 const& com_x,
               G1 const& com_y, Fr const& com_x_r, Fr const& com_y_r,
               pc::GetRefG const& get_gx, pc::GetRefG const& get_gy)
        : x(x),
          y(y),
          overlap(overlap),
          com_x(com_x),
          com_y(com_y),
          com_x_r(com_x_r),
          com_y_r(com_y_r),
          get_gx(get_gx),
          get_gy(get_gy) {}
    std::vector<Fr> const& x;
    std::vector<Fr> const& y;
    std::vector<OverlapPosition> const& overlap;
    G1 const& com_x;
    G1 const& com_y;
    Fr const& com_x_r;
    Fr const& com_y_r;
    pc::GetRefG const& get_gx;
    pc::GetRefG const& get_gy;
    int64_t xn() const { return (int64_t)x.size(); }
    int64_t yn() const { return (int64_t)y.size(); }
  };

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input) {
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

    typename EqualIp<HyraxA>::ProveInput eip_input(
        input.x, a, input.com_x, input.com_x_r, input.get_gx, input.y, b,
        input.com_y, input.com_y_r, input.get_gy, z);
    EqualIp<HyraxA>::Prove(proof, seed, eip_input);
  }

  struct VerifyInput {
    VerifyInput(int64_t xn, G1 const& com_x, pc::GetRefG const& get_gx, int64_t yn,
                G1 const& com_y, pc::GetRefG const& get_gy,
                std::vector<OverlapPosition> const& overlap)
        : xn(xn),
          com_x(com_x),
          get_gx(get_gx),
          yn(yn),
          com_y(com_y),
          get_gy(get_gy),
          overlap(overlap) {}
    int64_t const xn;
    G1 const& com_x;
    pc::GetRefG const& get_gx;
    int64_t const yn;
    G1 const& com_y;
    pc::GetRefG const& get_gy;
    std::vector<OverlapPosition> const& overlap;
  };

  static bool Verify(h256_t seed, VerifyInput const& input,
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

    typename EqualIp<HyraxA>::VerifyInput eip_input(
        a, input.com_x, input.get_gx, b, input.com_y, input.get_gy);
    return EqualIp<HyraxA>::Verify(seed, proof, eip_input);
  }

  static bool Test();
};

template <typename HyraxA>
bool Overlap<HyraxA>::Test() {
  auto seed = misc::RandH256();

  int64_t x_g_offset = 20;
  int64_t y_g_offset = 40;
  pc::GetRefG get_gx = [x_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[x_g_offset + i];
  };
  pc::GetRefG get_gy = [y_g_offset](int64_t i) -> G1 const& {
    return pc::PcG()[y_g_offset + i];
  };

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
  G1 com_x = pc::PcComputeCommitmentG(get_gx, x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = pc::PcComputeCommitmentG(get_gy, y, com_y_r);

  ProveInput prove_input(x, y, overlap, com_x, com_y, com_x_r, com_y_r,
                         get_gx, get_gy);

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

  auto xn = (int64_t)x.size();
  auto yn = (int64_t)y.size();
  VerifyInput verify_input(xn, com_x, get_gx, yn, com_y, get_gy,
                           overlap);
  bool success = Verify(seed, verify_input, proof);
  std::cout << __FILE__ << " " << __FN__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}
}  // namespace clink