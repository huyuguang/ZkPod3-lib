#pragma once

#include "./details.h"
#include "./equal_ip.h"

// for matrix n*s, open com(gr,row) and com(gc,col), prove consistency of them.
// That is, prove they were generated by same matrix.

namespace clink {

// HyraxA: hyrax::A2 or hyrax::A3
template <typename HyraxA>
struct Matrix {
  using Proof = typename EqualIp<HyraxA>::Proof;
  typedef std::function<Fr const&(int64_t i, int64_t j)> GetMatrix;

  struct ProveInput {
    ProveInput(GetMatrix get_m, std::vector<G1> const& com_rows,
               std::vector<Fr> const& com_row_rs, int64_t r_g_offset,
               std::vector<G1> const& com_cols,
               std::vector<Fr> const& com_col_rs, int64_t c_g_offset)
        : get_m(get_m),
          com_rows(com_rows),
          com_row_rs(com_row_rs),
          r_g_offset(r_g_offset),
          com_cols(com_cols),
          com_col_rs(com_col_rs),
          c_g_offset(c_g_offset) {
      assert(n() == (int64_t)com_row_rs.size());
      assert(s() == (int64_t)com_col_rs.size());
    }
    int64_t n() const { return (int64_t)com_rows.size(); }
    int64_t s() const { return (int64_t)com_cols.size(); }
    GetMatrix get_m;
    std::vector<G1> const& com_rows;
    std::vector<Fr> const& com_row_rs;
    int64_t const r_g_offset;
    std::vector<G1> const& com_cols;
    std::vector<Fr> const& com_col_rs;
    int64_t const c_g_offset;
  };

  static void UpdateSeed(h256_t& seed, std::vector<G1> const& c1,
                         std::vector<G1> const& c2) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, c1);
    HashUpdate(hash, c2);
    hash.Final(seed.data());
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input) {
    int64_t n = input.n();
    int64_t s = input.s();

    UpdateSeed(seed, input.com_rows, input.com_cols);
    std::vector<Fr> d(n);
    ComputeFst(seed, "consistency::matrix::d", d);
    std::vector<Fr> e(s);
    ComputeFst(seed, "consistency::matrix::e", e);

    std::array<parallel::Task, 2> tasks;
    std::vector<Fr> v(s, FrZero());
    G1 com_v;
    Fr com_v_r;
    tasks[0] = [n, s, &d, &v, &input, &com_v, &com_v_r]() {
      for (int64_t j = 0; j < s; ++j) {
        for (int64_t i = 0; i < n; ++i) {
          v[j] += d[i] * input.get_m(i, j);
        }
      }
      // compute com(v)
      com_v = MultiExpBdlo12(input.com_rows, d);
      com_v_r = InnerProduct(input.com_row_rs, d);

#ifdef _DEBUG
      G1 check_com_v = PcComputeCommitmentG(input.r_g_offset, v, com_v_r);
      assert(com_v == check_com_v);
#endif
    };

    std::vector<Fr> w(n, FrZero());
    G1 com_w;
    Fr com_w_r;
    tasks[1] = [n, s, &e, &w, &input, &com_w, &com_w_r]() {
      for (int64_t i = 0; i < n; ++i) {
        for (int64_t j = 0; j < s; ++j) {
          w[i] += e[j] * input.get_m(i, j);
        }
      }
      // compute com(w)
      com_w = MultiExpBdlo12(input.com_cols, e);
      com_w_r = InnerProduct(input.com_col_rs, e);

#ifdef _DEBUG
      G1 check_com_w = PcComputeCommitmentG(input.c_g_offset, w, com_w_r);
      assert(com_w == check_com_w);
#endif
    };

    parallel::Invoke(tasks);

    Fr z = InnerProduct(v, e);
    assert(z == InnerProduct(w, d));

    typename EqualIp<HyraxA>::ProveInput eip_input(
        v, e, com_v, com_v_r, input.r_g_offset, w, d, com_w, com_w_r,
        input.c_g_offset, z);
    EqualIp<HyraxA>::Prove(proof, seed, eip_input);
  }

  struct VerifyInput {
    VerifyInput(std::vector<G1> const& com_rows, int64_t r_g_offset,
                std::vector<G1> const& com_cols, int64_t c_g_offset)
        : com_rows(com_rows),
          r_g_offset(r_g_offset),
          com_cols(com_cols),
          c_g_offset(c_g_offset) {}
    int64_t n() const { return (int64_t)com_rows.size(); }
    int64_t s() const { return (int64_t)com_cols.size(); }
    std::vector<G1> const& com_rows;
    int64_t const r_g_offset;
    std::vector<G1> const& com_cols;
    int64_t const c_g_offset;
  };

  static bool Verify(h256_t seed, VerifyInput const& input,
                     Proof const& proof) {
    int64_t n = input.n();
    int64_t s = input.s();

    UpdateSeed(seed, input.com_rows, input.com_cols);

    std::vector<Fr> d(n);
    ComputeFst(seed, "consistency::matrix::d", d);
    std::vector<Fr> e(s);
    ComputeFst(seed, "consistency::matrix::e", e);

    std::array<parallel::Task, 2> tasks;

    G1 com_v;
    tasks[0] = [&proof, &input, &d, &com_v]() {
      com_v = MultiExpBdlo12(input.com_rows, d);
    };

    G1 com_w;
    tasks[1] = [&proof, &input, &e, &com_w]() {
      com_w = MultiExpBdlo12(input.com_cols, e);
    };

    parallel::Invoke(tasks);

    typename EqualIp<HyraxA>::VerifyInput eip_input(e, com_v, input.r_g_offset,
                                                    d, com_w, input.c_g_offset);

    return EqualIp<HyraxA>::Verify(seed, proof, eip_input);
  }

  static bool Test(int64_t n, int64_t s);
};

template <typename HyraxA>
bool Matrix<HyraxA>::Test(int64_t n, int64_t s) {
  auto seed = misc::RandH256();
  int64_t r_g_offset = 10;
  int64_t c_g_offset = 30;
  std::vector<Fr> m(n * s);
  std::vector<Fr> com_row_rs(n);
  std::vector<Fr> com_col_rs(s);
  FrRand(m);
  FrRand(com_row_rs);
  FrRand(com_col_rs);
  std::vector<G1> com_rows(n);
  std::vector<G1> com_cols(s);
  for (int64_t i = 0; i < n; ++i) {
    com_rows[i] = PcComputeCommitmentG(
        r_g_offset, s,
        [s, i, &m](int64_t j) -> Fr const& { return m[i * s + j]; },
        com_row_rs[i]);
  }
  for (int64_t j = 0; j < s; ++j) {
    com_cols[j] = PcComputeCommitmentG(
        c_g_offset, n,
        [s, j, &m](int64_t i) -> Fr const& { return m[i * s + j]; },
        com_col_rs[j]);
  }

  auto get_m = [&m, s](int64_t i, int64_t j) -> Fr const& {
    return m[i * s + j];
  };

  ProveInput prove_input(std::move(get_m), com_rows, com_row_rs, r_g_offset,
                         com_cols, com_col_rs, c_g_offset);

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

  VerifyInput verify_input(com_rows, r_g_offset, com_cols, c_g_offset);
  bool success = Verify(seed, verify_input, proof);
  std::cout << __FILE__ << " " << __FN__ << ": " << success
            << "\n\n\n\n\n\n";
  return success;
}

}  // namespace clink