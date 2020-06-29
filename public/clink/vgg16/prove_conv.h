#pragma once

#include "./context.h"
#include "./image_com.h"
#include "circuit/vgg16/vgg16.h"
#include "clink/equal_ip.h"

namespace clink::vgg16 {

inline void OneConvUpdateSeed(h256_t& seed, std::array<G1, 9> const& c) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  for (auto const& i : c) {
    HashUpdate(hash, i);
  }
  hash.Final(seed.data());
}

inline void OneConvUpdateSeed(h256_t& seed, G1 const& c) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, c);
  hash.Final(seed.data());
}

struct OneConvInputPub {
  std::array<G1, 9> cb;
  std::array<G1, 9> a3_cy;
  std::array<HyraxA::Proof, 9> rb_proofs;
  HyraxA::Proof xq_proof;
};

struct OneConvInputSec {
  std::vector<std::array<Fr const*, 9>> b;
  std::array<Fr, 9> rb;
};

inline void OneConvComputeFst(h256_t const& seed, std::string const& prefix,
                              size_t layer, size_t col, std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(layer);
  salt += " ";
  salt += std::to_string(col);
  ComputeFst(seed, salt, r);
}

inline void OneConvComputeFst(h256_t const& seed, std::string const& prefix,
                              std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(r.size());
  ComputeFst(seed, salt, r);
}

inline void OneConvComputeInputR(h256_t const& seed, size_t layer, size_t CDD,
                                 std::array<std::vector<Fr>, 9>& r) {
  auto parallel_f = [&seed, &r, layer, CDD](int64_t j) {
    r[j].resize(CDD);
    OneConvComputeFst(seed, "conv adapt input ", layer, j, r[j]);
  };
  parallel::For(9, parallel_f);
}

inline std::vector<Fr> OneConvComputeOutputR(h256_t const& seed, size_t layer,
                                             size_t KDD) {
  std::vector<Fr> r(KDD);
  std::string salt = "conv output ";
  salt += std::to_string(layer);
  OneConvComputeFst(seed, salt, r);
  return r;
}

inline std::vector<Fr> OneConvOutputR2S(size_t K, size_t C, size_t D,
                                        std::vector<Fr> const& r) {
  assert(r.size() == K * D * D);
  std::vector<Fr> s(K * C * D * D);
  for (size_t i = 0; i < K; ++i) {
    for (size_t j = 0; j < C; ++j) {
      for (size_t k = 0; k < D; ++k) {
        for (size_t l = 0; l < D; ++l) {
          s[i * C * D * D + j * D * D + k * D + l] = r[i * D * D + k * D + l];
        }
      }
    }
  }
  return s;
}

inline void OneConvInputProve(h256_t seed, ProveContext const& context,
                              size_t layer, OneConvInputPub& pub,
                              OneConvInputSec& sec) {
  Tick tick(__FN__);

#ifdef _DEBUG_CHECK
  if (kLayerTypeOrders[layer].first != kConv) {
    throw std::runtime_error("oops");
  }
#endif

  struct Ctx {
    std::array<Fr, 9> a3_y;
    std::array<Fr, 9> a3_ry;
    Fr a3_rz;
    G1 a3_cz;
    std::array<std::vector<Fr>, 9> r;
    std::vector<std::array<int64_t, 9>> B;
    std::vector<Fr> q;
    Fr ip_xq;  // <x,q>
  } ctx;

  size_t const order = kLayerTypeOrders[layer].second;
  auto const& input_image = *context.const_images()[layer];
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto const& x = input_image.pixels;
  Fr const& rx = context.image_com_sec().r[layer];
  G1 const& cx = context.image_com_pub().c[layer];
  FrRand(sec.rb.data(), sec.rb.size());

  // build sec.b and ctx.B
  ctx.B.resize(C * D * D);
  sec.b.resize(K * C * D * D);
  for (size_t i = 0; i < C * D * D; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      size_t m = j / 3;
      size_t n = j % 3;
      size_t r = i / (D * D);
      size_t p = i % (D * D);
      size_t q = p / D;
      size_t o = p % D;
      size_t ii = q + m;
      size_t jj = o + n;
      if (ii == 0 || jj == 0 || ii == (D + 1) || jj == (D + 1)) {
        sec.b[i][j] = &FrZero();
        ctx.B[i][j] = -1;
      } else {
        sec.b[i][j] = &x[r][ii - 1][jj - 1];
        ctx.B[i][j] = r * D * D + (ii - 1) * D + (jj - 1);
      }
    }
  }

  // commit every col of sec.b
  auto get_col_b_u = [&context, order](int64_t i) -> G1 const& {
    auto range = context.auxi().data_u_conv(order);
    if (i >= (range.second - range.first)) throw std::runtime_error("oops");
    return range.first[i];
  };

  auto parallel_f1 = [&sec, &pub, &get_col_b_u, C, D](int64_t j) {
    auto get_x = [&sec, j](int64_t i) -> Fr const& { return *sec.b[i][j]; };
    pub.cb[j] =
        pc::PcComputeCommitmentG(C * D * D, get_col_b_u, get_x, sec.rb[j]);
  };
  parallel::For(9, parallel_f1);

  // update seed by pub.cb
  OneConvUpdateSeed(seed, pub.cb);

  // compute challenge ctx.r base on fst
  OneConvComputeInputR(seed, layer, C * D * D, ctx.r);

  // build ctx.q base ctx.B and ctx.r
  ctx.q.resize(C * D * D);
  std::fill(ctx.q.begin(), ctx.q.end(), FrZero());
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < C * D * D; ++i) {
      auto const& Bij = ctx.B[i][j];
      if (Bij != -1) {
        ctx.q[Bij] += ctx.r[j][i];
      }
    }
  }

  // prove <ctx.r[j], vec_bj>
  FrRand(ctx.a3_ry.data(), ctx.a3_ry.size());
  auto parallel_f2 = [&ctx, &pub, &sec, &seed, C, D, layer,
                      &get_col_b_u](int64_t j) {
    std::vector<Fr> vec_bj(C * D * D);
    for (size_t i = 0; i < C * D * D; ++i) {
      vec_bj[i] = *sec.b[i][j];
    }
    ctx.a3_y[j] = InnerProduct(vec_bj, ctx.r[j]);
    pub.a3_cy[j] = pc::PcComputeCommitmentG(ctx.a3_y[j], ctx.a3_ry[j]);

    HyraxA::ProveInput a3_rb_prove_input(vec_bj, ctx.r[j], ctx.a3_y[j],
                                         get_col_b_u, pc::PcG(0));
    HyraxA::CommitmentSec a3_rb_com_sec(sec.rb[j], ctx.a3_ry[j]);
    HyraxA::CommitmentPub a3_rb_com_pub(pub.cb[j], pub.a3_cy[j]);
    HyraxA::Prove(pub.rb_proofs[j], seed, a3_rb_prove_input, a3_rb_com_pub,
                  a3_rb_com_sec);
  };
  parallel::For(9, parallel_f2);

  // prove <q, vec_x>
  std::vector<Fr> vec_x(C * D * D);
  for (size_t i = 0; i < C; ++i) {
    for (size_t j = 0; j < D; ++j) {
      for (size_t k = 0; k < D; ++k) {
        vec_x[i * D * D + j * D + k] = x[i][j][k];
      }
    }
  }

  ctx.ip_xq = std::accumulate(ctx.a3_y.begin(), ctx.a3_y.end(), FrZero());
  ctx.a3_rz = std::accumulate(ctx.a3_ry.begin(), ctx.a3_ry.end(), FrZero());
  G1 a3_cz = std::accumulate(pub.a3_cy.begin(), pub.a3_cy.end(), G1Zero());
#ifdef _DEBUG_CHECK
  if (InnerProduct(vec_x, ctx.q) != ctx.ip_xq) {
    throw std::runtime_error("oops");
  }
  if (a3_cz != pc::PcComputeCommitmentG(ctx.ip_xq, ctx.a3_rz)) {
    throw std::runtime_error("oops");
  }
#endif
  HyraxA::ProveInput a3_xq_prove_input(vec_x, ctx.q, ctx.ip_xq, pc::kGetRefG,
                                       pc::PcG(0));
  HyraxA::CommitmentSec a3_xq_com_sec(rx, ctx.a3_rz);
  HyraxA::CommitmentPub a3_xq_com_pub(cx, a3_cz);
  HyraxA::Prove(pub.xq_proof, seed, a3_xq_prove_input, a3_xq_com_pub,
                a3_xq_com_sec);

  // extend sec.b (to sec.b'), now sec.b is [KCDD][9]
  for (size_t i = C * D * D; i < K * C * D * D; ++i) {
    sec.b[i] = sec.b[i % (C * D * D)];
  }

#ifdef _DEBUG_CHECK
  auto parallel_f3 = [&pub, &sec, C, D, K](int64_t j) {
    auto get_x = [&sec, j](int64_t i) -> Fr const& { return *sec.b[i][j]; };
    auto c = pc::PcComputeCommitmentG(K * C * D * D, get_x, sec.rb[j]);
    if (pub.cb[j] != c) {
      throw std::runtime_error("oops");
    }
  };
  parallel::For(9, parallel_f3);
#endif
}

struct OneConvR1csPub {
  clink::ParallelR1cs<typename R1cs>::Proof r1cs_proof;
  std::vector<G1> com_w;
  size_t r1cs_ret_index;
};

struct OneConvR1csSec {
  std::vector<Fr> y;
  Fr ry;
};

// prove y=<x,para>
inline void OneConvR1csProve(h256_t seed, ProveContext const& context,
                             size_t layer, OneConvInputPub const& input_pub,
                             OneConvInputSec const& input_sec,
                             OneConvR1csPub& pub, OneConvR1csSec& sec) {
  Tick tick(__FN__);
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto order = kLayerTypeOrders[layer].second;

#ifdef _DEBUG_CHECK
  if (input_sec.b.size() != K * C * D * D) throw std::runtime_error("oops");
#endif

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::IpGadget gadget(pb, "vgg16 conv gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  pub.r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  auto s = r1cs_info->num_variables;
  std::vector<std::vector<Fr>> w(s);
  auto n = K * C * D * D;
  for (auto& i : w) i.resize(n);

  // convert para from K*C*3*3 to KCDD*9
  Para::ConvLayer const& para_conv = context.para().conv_layer(order);
  std::vector<std::array<Fr const*, 9>> p(K * C * D * D);
  for (size_t i = 0; i < K * C * D * D; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      auto ii = i / (D * D);
      p[i][j] = &para_conv.coefs[ii / C][ii % C][j / 3][j % 3];
    }
  }

  for (size_t j = 0; j < n; ++j) {
    gadget.Assign(input_sec.b[j], p[j]);
    assert(pb.is_satisfied());
    auto v = pb.full_variable_assignment();
    for (int64_t i = 0; i < s; ++i) {
      w[i][j] = v[i];
    }
    for (int64_t i = 0; i < 9; ++i) {
      assert(w[i][j] == *input_sec.b[j][i]);
    }
    for (int64_t i = 0; i < 9; ++i) {
      assert(w[i + 9][j] == *p[j][i]);
    }
  }

  pub.com_w.resize(s);
  std::vector<Fr> com_w_r(s);

  std::cout << "compute conv com(witness)\n";
  for (size_t i = 0; i < 9; ++i) {
    com_w_r[i] = input_sec.rb[i];
    pub.com_w[i] = input_pub.cb[i];
    com_w_r[i + 9] = context.para_com_sec().conv.coef_r[order][i];
    pub.com_w[i + 9] = context.para_com_pub().conv.coef[order][i];
  }

  auto parallel_f = [&com_w_r, &pub, &w, order](int64_t i) {
    std::vector<Fr> const& x = w[i + 18];
    Fr& r = com_w_r[i + 18];
    G1& c = pub.com_w[i + 18];
    r = FrRand();
    c = pc::PcComputeCommitmentG(x, r, true);
  };
  parallel::For<int64_t>(s - 18, parallel_f);

  // save output
  sec.y = w[pub.r1cs_ret_index];
  sec.ry = com_w_r[pub.r1cs_ret_index];

  typename R1cs::ProveInput r1cs_input(*r1cs_info, std::move(w), pub.com_w,
                                       com_w_r, pc::kGetRefG);

  R1cs::Prove(pub.r1cs_proof, seed, std::move(r1cs_input));
}

struct OneConvOutputPub {
  G1 cy;
  EqualIp<HyraxA>::Proof eip_proof;
};

inline void OneConvOutputProve(h256_t seed, ProveContext const& context,
                               size_t layer, OneConvR1csPub const& r1cs_pub,
                               OneConvR1csSec const& r1cs_sec,
                               OneConvOutputPub& pub) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  size_t const order = kLayerTypeOrders[layer].second;  
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;  

  G1 const& cx = r1cs_pub.com_w[r1cs_pub.r1cs_ret_index];
  Fr const& rx = r1cs_sec.ry;
  std::vector<Fr> const& x = r1cs_sec.y;

#ifdef _DEBUG_CHECK
  if (x.size() != K * C * D * D) throw std::runtime_error("oops");
#endif

  std::vector<Fr> y(K * D * D, FrZero());
  Fr ry;

  for (size_t i = 0; i < K; ++i) {
    for (size_t l = 0; l < C; ++l) {
      for (size_t j = 0; j < D; ++j) {
        for (size_t k = 0; k < D; ++k) {
          y[i * D * D + j * D + k] += x[i * C * D * D + l * D * D + j * D + k];
        }
      }
    }
  }

  Fr const& rz = context.image_com_sec().r[layer + 1];
  G1 const& cz = context.image_com_pub().c[layer + 1];
  Fr rb = context.para_com_sec().conv.bias_r[order] *
          fp::RationalConst<8, 24>().kFrN;
  G1 cb =
      context.para_com_pub().conv.bias[order] * fp::RationalConst<8, 24>().kFrN;

  ry = rz - rb;
  pub.cy = pc::PcComputeCommitmentG(y, ry);

  if (cz != pub.cy + cb) {
    throw std::runtime_error("oops");
  }

#ifdef _DEBUG_CHECK
  auto const& output_image = *context.const_images()[layer + 1];
  Para::ConvLayer const& para_conv = context.para().conv_layer(order);
  auto const& bias = para_conv.bias;  // vector<Fr>, size=K
  for (size_t i = 0; i < K; ++i) {
    for (size_t j = 0; j < D; ++j) {
      for (size_t k = 0; k < D; ++k) {
        auto offset = i * D * D + j * D + k;
        auto b = bias[i] * fp::RationalConst<8, 24>().kFrN;
        if (y[offset] + b != output_image.pixels[i][j][k]) {
          throw std::runtime_error("oops");
        }
      }
    }
  }
#endif

  OneConvUpdateSeed(seed, pub.cy);

  // r.size = KDD
  std::vector<Fr> r = OneConvComputeOutputR(seed, layer, K * D * D);

  // s.size = KCDD
  std::vector<Fr> s = OneConvOutputR2S(K, C, D, r);

  Fr ip = InnerProduct(y, r);

#ifdef _DEBUG_CHECK
  if (InnerProduct(x, s) != ip) {
    throw std::runtime_error("oops");
  }
#endif

  // prove <x,s>==<y,r>
  EqualIp<HyraxA>::ProveInput eip_input(x, s, cx, rx, pc::kGetRefG, y, r,
                                        pub.cy, ry, pc::kGetRefG, ip);
  EqualIp<HyraxA>::Prove(pub.eip_proof, seed, eip_input);
}

struct OneConvProof {
  OneConvInputPub input;
  OneConvR1csPub r1cs;
  OneConvOutputPub output;
};

inline void OneConvProve(h256_t seed, ProveContext const& context, size_t layer,
                         OneConvProof& proof) {
  Tick tick(__FN__);
  OneConvInputPub input_pub;
  OneConvInputSec input_sec;
  OneConvInputProve(seed, context, layer, input_pub, input_sec);

  OneConvR1csPub r1cs_pub;
  OneConvR1csSec r1cs_sec;
  OneConvR1csProve(seed, context, layer, input_pub, input_sec, r1cs_pub,
                   r1cs_sec);

  OneConvOutputPub output_pub;
  OneConvOutputProve(seed, context, layer, r1cs_pub, r1cs_sec, output_pub);

  proof.input = input_pub;
  proof.r1cs = r1cs_pub;
  proof.output = output_pub;
}
}  // namespace clink::vgg16