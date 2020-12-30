#pragma once

#include "./conv_pub.h"

namespace clink::vgg16 {

inline void OneConvInputProvePreprocess(h256_t seed,
                                        ProveContext const& context,
                                        size_t layer, OneConvProof& proof,
                                        OneConvInputSec& input_sec,
                                        SafeVec<AdaptProveItem>& adapt_man) {
  Tick tick(__FN__, std::to_string(layer));

  CHECK(kLayerTypeOrders[layer].first == kConv, "");

  struct Ctx {
    std::array<std::vector<Fr>, 9> r;
    std::vector<std::array<int64_t, 9>> B;
    std::vector<Fr> q;
  } ctx;

  OneConvInputPub& input_pub = proof.input_pub;
  size_t const order = kLayerTypeOrders[layer].second;
  auto const& input_image = *context.const_images()[layer];
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto DD = D * D;
  auto CDD = C * DD;
  auto KCDD = K * CDD;
  auto const& x = input_image.pixels;
  auto const& vec_x = input_image.data;
  Fr const& rx = context.image_com_sec().r[layer];
  G1 const& cx = context.image_com_pub().c[layer];
  FrRand(input_sec.rb.data(), input_sec.rb.size());

  // build input_sec.b and ctx.B
  ctx.B.resize(KCDD);
  input_sec.b.resize(KCDD);
  for (size_t i = 0; i < CDD; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      size_t m = j / 3;
      size_t n = j % 3;
      size_t r = i / DD;
      size_t p = i % DD;
      size_t q = p / D;
      size_t o = p % D;
      size_t ii = q + m;
      size_t jj = o + n;
      if (ii == 0 || jj == 0 || ii == (D + 1) || jj == (D + 1)) {
        input_sec.b[i][j] = &FrZero();
        ctx.B[i][j] = -1;
      } else {
        input_sec.b[i][j] = &x[r][ii - 1][jj - 1];
        ctx.B[i][j] = r * D * D + (ii - 1) * D + (jj - 1);
      }
    }
  }

  // commit every col of input_sec.b
  auto get_col_b_u = [&context, order](int64_t i) -> G1 const& {
    auto range = context.auxi().data_u_conv(order);
    CHECK(i < (range.second - range.first), "");
    return range.first[i];
  };

  auto parallel_f1 = [&input_sec, &input_pub, &get_col_b_u, C, D](int64_t j) {
    auto get_x = [&input_sec, j](int64_t i) -> Fr const& {
      return *input_sec.b[i][j];
    };
    input_pub.cb[j] =
        pc::ComputeCom(C * D * D, get_col_b_u, get_x, input_sec.rb[j]);
  };
  parallel::For(9, parallel_f1);

  // extend input_sec.b (to input_sec.b'), now input_sec.b is [KCDD][9]
  for (size_t i = CDD; i < KCDD; ++i) {
    input_sec.b[i] = input_sec.b[i % (CDD)];
    ctx.B[i] = ctx.B[i % (CDD)];
  }

  if (DEBUG_CHECK) {
    auto parallel_f3 = [&input_pub, &input_sec, KCDD](int64_t j) {
      auto get_x = [&input_sec, j](int64_t i) -> Fr const& {
        return *input_sec.b[i][j];
      };
      CHECK(input_pub.cb[j] == pc::ComputeCom(KCDD, get_x, input_sec.rb[j]),
            "");
    };
    parallel::For(9, parallel_f3);
  }

  // update seed by input_pub.cb
  OneConvUpdateSeed(seed, input_pub.cb);

  // compute challenge ctx.r base on fst
  OneConvComputeInputR(seed, layer, KCDD, ctx.r);

  // build ctx.q base ctx.B and ctx.r
  ctx.q.resize(CDD, FrZero());
  for (size_t j = 0; j < 9; ++j) {
    for (size_t i = 0; i < KCDD; ++i) {
      auto const& Bij = ctx.B[i][j];
      if (Bij != -1) {
        ctx.q[Bij] += ctx.r[j][i];
      }
    }
  }

  AdaptProveItem adapt_item;
  adapt_item.Init(10, "conv_input_" + std::to_string(layer), FrZero());
  for (size_t j = 0; j < 9; ++j) {
    adapt_item.x[j].resize(KCDD);
    for (size_t i = 0; i < KCDD; ++i) {
      adapt_item.x[j][i] = *input_sec.b[i][j];
    }
    adapt_item.a[j] = ctx.r[j];
    adapt_item.cx[j] = input_pub.cb[j];
    adapt_item.rx[j] = input_sec.rb[j];
  }
  adapt_item.x.back() = vec_x;
  adapt_item.a.back() = -ctx.q;
  adapt_item.cx.back() = cx;
  adapt_item.rx.back() = rx;

  adapt_man.emplace(std::move(adapt_item));
}

// prove y=<x,para>
inline void OneConvR1csProvePreprocess(
    h256_t seed, ProveContext const& context, size_t layer,
    OneConvInputSec const& input_sec, OneConvProof& proof,
    std::shared_ptr<OneConvR1csSec> pr1cs_sec,
    SafeVec<R1csProveItem>& r1cs_man) {
  Tick tick(__FN__, std::to_string(layer));
  (void)seed;
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;
  auto order = kLayerTypeOrders[layer].second;

  CHECK(input_sec.b.size() == K * C * D * D, "");

  auto const& input_pub = proof.input_pub;
  auto& r1cs_pub = proof.r1cs_pub;
  auto& r1cs_sec = *pr1cs_sec;
  libsnark::protoboard<Fr> pb;
  circuit::vgg16::IpGadget gadget(pb, "vgg16 conv gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  // see protoboard<FieldT>::val
  r1cs_pub.r1cs_ret_index = gadget.ret().index - 1;
  r1cs_sec.r1cs_info.reset(new R1csInfo(pb));
  auto s = r1cs_sec.r1cs_info->num_variables;
  r1cs_sec.w.resize(s);
  auto n = K * C * D * D;
  for (auto& i : r1cs_sec.w) i.resize(n);
  std::cout << Tick::GetIndentString() << r1cs_sec.r1cs_info->to_string()
            << ", repeat times: " << n << "\n";

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
      r1cs_sec.w[i][j] = v[i];
    }
    for (int64_t i = 0; i < 9; ++i) {
      assert(r1cs_sec.w[i][j] == *input_sec.b[j][i]);
    }
    for (int64_t i = 0; i < 9; ++i) {
      assert(r1cs_sec.w[i + 9][j] == *p[j][i]);
    }
  }

  r1cs_pub.com_w.resize(s);
  r1cs_sec.com_w_r.resize(s);

  // std::cout << "compute conv com(witness)\n";
  for (size_t i = 0; i < 9; ++i) {
    r1cs_sec.com_w_r[i] = input_sec.rb[i];
    r1cs_pub.com_w[i] = input_pub.cb[i];
    r1cs_sec.com_w_r[i + 9] = context.para_com_sec().conv.coef_r[order][i];
    r1cs_pub.com_w[i + 9] = context.para_com_pub().conv.coef[order][i];
  }

  auto parallel_f = [&r1cs_sec, &r1cs_pub](int64_t i) {
    std::vector<Fr> const& x = r1cs_sec.w[i + 18];
    Fr& r = r1cs_sec.com_w_r[i + 18];
    G1& c = r1cs_pub.com_w[i + 18];
    r = FrRand();
    c = pc::ComputeCom(x, r, true);
  };
  parallel::For<int64_t>(s - 18, parallel_f);

  // save output
  r1cs_sec.y = r1cs_sec.w[r1cs_pub.r1cs_ret_index];
  r1cs_sec.ry = r1cs_sec.com_w_r[r1cs_pub.r1cs_ret_index];

  R1csProveItem item;
  item.r1cs_sec = pr1cs_sec;  // save ref
  item.r1cs_input.reset(new R1cs::ProveInput(
      *pr1cs_sec->r1cs_info, ConvR1csTag(layer), std::move(pr1cs_sec->w),
      r1cs_pub.com_w, pr1cs_sec->com_w_r, pc::kGetRefG1));
  r1cs_man.emplace(std::move(item));
}

inline void OneConvOutputProvePreprocess(h256_t seed,
                                         ProveContext const& context,
                                         size_t layer,
                                         OneConvR1csSec const& r1cs_sec,
                                         OneConvProof& proof,
                                         SafeVec<AdaptProveItem>& item_man) {
  Tick tick(__FN__, std::to_string(layer));
  namespace fp = circuit::fp;
  size_t const order = kLayerTypeOrders[layer].second;
  auto K = kImageInfos[layer + 1].C;
  auto C = kImageInfos[layer].C;
  auto D = kImageInfos[layer].D;

  auto const& r1cs_pub = proof.r1cs_pub;
  auto& output_pub = proof.output_pub;
  G1 const& cx = r1cs_pub.com_w[r1cs_pub.r1cs_ret_index];
  Fr const& rx = r1cs_sec.ry;
  std::vector<Fr> const& x = r1cs_sec.y;

  CHECK(x.size() == K * C * D * D, "");

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
  output_pub.cy = pc::ComputeCom(y, ry);

  CHECK(cz == output_pub.cy + cb, "");

  if (DEBUG_CHECK) {
    auto const& output_image = *context.const_images()[layer + 1];
    Para::ConvLayer const& para_conv = context.para().conv_layer(order);
    auto const& bias = para_conv.bias;  // vector<Fr>, size=K
    for (size_t i = 0; i < K; ++i) {
      for (size_t j = 0; j < D; ++j) {
        for (size_t k = 0; k < D; ++k) {
          auto offset = i * D * D + j * D + k;
          auto b = bias[i] * fp::RationalConst<8, 24>().kFrN;
          CHECK(y[offset] + b == output_image.pixels[i][j][k], "");
        }
      }
    }
  }

  OneConvUpdateSeed(seed, output_pub.cy);

  // r.size = KDD
  std::vector<Fr> r = OneConvComputeOutputR(seed, layer, K * D * D);

  // s.size = KCDD
  std::vector<Fr> s = OneConvOutputR2S(K, C, D, r);

  DCHECK(InnerProduct(y, r) == InnerProduct(x, s), "");

  AdaptProveItem adapt_item;
  adapt_item.Init(2, "conv_output_" + std::to_string(layer), FrZero());
  adapt_item.x[0] = x;  // ref to r1cs_sec.y, so copy it
  adapt_item.a[0] = std::move(s);
  adapt_item.cx[0] = cx;
  adapt_item.rx[0] = rx;
  adapt_item.x[1] = std::move(y);
  adapt_item.a[1] = std::move(r);
  adapt_item.a[1] = -adapt_item.a[1];
  adapt_item.cx[1] = output_pub.cy;
  adapt_item.rx[1] = ry;

  item_man.emplace(std::move(adapt_item));
}

inline void OneConvProvePreprocess(h256_t seed, ProveContext const& context,
                                   size_t layer, OneConvProof& proof,
                                   SafeVec<AdaptProveItem>& adapt_man,
                                   SafeVec<R1csProveItem>& r1cs_man) {
  Tick tick(__FN__, std::to_string(layer));
  OneConvInputSec input_sec;
  OneConvInputProvePreprocess(seed, context, layer, proof, input_sec,
                              adapt_man);

  std::shared_ptr<OneConvR1csSec> r1cs_sec(new OneConvR1csSec);
  OneConvR1csProvePreprocess(seed, context, layer, input_sec, proof, r1cs_sec,
                             r1cs_man);

  OneConvOutputProvePreprocess(seed, context, layer, *r1cs_sec, proof,
                               adapt_man);
}

}  // namespace clink::vgg16