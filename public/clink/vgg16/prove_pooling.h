#pragma once
#
#include "./pooling_pub.h"

namespace clink::vgg16 {


inline void PoolingInputProvePreprocess(h256_t seed,
                                        ProveContext const& context,
                                        PoolingProof& proof,
                                        PoolingInputSec& input_sec,
                                        AdaptProveItemMan& item_man) {
  Tick tick(__FN__);
  PoolingInputPub& input_pub = proof.input_pub;
  std::array<std::vector<Fr>, 9> x;
  std::array<Fr, 9> com_x_r;

  for (size_t i = 0; i < 5; ++i) {
    auto layer = kPoolingLayers[i];
    input_pub.cx[i] = context.image_com_pub().c[layer];
    com_x_r[i] = context.image_com_sec().r[layer];
    x[i] = context.const_images()[layer]->data;
  }

  // a=x[5],b=x[6],c=x[7],d=x[8]
  PoolingBuildAbcd(x);

  auto parallel_f = [&com_x_r, &input_pub, &x](int64_t i) {
    com_x_r[i] = FrRand();
    input_pub.cx[i] = pc::PcComputeCommitmentG(x[i], com_x_r[i]);
  };
  parallel::For<int64_t>(5, 9, parallel_f);

  PoolingUpdateSeed(seed, input_pub.cx[5], input_pub.cx[6], input_pub.cx[7],
                    input_pub.cx[8]);

  std::array<std::vector<Fr>, 9> q;
  PoolingBuildInputQ(seed, q);

  // now we have 9 x and 9 q
  std::array<Fr, 9> xq;
  std::array<Fr, 9> com_xq_r;
  std::array<G1, 9> com_xq;

  SelectInputComIpR(com_xq_r);
  auto parallel_f2 = [&x, &q, &xq, &com_xq_r, &com_xq](int64_t i) {
    xq[i] = InnerProduct(x[i], q[i]);
    com_xq[i] = pc::PcComputeCommitmentG(xq[i], com_xq_r[i]);
  };
  parallel::For<int64_t>(9, parallel_f2);

#ifdef _DEBUG_CHECK
  Fr sum_xq_5 = std::accumulate(xq.begin(), xq.begin() + 5, FrZero());
  Fr sum_xq_4 = std::accumulate(xq.begin() + 5, xq.end(), FrZero());
  if (sum_xq_5 != sum_xq_4) throw std::runtime_error("oops");
  G1 sum_com_xq_5 =
      std::accumulate(com_xq.begin(), com_xq.begin() + 5, G1Zero());
  G1 sum_com_xq_4 = std::accumulate(com_xq.begin() + 5, com_xq.end(), G1Zero());
  if (sum_com_xq_5 != sum_com_xq_4) throw std::runtime_error("oops");
#endif

  input_sec.a = x[5];
  input_sec.b = x[6];
  input_sec.c = x[7];
  input_sec.d = x[8];
  input_sec.r_a = com_x_r[5];
  input_sec.r_b = com_x_r[6];
  input_sec.r_c = com_x_r[7];
  input_sec.r_d = com_x_r[8];

  AdaptProveItem adapt_item;
  adapt_item.Init(9, "pooling_input");
  for (size_t j = 0; j < 9; ++j) {
    adapt_item.x[j] = std::move(x[j]);
    adapt_item.a[j] = std::move(q[j]);
    adapt_item.cx[j] = input_pub.cx[j];
    adapt_item.rx[j] = com_x_r[j];
    if (j < 5) adapt_item.a[j] = -adapt_item.a[j];
  }
  item_man.emplace(std::move(adapt_item));
}

inline void PoolingR1csProvePreprocess(
    h256_t seed, ProveContext const& /*context*/,
    PoolingInputSec const& input_sec, PoolingProof& proof,
    std::shared_ptr<PoolingR1csSec> pr1cs_sec, ParallelVoidTaskMan& task_man) {
  Tick tick(__FN__);
  auto const& input_pub = proof.input_pub;
  auto& r1cs_pub = proof.r1cs_pub;
  auto& r1cs_proof = proof.r1cs_proof;
  auto& r1cs_sec = *pr1cs_sec;
  libsnark::protoboard<Fr> pb;
  circuit::vgg16::PoolingGadget<8, 24> gadget(pb, "vgg16 pooling gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  // see protoboard<FieldT>::val
  r1cs_pub.r1cs_ret_index = gadget.ret().index - 1;
  r1cs_sec.r1cs_info.reset(new R1csInfo(pb));
  auto s = r1cs_sec.r1cs_info->num_variables;
  r1cs_sec.w.resize(s);
  auto n = input_sec.a.size();
  for (auto& i : r1cs_sec.w) i.resize(n);
  std::cout << "PoolingR1csProvePreprocess: " << r1cs_sec.r1cs_info->to_string()
            << ", repeat times: " << n << "\n";

#ifdef _DEBUG_CHECK
  if (input_sec.a.size() != n || input_sec.b.size() != n ||
      input_sec.c.size() != n || input_sec.d.size() != n) {
    throw std::runtime_error("oops");
  }
  if (n != PoolingGetCircuitCount()) {
    throw std::runtime_error("oops");
  }
#endif

  for (size_t j = 0; j < n; ++j) {
    std::array<Fr const*, 4> data = {&input_sec.a[j], &input_sec.b[j],
                                     &input_sec.c[j], &input_sec.d[j]};
    gadget.Assign(data);
    assert(pb.is_satisfied());
    auto v = pb.full_variable_assignment();
    for (int64_t i = 0; i < s; ++i) {
      r1cs_sec.w[i][j] = v[i];
    }
  }

  r1cs_pub.com_w.resize(s);
  r1cs_sec.com_w_r.resize(s);

  std::cout << "compute pooling com(witness)\n";
  r1cs_sec.com_w_r[0] = input_sec.r_a;
  r1cs_pub.com_w[0] = input_pub.cx[5];
  r1cs_sec.com_w_r[1] = input_sec.r_b;
  r1cs_pub.com_w[1] = input_pub.cx[6];
  r1cs_sec.com_w_r[2] = input_sec.r_c;
  r1cs_pub.com_w[2] = input_pub.cx[7];
  r1cs_sec.com_w_r[3] = input_sec.r_d;
  r1cs_pub.com_w[3] = input_pub.cx[8];

  auto parallel_f = [&r1cs_sec, &r1cs_pub](int64_t i) {
    std::vector<Fr> const& x = r1cs_sec.w[i + 4];
    Fr& r = r1cs_sec.com_w_r[i + 4];
    G1& c = r1cs_pub.com_w[i + 4];
    r = FrRand();
    c = pc::PcComputeCommitmentG(x, r, true);
  };
  parallel::For<int64_t>(s - 4, parallel_f);

  // save output
  r1cs_sec.y = r1cs_sec.w[r1cs_pub.r1cs_ret_index];
  r1cs_sec.ry = r1cs_sec.com_w_r[r1cs_pub.r1cs_ret_index];

  parallel::VoidTask task = [seed, &r1cs_pub, pr1cs_sec, &r1cs_proof]() {
    typename R1cs::ProveInput r1cs_input(
        *pr1cs_sec->r1cs_info, std::move(pr1cs_sec->w), r1cs_pub.com_w,
        pr1cs_sec->com_w_r, pc::kGetRefG);

    R1cs::Prove(r1cs_proof, seed, std::move(r1cs_input));
  };
  task_man.emplace(std::move(task));
}

inline void PoolingOutputProvePreprocess(h256_t seed,
                                         ProveContext const& context,                                         
                                         PoolingR1csSec const& r1cs_sec,                                         
                                         PoolingProof& proof,
                                         AdaptProveItemMan& item_man) {
  Tick tick(__FN__);
  auto const& r1cs_pub = proof.r1cs_pub;
  auto& output_pub = proof.output_pub;
  std::array<std::vector<Fr>, 6> x;
  std::array<Fr, 6> com_x_r;
  x[0] = r1cs_sec.y;
  output_pub.cx[0] = r1cs_pub.com_w[r1cs_pub.r1cs_ret_index];
  com_x_r[0] = r1cs_sec.ry;

  for (size_t i = 0; i < kPoolingLayers.size(); ++i) {
    auto layer = kPoolingLayers[i] + 1;  // next
    auto image = context.const_images()[layer];
    x[i + 1] = image->data;
    output_pub.cx[i + 1] = context.image_com_pub().c[layer];
    com_x_r[i + 1] = context.image_com_sec().r[layer];
  }

#ifdef _DEBUG_CHECK
  std::vector<Fr> check_x;
  for (size_t i = 1; i < x.size(); ++i) {
    check_x.insert(check_x.end(), x[i].begin(), x[i].end());
  }
  if (x[0] != check_x) throw std::runtime_error("oops");
#endif

  PoolingUpdateSeed(seed, r1cs_pub);

  std::array<std::vector<Fr>, 6> q;
  PoolingBuildOutputQ(seed, q);

#ifdef _DEBUG_CHECK
  // now we have 6 x and 6 q
  std::array<Fr, 6> xq;
  std::array<Fr, 6> com_xq_r;
  std::array<G1, 6> com_xq;
  SelectOutputComIpR(com_xq_r);
  auto parallel_f = [&x, &q, &xq, &com_xq_r, &com_xq](int64_t i) {
    xq[i] = InnerProduct(x[i], q[i]);
    com_xq[i] = pc::PcComputeCommitmentG(xq[i], com_xq_r[i]);
  };
  parallel::For<int64_t>(xq.size(), parallel_f);

  Fr sum_xq = std::accumulate(xq.begin() + 1, xq.end(), FrZero());
  if (sum_xq != xq[0]) throw std::runtime_error("oops");
  G1 sum_com_xq = std::accumulate(com_xq.begin() + 1, com_xq.end(), G1Zero());
  if (sum_com_xq != com_xq[0]) throw std::runtime_error("oops");
#endif

  // auto parallel_f2 = [&seed, &x, &q, &xq, &com_xq, &com_x_r, &com_xq_r,
  //                    &output_pub](int64_t i) {
  //  HyraxA::ProveInput input(x[i], q[i], xq[i], pc::kGetRefG, pc::PcG(0));
  //  output_pub.ip_com_pubs[i].tau = com_xq[i];
  //  HyraxA::CommitmentSec ip_com_sec(com_x_r[i], com_xq_r[i]);
  //  HyraxA::Prove(output_pub.ip_proofs[i], seed, input,
  //  output_pub.ip_com_pubs[i],
  //                ip_com_sec);
  //};
  // parallel::For(6, parallel_f2);

  AdaptProveItem adapt_item;
  adapt_item.Init(6, "pooling_output");
  for (size_t i = 0; i < 6; ++i) {
    adapt_item.x[i] = std::move(x[i]);
    adapt_item.a[i] = std::move(q[i]);
    adapt_item.cx[i] = output_pub.cx[i];
    adapt_item.rx[i] = com_x_r[i];
    if (i == 0) {
      adapt_item.a[i] = -adapt_item.a[i];
    }
  }
  item_man.emplace(std::move(adapt_item));
}

inline void PoolingProvePreprocess(h256_t seed, ProveContext const& context,
                                   PoolingProof& proof,
                                   AdaptProveItemMan& item_man,
                                   ParallelVoidTaskMan& task_man) {
  Tick tick(__FN__);
  PoolingInputSec input_sec;
  PoolingInputProvePreprocess(seed, context, proof, input_sec, item_man);

  std::shared_ptr<PoolingR1csSec> r1cs_sec(new PoolingR1csSec);
  PoolingR1csProvePreprocess(seed, context, input_sec, proof, r1cs_sec,
                             task_man);

  PoolingOutputProvePreprocess(seed, context, *r1cs_sec, proof, item_man);
}
}  // namespace clink::vgg16