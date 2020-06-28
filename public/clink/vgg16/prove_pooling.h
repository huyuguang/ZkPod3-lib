#pragma once

#include "./image_com.h"
#include "./para_com.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

inline void PoolingUpdateSeed(h256_t& seed, G1 const& a, G1 const& b,
                              G1 const& c, G1 const& d) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  HashUpdate(hash, a);
  HashUpdate(hash, b);
  HashUpdate(hash, c);
  HashUpdate(hash, d);
  hash.Final(seed.data());
}

inline void PoolingComputeFst(h256_t const& seed, std::string const& prefix,
                              size_t order, std::vector<Fr>& r) {
  std::string salt = prefix;
  salt += std::to_string(order);
  ComputeFst(seed, salt, r);
}

std::array<size_t, 5> PoolingGetLayers() {
  std::array<size_t, 5> layers;
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first == kPooling) {
      layers[kLayerTypeOrders[i].second] = i;
    }
  }
  return layers;
}

inline size_t PoolingGetCircuitCount() {
  auto layers = PoolingGetLayers();
  size_t total_size = 0;
  for (size_t i = 0; i < layers.size(); ++i) {
    auto layer = layers[i];
    auto C = kImageInfos[layer].channel_count;
    auto D = kImageInfos[layer].dimension;
    total_size += C * D * D;
  }
  return total_size / 4;
}

struct PoolingInputImage {
  std::vector<Fr> const* x;
  size_t C;
  size_t D;
};

void PoolingImageToAbcd(std::array<PoolingInputImage, 5> const& images,
                    std::vector<Fr>& a, std::vector<Fr>& b, std::vector<Fr>& c,
                    std::vector<Fr>& d) {
  size_t total_size = 0;
  for (size_t i = 0; i < images.size(); ++i) {
    total_size += images[i].x->size();
  }

  a.resize(total_size / 4);
  b.resize(total_size / 4);
  c.resize(total_size / 4);
  d.resize(total_size / 4);

  size_t offset = 0;
  for (size_t l = 0; l < images.size(); ++l) {
    auto C = images[l].C;
    auto D = images[l].D;
    auto DD = D * D;
    auto D_2 = D / 2;
    auto DD_4 = DD / 4;
    auto const& x = *images[l].x;
    for (size_t i = 0; i < C; ++i) {
      for (size_t j = 0; j < D / 2; ++j) {
        for (size_t k = 0; k < D / 2; ++k) {
          auto idx = offset + i * DD_4 + j * D_2 + k;
          a[idx] = x[i * DD + 2 * j * D + 2 * k];
          b[idx] = x[i * DD + 2 * j * D + 2 * k + 1];
          c[idx] = x[i * DD + (2 * j + 1) * D + 2 * k];
          d[idx] = x[i * DD + (2 * j + 1) * D + 2 * k + 1];
        }
      }
    }
    offset += C * DD_4;
  }
}

// 5 pooling layer + a + b + c + d
struct PoolingInputPub {
  std::array<HyraxA::CommitmentPub, 9> ip_com_pubs;
  std::array<HyraxA::Proof, 9> ip_proofs;
};

struct PoolingInputSec {
  std::vector<Fr> a;
  std::vector<Fr> b;
  std::vector<Fr> c;
  std::vector<Fr> d;
  Fr r_a;
  Fr r_b;
  Fr r_c;
  Fr r_d;
};

inline void SelectInputComIpR(std::array<Fr, 9>& com_ip_r) {
  FrRand(com_ip_r.data(), 5);
  auto sum = std::accumulate(com_ip_r.begin(), com_ip_r.begin() + 5, FrZero());
  auto part = SplitFr(sum, 4);
  for (size_t i = 0; i < 4; ++i) {
    com_ip_r[i + 5] = part[i];
  }
}

inline void SelectOutputComIpR(std::array<Fr, 6>& com_ip_r) {
  com_ip_r[0] = FrRand();
  auto part = SplitFr(com_ip_r[0], 5);
  for (size_t i = 0; i < 5; ++i) {
    com_ip_r[i + 1] = part[i];
  }
}

inline void PoolingBuildInputQ(h256_t const& seed,
                               std::array<std::vector<Fr>, 9>& q) {
  std::array<PoolingInputImage, 5> images;
  auto layers = PoolingGetLayers();
  for (size_t i = 0; i < 5; ++i) {
    images[i].C = kImageInfos[layers[i]].channel_count;
    images[i].D = kImageInfos[layers[i]].dimension;
    q[i].resize(images[i].C * images[i].D * images[i].D);
    PoolingComputeFst(seed, "pooling input q", i, q[i]);
    images[i].x = &q[i];
  }

  PoolingImageToAbcd(images, q[5], q[6], q[7], q[8]);
}

inline void PoolingBuildAbcd(std::array<std::vector<Fr>, 9>& x) {
  auto layers = PoolingGetLayers();
  std::array<PoolingInputImage, 5> images;
  for (size_t i = 0; i < 5; ++i) {
    auto layer = layers[i];
    images[i].x = &x[i];
    images[i].C = kImageInfos[layer].channel_count;
    images[i].D = kImageInfos[layer].dimension;
  }
  PoolingImageToAbcd(images, x[5], x[6], x[7], x[8]);
}

inline void PoolingInputProve(h256_t seed, ProveContext const& context,
                              PoolingInputPub& pub, PoolingInputSec& sec) {
  Tick tick(__FN__);
  auto layers = PoolingGetLayers();
  assert(layers.size() == 5);

  std::array<std::vector<Fr>, 9> x;
  std::array<Fr, 9> com_x_r;

  for (size_t i = 0; i < 5; ++i) {
    auto layer = layers[i];
    pub.ip_com_pubs[i].xi = context.image_com_pub().c[layer];
    com_x_r[i] = context.image_com_sec().r[layer];
    x[i] = context.const_images()[layer]->copy_vec();
  }

  // a=x[5],b=x[6],c=x[7],d=x[8]
  PoolingBuildAbcd(x);

  auto parallel_f = [&com_x_r, &pub, &x](int64_t i) {
    com_x_r[i] = FrRand();
    pub.ip_com_pubs[i].xi = pc::PcComputeCommitmentG(x[i], com_x_r[i]);
  };
  parallel::For<int64_t>(5, 9, parallel_f);

  PoolingUpdateSeed(seed, pub.ip_com_pubs[5].xi, pub.ip_com_pubs[6].xi,
                    pub.ip_com_pubs[7].xi, pub.ip_com_pubs[8].xi);

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

  auto parallel_f3 = [&seed, &x, &q, &xq, &com_xq, &com_x_r, &com_xq_r,
                      &pub](int64_t i) {
    HyraxA::ProveInput input(x[i], q[i], xq[i], pc::kGetRefG, pc::PcG(0));
    pub.ip_com_pubs[i].tau = com_xq[i];
    HyraxA::CommitmentSec ip_com_sec(com_x_r[i], com_xq_r[i]);
    HyraxA::Prove(pub.ip_proofs[i], seed, input, pub.ip_com_pubs[i],
                  ip_com_sec);
  };
  parallel::For(9, parallel_f3);

  sec.a = std::move(x[5]);
  sec.b = std::move(x[6]);
  sec.c = std::move(x[7]);
  sec.d = std::move(x[8]);
  sec.r_a = com_x_r[5];
  sec.r_b = com_x_r[6];
  sec.r_c = com_x_r[7];
  sec.r_d = com_x_r[8];
}

struct PoolingR1csPub {
  clink::ParallelR1cs<typename R1cs>::Proof r1cs_proof;
  std::vector<G1> com_w;
  size_t r1cs_ret_index;

  bool operator==(PoolingR1csPub const& b) const {
    return r1cs_proof == b.r1cs_proof && com_w == b.com_w &&
           r1cs_ret_index == b.r1cs_ret_index;
  }

  bool operator!=(PoolingR1csPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.PoolingR1csPub", ("p", r1cs_proof), ("c", com_w),
                       ("r", r1cs_ret_index));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.PoolingR1csPub", ("p", r1cs_proof), ("c", com_w),
                       ("r", r1cs_ret_index));
  }
};

struct PoolingR1csSec {
  std::vector<Fr> y;
  Fr ry;
};

inline void PoolingR1csProve(h256_t seed, ProveContext const& /*context*/,
                             PoolingInputPub const& input_pub,
                             PoolingInputSec const& input_sec,
                             PoolingR1csPub& pub, PoolingR1csSec& sec) {
  Tick tick(__FN__);

  libsnark::protoboard<Fr> pb;
  circuit::vgg16::PoolingGadget<8, 24> gadget(pb, "vgg16 pooling gadget");
  int64_t const primary_input_size = 0;
  pb.set_input_sizes(primary_input_size);
  pub.r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
  std::unique_ptr<R1csInfo> r1cs_info(new R1csInfo(pb));
  auto s = r1cs_info->num_variables;
  std::vector<std::vector<Fr>> w(s);
  auto n = input_sec.a.size();
  for (auto& i : w) i.resize(n);

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
      w[i][j] = v[i];
    }
  }

  pub.com_w.resize(s);
  std::vector<Fr> com_w_r(s);

  std::cout << "compute pooling com(witness)\n";
  com_w_r[0] = input_sec.r_a;
  pub.com_w[0] = input_pub.ip_com_pubs[5].xi;
  com_w_r[1] = input_sec.r_b;
  pub.com_w[1] = input_pub.ip_com_pubs[6].xi;
  com_w_r[2] = input_sec.r_c;
  pub.com_w[2] = input_pub.ip_com_pubs[7].xi;
  com_w_r[3] = input_sec.r_d;
  pub.com_w[3] = input_pub.ip_com_pubs[8].xi;

  auto parallel_f = [&com_w_r, &pub, &w](int64_t i) {
    std::vector<Fr> const& x = w[i + 4];
    Fr& r = com_w_r[i + 4];
    G1& c = pub.com_w[i + 4];
    r = FrRand();
    c = pc::PcComputeCommitmentG(x, r, true);
  };
  parallel::For<int64_t>(s - 4, parallel_f);

  // save output
  sec.y = w[pub.r1cs_ret_index];
  sec.ry = com_w_r[pub.r1cs_ret_index];

  typename R1cs::ProveInput r1cs_input(*r1cs_info, std::move(w), pub.com_w,
                                       com_w_r, pc::kGetRefG);

  R1cs::Prove(pub.r1cs_proof, seed, std::move(r1cs_input));
}

struct PoolingOutputPub {
  std::array<HyraxA::CommitmentPub, 6> ip_com_pubs;
  std::array<HyraxA::Proof, 6> ip_proofs;
};

inline void PoolingUpdateSeed(h256_t& seed, PoolingR1csPub const& pub) {
  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(pub);
  auto buf = os.get_shared_buffer();
  HashUpdate(hash, buf.data.get(), buf.size);
  hash.Final(seed.data());
}

inline void PoolingBuildOutputQ(h256_t const& seed,
                                std::array<std::vector<Fr>, 6>& q) {
  std::vector<size_t> size;
  size_t total_size = 0;
  for (size_t i = 0; i < kLayerTypeOrders.size(); ++i) {
    if (kLayerTypeOrders[i].first != kPooling) continue;
    auto const& info = kImageInfos[i + 1];
    size.push_back(info.channel_count * info.dimension * info.dimension);
    total_size += size.back();
  }
  q[0].resize(total_size);
  for (size_t i = 0; i < q.size() - 1; ++i) {
    q[i + 1].resize(size[i]);
  }

  ComputeFst(seed, "pooling output q", q[0]);
  std::vector<Fr>::const_iterator q_cur = q[0].begin();
  for (size_t i = 1; i < q.size(); ++i) {
    q[i].assign(q_cur, q_cur + q[i].size());
    q_cur += q[i].size();
  }
}

inline void PoolingOutputProve(h256_t seed, ProveContext const& context,
                               PoolingR1csPub const& r1cs_pub,
                               PoolingR1csSec const& r1cs_sec,
                               PoolingOutputPub& pub) {
  Tick tick(__FN__);
  std::array<std::vector<Fr>, 6> x;
  std::array<Fr, 6> com_x_r;
  x[0] = r1cs_sec.y;
  pub.ip_com_pubs[0].xi = r1cs_pub.com_w[r1cs_pub.r1cs_ret_index];
  com_x_r[0] = r1cs_sec.ry;

  auto layers = PoolingGetLayers();
  for (size_t i = 0; i < layers.size(); ++i) {
    auto layer = layers[i] + 1;
    auto image = context.const_images()[layer];
    x[i + 1] = image->copy_vec();
    pub.ip_com_pubs[i + 1].xi = context.image_com_pub().c[layer];
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

#ifdef _DEBUG_CHECK
  Fr sum_xq = std::accumulate(xq.begin() + 1, xq.end(), FrZero());
  if (sum_xq != xq[0]) throw std::runtime_error("oops");
  G1 sum_com_xq = std::accumulate(com_xq.begin() + 1, com_xq.end(), G1Zero());
  if (sum_com_xq != com_xq[0]) throw std::runtime_error("oops");
#endif

  auto parallel_f2 = [&seed, &x, &q, &xq, &com_xq, &com_x_r, &com_xq_r,
                      &pub](int64_t i) {
    HyraxA::ProveInput input(x[i], q[i], xq[i], pc::kGetRefG, pc::PcG(0));
    pub.ip_com_pubs[i].tau = com_xq[i];
    HyraxA::CommitmentSec ip_com_sec(com_x_r[i], com_xq_r[i]);
    HyraxA::Prove(pub.ip_proofs[i], seed, input, pub.ip_com_pubs[i],
                  ip_com_sec);
  };
  parallel::For(6, parallel_f2);
}

struct PoolingProof {
  PoolingInputPub input;
  PoolingR1csPub r1cs;
  PoolingOutputPub output;
};

inline void PoolingProve(h256_t seed, ProveContext const& context,
                         PoolingProof& proof) {
  Tick tick(__FN__);
  PoolingInputPub input_pub;
  PoolingInputSec input_sec;
  PoolingInputProve(seed, context, input_pub, input_sec);

  PoolingR1csPub r1cs_pub;
  PoolingR1csSec r1cs_sec;
  PoolingR1csProve(seed, context, input_pub, input_sec, r1cs_pub, r1cs_sec);

  PoolingOutputPub output_pub;
  PoolingOutputProve(seed, context, r1cs_pub, r1cs_sec, output_pub);

  proof.input = input_pub;
  proof.r1cs = r1cs_pub;
  proof.output = output_pub;
}
}  // namespace clink::vgg16