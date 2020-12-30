#pragma once

#include "../details.h"
#include "../parallel_r1cs.h"
#include "./demo.h"
#include "circuit/mnist/mnist.h"

// verifiable mnist prediction

namespace clink {

template <typename Policy>
struct Mnist {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using R1cs = typename clink::ParallelR1cs<Policy>;

  struct Para {
    std::array<std::array<Fr, 10>, 5> conv;
    std::array<std::array<Fr, 13 * 13 * 5 + 1>, 10> dense1;
    std::array<std::array<Fr, 10 + 1>, 10> dense2;
  };

  struct ParaCommitmentPub {
    // G1 all;
    std::array<G1, 10> conv;
    std::array<G1, 10> dense1;
    std::array<G1, 10> dense2;
  };

  struct ParaCommitmentSec {
    // Fr all;
    std::array<Fr, 10> conv;
    std::array<Fr, 10> dense1;
    std::array<Fr, 10> dense2;
  };

  static void ComputeParaCom(ParaCommitmentPub& com_pub,
                             ParaCommitmentSec& com_sec, Para const& para) {
    std::array<G1, 5> conv_g1;
    auto parallel_g = [&conv_g1](int64_t i) {
      conv_g1[i] = G1Zero();
      auto const& g = pc::PcG();
      for (size_t j = 0; j < 13 * 13; ++j) {
        conv_g1[i] += g[i + j * 5];
      }
    };
    parallel::For(5, parallel_g);

    auto parallel_f = [&para, &com_sec, &com_pub, &conv_g1](int64_t i) {
      com_sec.conv[i] = FrRand();
      auto get_g = [&conv_g1](size_t j) -> G1 const& {
        return j == 0 ? pc::PcH() : conv_g1[j - 1];
      };
      auto get_f = [&com_sec, i, &para](size_t j) -> Fr const& {
        return j == 0 ? com_sec.conv[i] : para.conv[j - 1][i];
      };
      com_pub.conv[i] = MultiExpBdlo12<G1>(get_g, get_f, 6);

      com_sec.dense1[i] = FrRand();
      com_pub.dense1[i] = pc::ComputeCom(13 * 13 * 5 + 1, para.dense1[i].data(),
                                         com_sec.dense1[i]);

      com_sec.dense2[i] = FrRand();
      com_pub.dense2[i] =
          pc::ComputeCom(10 + 1, para.dense2[i].data(), com_sec.dense2[i]);
    };
    parallel::For(10, parallel_f);

    // TODO: combine to one commitment
    // TODO: prove all parameters fit FixPoint<D,N>
  }

  struct ProveConvInput {
    typedef circuit::cnn::ConvGadget<6, 24, 4, 4, 3, 3> ConvGadget;
    ProveConvInput(std::array<Fr, 28 * 28> const& data, Para const& para,
                   ParaCommitmentPub const& para_com_pub,
                   ParaCommitmentSec const& para_com_sec)
        : data(data),
          para(para),
          para_com_pub(para_com_pub),
          para_com_sec(para_com_sec) {
      namespace fp = circuit::fp;
      libsnark::protoboard<Fr> pb;
      ConvGadget gadget(pb, "Mnist Gadget");
      int64_t const primary_input_size = 4 * 4;
      pb.set_input_sizes(primary_input_size);
      r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
      r1cs_info.reset(new R1csInfo(pb));
      s = r1cs_info->num_variables;
      w.resize(s);
      for (auto& i : w) i.resize(n);

      // public input
      std::array<std::array<Fr, 4 * 4>, n> public_inputs;
      for (size_t k = 0; k < 13 * 13; ++k) {
        auto& view = public_inputs[k * 5];
        size_t i = (k / 13) * 2;
        size_t j = (k % 13) * 2;

        for (size_t v = 0; v < view.size(); ++v) {
          size_t r = v / (3 + 1);
          size_t c = v % (3 + 1);
          view[v] = data[i * 28 + j + r * 28 + c];
        }

        for (size_t ii = 1; ii < 5; ++ii) {
          public_inputs[k * 5 + ii] = view;
        }
      }

      // secret input
      std::array<std::array<Fr, 3 * 3 + 1>, n> secret_inputs;
      for (size_t k = 0; k < 13 * 13; ++k) {
        for (size_t i = 0; i < 5; ++i) {
          secret_inputs[k * 5 + i] = para.conv[i];
        }
      }
#ifdef _DEBUG
      auto parallel_sec_com = [this, &secret_inputs, &para_com_sec,
                               &para_com_pub](int64_t i) {
        auto get_x = [&secret_inputs, i](int64_t j) -> Fr const& {
          return secret_inputs[j][i];
        };
        G1 com = pc::ComputeCom(n, get_x, para_com_sec.conv[i]);
        assert(com == para_com_pub.conv[i]);
      };
      parallel::For(3 * 3 + 1, parallel_sec_com);
#endif

      for (int64_t j = 0; j < n; ++j) {
        gadget.Assign(public_inputs[j], secret_inputs[j]);
        assert(pb.is_satisfied());
        auto v = pb.full_variable_assignment();
        for (int64_t i = 0; i < s; ++i) {
          w[i][j] = v[i];
        }
        for (int64_t i = 0; i < 4 * 4; ++i) {
          assert(w[i][j] == public_inputs[j][i]);
        }
        for (int64_t i = 0; i < 3 * 3 + 1; ++i) {
          assert(w[i + 4 * 4][j] == secret_inputs[j][i]);
        }
      }
    }
    std::array<Fr, 28 * 28> const& data;
    Para const& para;
    ParaCommitmentPub const& para_com_pub;
    ParaCommitmentSec const& para_com_sec;
    static int64_t constexpr n = 13 * 13 * 5;

    std::unique_ptr<R1csInfo> r1cs_info;
    int64_t s;
    std::vector<std::vector<Fr>> mutable w;
    int64_t r1cs_ret_index;
  };

  struct ConvProof {
    typename R1cs::Proof r1cs_proof;
    std::vector<G1> com_w;

    bool operator==(ConvProof const& b) const {
      return r1cs_proof == b.r1cs_proof && com_w == b.com_w;
    }

    bool operator!=(ConvProof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("conv.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("conv.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
  };

  struct ProveOutput {
    Fr com_r;
    G1 com;
    std::vector<Fr> data;
  };

  static void ProveConv(ConvProof& proof, ProveOutput& output, h256_t seed,
                        ProveConvInput const& input) {
    Tick tick(__FN__);
    namespace fp = circuit::fp;
    std::vector<G1> com_w(input.s);
    std::vector<Fr> com_w_r(input.s);

    std::cout << "compute com(witness)\n";
    auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
      if (i < 4 * 4) {
        com_w_r[i] = FrZero();
        com_w[i] = pc::ComputeCom(input.w[i], com_w_r[i], true);
      } else if (i >= 4 * 4 && i < 4 * 4 + 3 * 3 + 1) {
        com_w[i] = input.para_com_pub.conv[i - 4 * 4];
        com_w_r[i] = input.para_com_sec.conv[i - 4 * 4];
      } else {
        com_w_r[i] = FrRand();
        com_w[i] = pc::ComputeCom(input.w[i], com_w_r[i], true);
      }
    };
    parallel::For<int64_t>(input.s, parallel_f);

    // save output for next(dense)
    output.com_r = com_w_r[input.r1cs_ret_index];
    output.com = com_w[input.r1cs_ret_index];
    output.data = input.w[input.r1cs_ret_index];

#ifdef _DEBUG
    for (size_t j = 0; j < output.data.size(); ++j) {
      if (j % 5 == 0) {
        std::cout << std::left << std::setw(8) << std::setfill(' ') << j / 5;
      }
      double dret = fp::RationalToDouble<6, 24>(output.data[j]);
      std::cout << std::right << std::setw(12) << std::setfill(' ') << dret;
      if ((j + 1) % 5 == 0) std::cout << "\n";
    }
#endif

    typename R1cs::ProveInput r1cs_input(*input.r1cs_info, "mnist",
                                         std::move(input.w), com_w, com_w_r,
                                         pc::kGetRefG1);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    proof.com_w = std::move(com_w);
  }

  struct VerifyConvInput {
    VerifyConvInput(std::array<Fr, 28 * 28> const& data,
                    ParaCommitmentPub const& para_com_pub)
        : data(data), para_com_pub(para_com_pub) {
      typedef circuit::cnn::ConvGadget<6, 24, 4, 4, 3, 3> ConvGadget;
      libsnark::protoboard<Fr> pb;
      ConvGadget gadget(pb, "Mnist Gadget");
      int64_t const primary_input_size = 4 * 4;
      pb.set_input_sizes(primary_input_size);
      // see FieldT protoboard<FieldT>::val
      r1cs_ret_index = gadget.ret().index - 1;
      r1cs_info.reset(new R1csInfo(pb));
      m = r1cs_info->num_constraints;
      s = r1cs_info->num_variables;

      // public input
      std::array<std::array<Fr, 4 * 4>, n> public_inputs;
      for (size_t k = 0; k < 13 * 13; ++k) {
        auto& view = public_inputs[k * 5];
        size_t i = (k / 13) * 2;
        size_t j = (k % 13) * 2;

        for (size_t v = 0; v < view.size(); ++v) {
          size_t r = v / (3 + 1);
          size_t c = v % (3 + 1);
          view[v] = data[i * 28 + j + r * 28 + c];
        }

        for (size_t ii = 1; ii < 5; ++ii) {
          public_inputs[k * 5 + ii] = view;
        }
      }

      // public_w
      public_w.resize(4 * 4);
      for (size_t k = 0; k < 4 * 4; ++k) {
        public_w[k].resize(n);
        for (size_t i = 0; i < n; ++i) {
          public_w[k][i] = public_inputs[i][k];
        }
      }
    }
    static constexpr int64_t n = 13 * 13 * 5;
    std::array<Fr, 28 * 28> const& data;
    ParaCommitmentPub const& para_com_pub;
    std::unique_ptr<R1csInfo> r1cs_info;
    size_t r1cs_ret_index;
    int64_t m;
    int64_t s;
    std::vector<std::vector<Fr>> public_w;
  };

  static bool VerifyConv(ConvProof const& proof, h256_t seed,
                         VerifyConvInput const& input) {
    Tick tick(__FN__);
    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }

    // Check com of secret input
    for (size_t i = 0; i < 3 * 3 + 1; ++i) {
      if (proof.com_w[i + 4 * 4] != input.para_com_pub.conv[i]) {
        assert(false);
        return false;
      }
    }

    // Not need to check com of public input since R1cs::VerifyConv will check
    // it.

    typename R1cs::VerifyInput pr_input(input.n, *input.r1cs_info, "mnist",
                                        proof.com_w, input.public_w,
                                        pc::kGetRefG1);
    return R1cs::Verify(proof.r1cs_proof, seed, pr_input);
  }

  template <size_t M, size_t N>
  struct ProveDenseInput {
    ProveDenseInput(std::array<std::array<Fr, M + 1>, N> const& para_dense,
                    std::array<G1, N> const& com_para_dense,
                    std::array<Fr, N> const& com_para_dense_r,
                    ProveOutput&& last_output)
        : para_dense(para_dense),
          com_para_dense(com_para_dense),
          com_para_dense_r(com_para_dense_r),
          com_data_r(last_output.com_r),
          com_data(last_output.com),
          data(std::move(last_output.data)) {
      namespace fp = circuit::fp;
      CHECK(data.size() == M, "");
      this->data.push_back(fp::RationalConst<6, 24>().kFrN);
      this->com_data += pc::PcG(M) * data.back();
    }
    std::array<std::array<Fr, M + 1>, N> const& para_dense;
    std::array<G1, N> const& com_para_dense;
    std::array<Fr, N> const& com_para_dense_r;
    Fr com_data_r;
    G1 com_data;
    std::vector<Fr> data;
  };

  struct DenseProof {
    G1 com;
    G1 com_z;
    groth09::Sec51a::Proof proof_51;
    hyrax::A2::Proof proof_hy;

    bool operator==(DenseProof const& b) const {
      return com == b.com && com_z == b.com_z && proof_51 == b.proof_51 &&
             proof_hy == b.proof_hy;
    }

    bool operator!=(DenseProof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("d.p", ("cd", com), ("cz", com_z), ("51", proof_51),
                         ("hy", proof_hy));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("d.p", ("cd", com), ("cz", com_z), ("51", proof_51),
                         ("hy", proof_hy));
    }
  };

  template <size_t M, size_t N>
  static std::vector<Fr> ComputeDenseFst(h256_t const& seed) {
    std::vector<Fr> x(N);
    std::string str = "dense";
    str += std::to_string(M);
    str += "_";
    str += std::to_string(N);
    ComputeFst(seed, str, x);
    return x;
  }

  template <size_t M, size_t N>
  static void ProveDense(DenseProof& proof, ProveOutput& output, h256_t seed,
                         ProveDenseInput<M, N> const& input) {
    Tick tick(__FN__);
    namespace fp = circuit::fp;
    assert(input.para_dense.size() == N);
    assert(input.data.size() == M + 1);

    if (M == 10) {
      std::cout << "data:\n";
      for (auto const& i : input.data) {
        double ret = fp::RationalToDouble<6, 24>(i);
        std::cout << std::right << std::setw(12) << std::setfill(' ') << ret;
      }
      std::cout << "\n";
    }

    // build output
    output.data.resize(N);
    auto parallel_f = [&input, &output](int64_t i) {
      assert(input.para_dense[i].size() == M + 1);
      output.data[i] =
          std::inner_product(input.data.begin(), input.data.end(),
                             input.para_dense[i].begin(), FrZero());
    };
    parallel::For(N, parallel_f);

    //#ifdef _DEBUG
    std::cout << "dense before relu:\n";
    for (size_t i = 0; i < N; ++i) {
      double ret = fp::RationalToDouble<6, 24 + 24>(output.data[i]);
      std::cout << std::right << std::setw(12) << std::setfill(' ') << ret;
    }
    std::cout << "\n";
    //#endif

    output.com_r = FrRand();
    output.com = pc::ComputeCom(output.data, output.com_r);

    // prove
    std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
    G1 com_e = G1Zero();
    Fr com_e_r = FrZero();
    // TODO: change to multiexp and inner_product
    for (size_t i = 0; i < N; ++i) {
      com_e += input.com_para_dense[i] * x[i];
      com_e_r += input.com_para_dense_r[i] * x[i];
    }
    std::vector<Fr> e(M + 1);
    for (size_t i = 0; i < M + 1; ++i) {
      e[i] = FrZero();
      for (size_t j = 0; j < N; ++j) {
        e[i] += input.para_dense[j][i] * x[j];
      }
    }

    Fr z =
        std::inner_product(x.begin(), x.end(), output.data.begin(), FrZero());
    Fr com_z_r = FrRand();
    G1 com_z = pc::ComputeCom(z, com_z_r);
    proof.com_z = com_z;
    proof.com = output.com;

#ifdef _DEBUG
    assert(com_e == pc::ComputeCom(e, com_e_r));
    assert(z == std::inner_product(e.begin(), e.end(), input.data.begin(),
                                   FrZero()));
#endif

    // prove left
    hyrax::A2::ProveInput input_hy("mnist", output.data, x, z, pc::kGetRefG1,
                                   pc::PcG(0));
    hyrax::A2::CommitmentPub com_pub_hy(output.com, com_z);
    hyrax::A2::CommitmentSec com_sec_hy(output.com_r, com_z_r);
    hyrax::A2::Prove(proof.proof_hy, seed, input_hy, com_pub_hy, com_sec_hy);

    // prove right
    // std::cout << "prove, seed: " << misc::HexToStr(seed) << "\n";
    // std::cout << "prove, com_e: " << com_e << "\n";
    // std::cout << "prove, com_data: " << input.com_data << "\n";
    // std::cout << "prove, com_z: " << com_z << "\n";
    groth09::Sec51a::ProveInput input_51(e, input.data, z, pc::kGetRefG1,
                                         pc::kGetRefG1, pc::PcG(0));
    groth09::Sec51a::CommitmentPub com_pub_51(com_e, input.com_data, com_z);
    groth09::Sec51a::CommitmentSec com_sec_51(com_e_r, input.com_data_r,
                                              com_z_r);
    groth09::Sec51a::Prove(proof.proof_51, seed, input_51, com_pub_51,
                           com_sec_51);
  }

  template <size_t M, size_t N>
  struct VerifyDenseInput {
    VerifyDenseInput(G1 const& com_data,
                     std::array<G1, N> const& com_para_dense)
        : com_data(com_data), com_para_dense(com_para_dense) {
      namespace fp = circuit::fp;
      this->com_data += pc::PcG(M) * fp::RationalConst<6, 24>().kFrN;
    }
    G1 com_data;
    std::array<G1, N> const& com_para_dense;
  };

  template <size_t M, size_t N>
  static bool VerifyDense(DenseProof const& proof, h256_t seed,
                          VerifyDenseInput<M, N> const& input) {
    Tick tick(__FN__);
    std::vector<Fr> x = ComputeDenseFst<M, N>(seed);
    G1 com_e = G1Zero();
    // TODO: change to multiexp
    for (size_t i = 0; i < N; ++i) {
      com_e += input.com_para_dense[i] * x[i];
    }

    hyrax::A2::CommitmentPub com_pub_hy(proof.com, proof.com_z);
    hyrax::A2::VerifyInput input_hy("mnist", x, com_pub_hy, pc::kGetRefG1,
                                    pc::PcG(0));
    if (!hyrax::A2::Verify(proof.proof_hy, seed, input_hy)) {
      assert(false);
      return false;
    }

    // std::cout << "verify, seed: " << misc::HexToStr(seed) << "\n";
    // std::cout << "verify, com_e: " << com_e << "\n";
    // std::cout << "verify, com: " << input.com << "\n";
    // std::cout << "verify, com_z: " << proof.com_z << "\n";
    groth09::Sec51a::CommitmentPub com_pub_51(com_e, input.com_data,
                                              proof.com_z);
    groth09::Sec51a::VerifyInput input_51(com_pub_51, pc::kGetRefG1,
                                          pc::kGetRefG1, pc::PcG(0));
    if (!groth09::Sec51a::Verify(proof.proof_51, seed, input_51)) {
      assert(false);
      return false;
    }

    return true;
  }

  // the dense output is FixedPoint<D, 2N>
  template <size_t N>
  struct ProveRelu2Input {
    typedef circuit::fixed_point::Relu2Gadget<6, 24 * 2, 24> Relu2Gadget;
    ProveRelu2Input(ProveOutput&& last_output)
        : data(std::move(last_output.data)),
          com(last_output.com),
          com_r(last_output.com_r) {
      namespace fp = circuit::fp;
      assert(data.size() == N);

      libsnark::protoboard<Fr> pb;
      Relu2Gadget gadget(pb, "Mnist Relu2Gadget");
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);
      r1cs_ret_index = gadget.ret().index - 1;  // see protoboard<FieldT>::val
      r1cs_info.reset(new R1csInfo(pb));
      s = r1cs_info->num_variables;
      w.resize(s);
      for (auto& i : w) i.resize(N);

#ifdef _DEBUG
      assert(com == pc::ComputeCom(data, com_r));
#endif

      for (int64_t j = 0; j < (int64_t)N; ++j) {
        gadget.Assign(data[j]);
        assert(pb.is_satisfied());
        auto v = pb.full_variable_assignment();
        for (int64_t i = 0; i < s; ++i) {
          w[i][j] = v[i];
        }
        assert(w[0][j] == data[j]);
      }
    }

    std::vector<Fr> data;
    G1 com;
    Fr com_r;

    std::unique_ptr<R1csInfo> r1cs_info;
    int64_t s;
    std::vector<std::vector<Fr>> mutable w;
    int64_t r1cs_ret_index;
  };

  struct Relu2Proof {
    typename R1cs::Proof r1cs_proof;
    std::vector<G1> com_w;
    bool operator==(Relu2Proof const& b) const {
      return r1cs_proof == b.r1cs_proof && com_w == b.com_w;
    }

    bool operator!=(Relu2Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("dr.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("dr.p", ("r1cs", r1cs_proof), ("w", com_w));
    }
  };

  template <size_t N>
  static void ProveRelu2(Relu2Proof& proof, ProveOutput& output, h256_t seed,
                         ProveRelu2Input<N> const& input) {
    Tick tick(__FN__);
    namespace fp = circuit::fp;
    std::vector<G1> com_w(input.s);
    std::vector<Fr> com_w_r(input.s);

    std::cout << "compute com(witness)\n";
    auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
      if (i == 0) {
        com_w_r[i] = input.com_r;
        com_w[i] = input.com;
      } else {
        com_w_r[i] = FrRand();
        com_w[i] = pc::ComputeCom(input.w[i], com_w_r[i], true);
      }
    };
    parallel::For<int64_t>(input.s, parallel_f);

    // save output for next step
    output.com_r = com_w_r[input.r1cs_ret_index];
    output.com = com_w[input.r1cs_ret_index];
    output.data = input.w[input.r1cs_ret_index];

#ifdef _DEBUG
    std::cout << "dense relu:\n";
    for (size_t j = 0; j < output.data.size(); ++j) {
      double dret = fp::RationalToDouble<6, 24>(output.data[j]);
      std::cout << std::right << std::setw(12) << std::setfill(' ') << dret;
    }
    std::cout << "\n";
#endif

    // prove
    typename R1cs::ProveInput r1cs_input(*input.r1cs_info, "mnist",
                                         std::move(input.w), com_w, com_w_r,
                                         pc::kGetRefG1);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    proof.com_w = std::move(com_w);
  }

  struct VerifyRelu2Input {
    VerifyRelu2Input(G1 const& com) : com(com) {
      typedef circuit::fixed_point::Relu2Gadget<6, 24 * 2, 24> Relu2Gadget;
      libsnark::protoboard<Fr> pb;
      Relu2Gadget gadget(pb, "Mnist Relu2Gadget");
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);
      r1cs_ret_index =
          gadget.ret().index - 1;  // see FieldT protoboard<FieldT>::val
      r1cs_info.reset(new R1csInfo(pb));
      m = r1cs_info->num_constraints;
      s = r1cs_info->num_variables;
    }
    static constexpr int64_t n = 10;
    G1 const& com;
    std::unique_ptr<R1csInfo> r1cs_info;
    size_t r1cs_ret_index;
    int64_t m;
    int64_t s;
    std::vector<std::vector<Fr>> public_w;  // empty
  };

  static bool VerifyRelu2(Relu2Proof const& proof, h256_t seed,
                          VerifyRelu2Input const& input) {
    Tick tick(__FN__);
    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }

    // Check com of secret input
    if (proof.com_w[0] != input.com) {
      assert(false);
      return false;
    }

    typename R1cs::VerifyInput pr_input(input.n, *input.r1cs_info, "mnist",
                                        proof.com_w, input.public_w,
                                        pc::kGetRefG1);
    return R1cs::Verify(proof.r1cs_proof, seed, pr_input);
  }

  template <typename ProofT>
  static void UpdateSeed(h256_t& seed, ProofT const& proof) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(proof);
    auto buf = os.get_shared_buffer();
    HashUpdate(hash, buf.data.get(), buf.size);
    hash.Final(seed.data());
  }

  struct Proof {
    ConvProof conv;          // conv,relu,flatten
    DenseProof dense1;       // dense1
    Relu2Proof dense1_relu;  // relu2 for dense1
    DenseProof dense2;       // dense2

    bool operator==(Proof const& b) const {
      return conv == b.conv && dense1 == b.dense1 &&
             dense1_relu == b.dense1_relu && dense2 == b.dense2;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("mnist.p", ("c", conv), ("d1", dense1),
                         ("dr", dense1_relu), ("d2", dense2));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("mnist.p", ("c", conv), ("d1", dense1),
                         ("dr", dense1_relu), ("d2", dense2));
    }
  };

  static void Prove(h256_t seed, Proof& proof,
                    std::array<Fr, 28 * 28> const& data, Para const& para,
                    ParaCommitmentPub const& para_com_pub,
                    ParaCommitmentSec const& para_com_sec) {
    Tick tick(__FN__);
    // prove conv
    ProveConvInput conv_input(data, para, para_com_pub, para_com_sec);
    ProveOutput conv_output;
    ProveConv(proof.conv, conv_output, seed, conv_input);
    UpdateSeed(seed, proof.conv);

    // prove dense1
    ProveDenseInput<13 * 13 * 5, 10> dense1_input(
        para.dense1, para_com_pub.dense1, para_com_sec.dense1,
        std::move(conv_output));
    ProveOutput dense1_output;
    ProveDense<13 * 13 * 5, 10>(proof.dense1, dense1_output, seed,
                                dense1_input);
    UpdateSeed(seed, proof.dense1);

    // prove dense1 relu
    ProveRelu2Input<10> dense1_relu_input(std::move(dense1_output));
    ProveOutput dense1_relu_output;
    ProveRelu2<10>(proof.dense1_relu, dense1_relu_output, seed,
                   dense1_relu_input);
    UpdateSeed(seed, proof.dense1_relu);

    // prove dense2
    ProveDenseInput<10, 10> dense2_input(para.dense2, para_com_pub.dense2,
                                         para_com_sec.dense2,
                                         std::move(dense1_relu_output));
    ProveOutput dense2_output;
    ProveDense<10, 10>(proof.dense2, dense2_output, seed, dense2_input);
  }

  static bool Verify(h256_t seed, Proof const& proof,
                     std::array<Fr, 28 * 28> const& data,
                     ParaCommitmentPub const& para_com_pub) {
    Tick tick(__FN__);
    // verify conv
    VerifyConvInput conv_input(data, para_com_pub);
    if (!VerifyConv(proof.conv, seed, conv_input)) {
      assert(false);
      return false;
    }
    UpdateSeed(seed, proof.conv);

    // verify dense1
    VerifyDenseInput<13 * 13 * 5, 10> dense1_input(
        proof.conv.com_w[conv_input.r1cs_ret_index], para_com_pub.dense1);
    if (!VerifyDense<13 * 13 * 5, 10>(proof.dense1, seed, dense1_input)) {
      assert(false);
      return false;
    }
    UpdateSeed(seed, proof.dense1);

    // verify dense1 relu
    VerifyRelu2Input dense1_relu_input(proof.dense1.com);
    if (!VerifyRelu2(proof.dense1_relu, seed, dense1_relu_input)) {
      assert(false);
      return false;
    }
    UpdateSeed(seed, proof.dense1_relu);

    // verify dense2
    VerifyDenseInput<10, 10> dense2_input(
        proof.dense1_relu.com_w[dense1_relu_input.r1cs_ret_index],
        para_com_pub.dense2);
    if (!VerifyDense<10, 10>(proof.dense2, seed, dense2_input)) {
      assert(false);
      return false;
    }

    // now we need to pod the data which commit by proof.dense2.com
    return true;
  }

  static void LoadPara(mnist::dbl::Para const& dbl_para, Para& para) {
    namespace fp = circuit::fp;
    for (size_t k = 0; k < dbl_para.conv.size(); ++k) {
      for (size_t i = 0; i < 9; ++i) {
        para.conv[k][i] = fp::DoubleToRational<6, 24>(dbl_para.conv[k].coef[i]);
      }
      para.conv[k].back() = fp::DoubleToRational<6, 24>(dbl_para.conv[k].bias);
    }
    for (size_t k = 0; k < dbl_para.dense1.size(); ++k) {
      for (size_t i = 0; i < 13 * 13 * 5; ++i) {
        para.dense1[k][i] =
            fp::DoubleToRational<6, 24>(dbl_para.dense1[k].coef[i]);
      }
      para.dense1[k].back() =
          fp::DoubleToRational<6, 24>(dbl_para.dense1[k].bias);
    }
    for (size_t k = 0; k < dbl_para.dense2.size(); ++k) {
      for (size_t i = 0; i < 10; ++i) {
        para.dense2[k][i] =
            fp::DoubleToRational<6, 24>(dbl_para.dense2[k].coef[i]);
      }
      para.dense2[k].back() =
          fp::DoubleToRational<6, 24>(dbl_para.dense2[k].bias);
    }
  }

  // convert double to fr
  static void LoadImage(mnist::dbl::UniData const& dbl_uni_data,
                        std::array<Fr, 28 * 28>& data) {
    namespace fp = circuit::fp;
    for (size_t i = 0; i < 28 * 28; ++i) {
      data[i] = fp::DoubleToRational<6, 24>(dbl_uni_data[i]);
    }
  }

  static bool Test();

  static bool TestSerialize(Proof const& proof);
};

template <typename Policy>
bool Mnist<Policy>::TestSerialize(Proof const& proof) {
  {
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(proof.conv);
    std::cout << "proof conv size: " << os.get_shared_buffer().size << "\n";
  }
  {
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(proof.dense1);
    std::cout << "proof dense1 size: " << os.get_shared_buffer().size << "\n";
  }
  {
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(proof.dense1_relu);
    std::cout << "proof dense1_relu size: " << os.get_shared_buffer().size
              << "\n";
  }
  {
    yas::mem_ostream os;
    yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
    oa.serialize(proof.dense2);
    std::cout << "proof dense2 size: " << os.get_shared_buffer().size << "\n";
  }

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
  return true;
}

template <typename Policy>
bool Mnist<Policy>::Test() {
  std::unique_ptr<Para> para(new Para);
  LoadPara(mnist::demo::GetDblPara(), *para);

  std::array<Fr, 28 * 28> data;
  LoadImage(mnist::demo::GetDblUniData(), data);

  // compute com of para
  std::unique_ptr<ParaCommitmentPub> para_com_pub(new ParaCommitmentPub);
  std::unique_ptr<ParaCommitmentSec> para_com_sec(new ParaCommitmentSec);
  ComputeParaCom(*para_com_pub, *para_com_sec, *para);

  Tick tick(__FN__);

  auto seed = misc::RandH256();

  Proof proof;
  Prove(seed, proof, data, *para, *para_com_pub, *para_com_sec);

  bool success = Verify(seed, proof, data, *para_com_pub);

#ifndef DISABLE_SERIALIZE_CHECK
  success = success && TestSerialize(proof);
#endif

  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";

  return success;
}
}  // namespace clink