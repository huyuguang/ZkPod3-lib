#pragma once

#include "./para_com.h"
#include "./para_fr.h"
#include "./para_pub.h"
#include "circuit/fixed_point/fixed_point.h"
#include "circuit/vgg16/vgg16.h"

namespace clink::vgg16 {

// input type: <D,N>
// output type: <D, 2N>
inline void InferConv(Para::ConvLayer const& layer,
                      Image const& input_image, Image& output_image) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  size_t const C = layer.channel_count();
  size_t const D = layer.dimension;
  size_t const K = layer.kernel_count();

  assert(input_image.dimension() == D);
  assert(input_image.channel_count() == C);
  assert(output_image.dimension() == D);
  assert(output_image.channel_count() == K);

  auto get_image = [&input_image](size_t h, size_t i, size_t j) -> Fr const& {
    auto d = input_image.dimension();
    if (i == 0 || j == 0 || i == (d + 1) || j == (d + 1)) return FrZero();
    return input_image.pixels[h][i - 1][j - 1];
  };

  auto DD = D * D;
  auto CDD = C * DD;
  // auto KDD = K * DD;
  auto KCDD = K * CDD;

  boost::multi_array<Fr, 2> b(boost::extents[CDD][9]);
  for (size_t i = 0; i < C * D * D; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      size_t m = j / 3;
      size_t n = j % 3;
      size_t r = i / DD;
      size_t p = i % DD;
      size_t q = p / D;
      size_t o = p % D;
      b[i][j] = get_image(r, q + m, o + n);
    }
  }

  boost::multi_array<Fr, 2> c(boost::extents[KCDD][9]);
  for (size_t i = 0; i < KCDD; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      c[i][j] = b[i % CDD][j];
    }
  }

  boost::multi_array<Fr, 2> p(boost::extents[KCDD][9]);
  for (size_t i = 0; i < KCDD; ++i) {
    size_t coef_offset = i / DD;
    size_t coef_k = coef_offset / C;
    size_t coef_c = coef_offset % C;
    for (size_t j = 0; j < 9; ++j) {
      p[i][j] = layer.coefs[coef_k][coef_c][j / 3][j % 3];
    }
  }

  std::vector<Fr> x(KCDD);

  //for (size_t i = 0; i < KCDD; ++i) {
  //  x[i] = std::inner_product(c[i].begin(), c[i].end(), p[i].begin(), FrZero());
  //}
  auto parallel_f1 = [&x, &c, &p](int64_t i) {
    x[i] = std::inner_product(c[i].begin(), c[i].end(), p[i].begin(), FrZero());
  };
  parallel::For((int64_t)KCDD, parallel_f1);

  auto& output = output_image.pixels;
  for (size_t i = 0; i < K; ++i) {
    auto bias = layer.bias[i] * fp::RationalConst<8, 24>().kFrN;
    for (size_t j = 0; j < D; ++j) {
      for (size_t k = 0; k < D; ++k) {
        for (size_t l = 0; l < C; ++l) {
          output[i][j][k] += x[i * CDD + l * DD + j * D + k];
        }
        output[i][j][k] += bias;
      }
    }
  }

#ifdef _DEBUG
  std::vector<std::vector<std::vector<Fr>>> debug_output(K);
  for (size_t i = 0; i < K; ++i) {
    debug_output[i].resize(D);
    for (size_t j = 0; j < D; ++j) {
      debug_output[i][j].resize(D);
      for (size_t k = 0; k < D; ++k) {
        debug_output[i][j][k] = output[i][j][k];
      }
    }
  }
#endif
  return;
}

// input type: <D,2N>
// output type: <D,N>
inline void InferReluBn(Para::BnLayer const& layer, Image const& input_image,
                        Image& output_image) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  for (size_t i = 0; i < input_data.size(); ++i) {
    auto const& mu = layer.mu[i];
    auto const& beta = layer.beta[i];
    auto const& alpha = layer.alpha[i];
    for (size_t j = 0; j < input_data[0].size(); ++j) {
      for (size_t k = 0; k < input_data[0][0].size(); ++k) {
        auto& in_data = input_data[i][j][k];
        auto& out_data = output_data[i][j][k];
        //std::cout << "infer alpha: " << alpha << "\n";
        //std::cout << "infer beta: " << beta << "\n";
        //std::cout << "infer mu: " << mu << "\n";

        // relu
        //std::cout << "infer a : " << in_data << "\n";

        out_data = fp::ReducePrecision<8, 24 * 2, 24>(in_data);
        //std::cout << "infer a2 : " << out_data << "\n";
        
        out_data = out_data.isNegative() ? 0 : out_data;
        //std::cout << "infer b : " << out_data << "\n";

        // bn
        out_data =
            alpha * (out_data - mu) + beta * fp::RationalConst<8, 24>().kFrN;
        //std::cout << "infer c : " << out_data << "\n";

        // NOTE: output_data[i][j][k] is <D,2N>                                 
        // now we need to convert <D,2N> to <D,N>
        out_data = fp::ReducePrecision<8, 24 * 2, 24>(out_data);
        //std::cout << "infer ret : " << out_data << "\n\n";
        //std::cout << fp::RationalToDouble<8, 24>(output_data[i][j][k]) << "\n";

#ifdef _DEBUG_CHECK
        libsnark::protoboard<Fr> pb;
        circuit::vgg16::ReluBnGadget<8, 24> gadget(pb, "ReluBnGadget");
        gadget.Assign(in_data, alpha, beta, mu);
        if (!pb.is_satisfied()) {
          std::cout << __LINE__ << " oops";
          throw std::runtime_error("oops");
        }
        Fr gadget_data = pb.val(gadget.ret());
        if (gadget_data != out_data) {
          std::cout << "\n\noops:\n";
          throw std::runtime_error("oops");
        }
#endif 
      }
    }
  }  

#ifdef _DEBUG
  std::vector<std::vector<std::vector<Fr>>> debug_data(output_data.size());
  for (size_t i = 0; i < output_data.size(); ++i) {
    debug_data[i].resize(output_data[0].size());
    for (size_t j = 0; j < output_data[0].size(); ++j) {
      debug_data[i][j].resize(output_data[0][0].size());
      for (size_t k = 0; k < output_data[0][0].size(); ++k) {
        debug_data[i][j][k] = output_data[i][j][k];
      }
    }
  }
#endif
  return;
}

// input type: <D,N>
// output type: <D,N>
inline void InferMaxPooling(Image const& input_image,
                            Image& output_image) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  auto FrDN = fp::RationalConst<8, 24>().kFrDN;
  for (size_t i = 0; i < output_data.size(); ++i) {
    for (size_t j = 0; j < output_data[0].size(); ++j) {
      for (size_t k = 0; k < output_data[0][0].size(); ++k) {
        std::array<Fr, 4> rect_fr;
        rect_fr[0] = input_data[i][j * 2][k * 2] + FrDN;
        rect_fr[1] = input_data[i][j * 2][k * 2 + 1] + FrDN;
        rect_fr[2] = input_data[i][j * 2 + 1][k * 2] + FrDN;
        rect_fr[3] = input_data[i][j * 2 + 1][k * 2 + 1] + FrDN;
        std::array<mpz_class, 4> rect_mpz;
        rect_mpz[0] = rect_fr[0].getMpz();
        rect_mpz[1] = rect_fr[1].getMpz();
        rect_mpz[2] = rect_fr[2].getMpz();
        rect_mpz[3] = rect_fr[3].getMpz();
        mpz_class max_value =
            *std::max_element(rect_mpz.begin(), rect_mpz.end());        
        output_data[i][j][k].setMpz(max_value);
        output_data[i][j][k] = output_data[i][j][k] - FrDN;            
      }
    }
  }
#ifdef _DEBUG
  std::vector<std::vector<std::vector<Fr>>> debug_data(output_data.size());
  for (size_t i = 0; i < output_data.size(); ++i) {
    debug_data[i].resize(output_data[0].size());
    for (size_t j = 0; j < output_data[0].size(); ++j) {
      debug_data[i][j].resize(output_data[0][0].size());
      for (size_t k = 0; k < output_data[0][0].size(); ++k) {
        debug_data[i][j][k] = output_data[i][j][k];
      }
    }
  }
#endif
  return;
}

// input type: <D,N>
// output type: <D,2N>
inline void InferDense(Para::DenseLayer const& layer,
                       Image const& input_image, Image& output_image) {
  Tick tick(__FN__);
  namespace fp = circuit::fp;
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  assert(input_data[0].size() == 1 && input_data[0][0].size() == 1);
  assert(output_data[0].size() == 1 && output_data[0][0].size() == 1);
  std::vector<Fr> input(input_data.size() + 1);
  for (size_t i = 0; i < input_data.size(); ++i) {
    input[i] = input_data[i][0][0];
  }
  input.back() = fp::RationalConst<8, 24>().kFrN;

  for (size_t i = 0; i < output_data.size(); ++i) {
    //for (size_t j = 0; j < input.size(); ++j) {
    //  std::cout << fp::RationalToDouble<8, 24>(input[j]) << ";";
    //  std::cout << fp::RationalToDouble<8, 24>(layer.weight[i][j]) << "\n";
    //}
    output_data[i][0][0] = std::inner_product(
        input.begin(), input.end(), layer.weight[i].begin(), FrZero());
    //std::cout << fp::RationalToDouble<8, 24 + 24>(output_data[i][0][0]) << "\n";
  }
}

inline void Infer(Para const& para, dbl::Image const& dbl_image,
                  std::array<std::unique_ptr<Image>, 35>& images) {
  Tick tick(__FN__);
  images[0].reset(new Image(dbl_image));
  for (size_t i = 1; i < images.size(); ++i) {
    images[i].reset(new Image(kImageInfos[i]));
  }

  InferConv(para.conv_layer(0), *images[0], *images[1]);

  InferReluBn(para.bn_layer(0), *images[1], *images[2]);

  InferConv(para.conv_layer(1), *images[2], *images[3]);

  InferReluBn(para.bn_layer(1), *images[3], *images[4]);

  InferMaxPooling(*images[4], *images[5]);

  InferConv(para.conv_layer(2), *images[5], *images[6]);

  InferReluBn(para.bn_layer(2), *images[6], *images[7]);

  InferConv(para.conv_layer(3), *images[7], *images[8]);
  
  InferReluBn(para.bn_layer(3), *images[8], *images[9]);

  InferMaxPooling(*images[9], *images[10]);

  InferConv(para.conv_layer(4), *images[10], *images[11]);

  InferReluBn(para.bn_layer(4), *images[11], *images[12]);

  InferConv(para.conv_layer(5), *images[12], *images[13]);

  InferReluBn(para.bn_layer(5), *images[13], *images[14]);

  InferConv(para.conv_layer(6), *images[14], *images[15]);

  InferReluBn(para.bn_layer(6), *images[15], *images[16]);

  InferMaxPooling(*images[16], *images[17]);

  InferConv(para.conv_layer(7), *images[17], *images[18]);

  InferReluBn(para.bn_layer(7), *images[18], *images[19]);

  InferConv(para.conv_layer(8), *images[19], *images[20]);

  InferReluBn(para.bn_layer(8), *images[20], *images[21]);

  InferConv(para.conv_layer(9), *images[21], *images[22]);

  InferReluBn(para.bn_layer(9), *images[22], *images[23]);

  InferMaxPooling(*images[23], *images[24]);

  InferConv(para.conv_layer(10), *images[24], *images[25]);
  
  InferReluBn(para.bn_layer(10), *images[25], *images[26]);
  
  InferConv(para.conv_layer(11), *images[26], *images[27]);

  InferReluBn(para.bn_layer(11), *images[27], *images[28]);

  InferConv(para.conv_layer(12), *images[28], *images[29]);

  InferReluBn(para.bn_layer(12), *images[29], *images[30]);

  InferMaxPooling(*images[30], *images[31]);

  InferDense(para.dense_layer(0), *images[31], *images[32]);

  InferReluBn(para.bn_layer(13), *images[32], *images[33]);

  InferDense(para.dense_layer(1), *images[33], *images[34]);

  images[34]->dump<8, 24 + 24>();
}

};  // namespace clink::vgg16