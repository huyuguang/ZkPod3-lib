#pragma once

#include "para_pub.h"
#include "public.h"

namespace clink {
namespace vgg16 {
namespace dbl {
class Para {
 public:
  Para(std::string const& path) : path_(path) {
    Tick tick(__FN__);
    if (!LoadConvLayers()) throw std::runtime_error("load conv para failed");
    if (!LoadBnLayers()) throw std::runtime_error("load bn para failed");
    if (!LoadDenseLayers()) throw std::runtime_error("load dense para failed");
  }

  struct ConvLayer {
    ConvLayer(ConvLayerInfo const& info)
        : order(info.order), dimension(info.dimension) {
      coefs.resize(info.kernel_count);
      for (auto& i : coefs) {
        i.resize(info.channel_count);
      }
      bias.resize(info.kernel_count);
    }
    size_t const order;
    size_t const dimension;
    // size=K*C*3*3
    std::vector<std::vector<std::array<std::array<double, 3>, 3>>> coefs;
    std::vector<double> bias;  // size=K
    size_t kernel_count() const { return bias.size(); }
    size_t channel_count() const { return coefs[0].size(); }
  };

  struct BnLayer {
    BnLayer(BnLayerInfo const& info) : order(info.order) {
      mu.resize(info.channel_count);
      alpha.resize(info.channel_count);
      beta.resize(info.channel_count);
    }
    size_t const order;
    std::vector<double> mu;
    std::vector<double> alpha;
    std::vector<double> beta;
  };

  struct DenseLayer {
    DenseLayer(DenseLayerInfo const& info) : order(info.order) {
      weight.resize(info.type_count);
      for (auto& i : weight) {
        i.resize(info.input_count + 1);
      }
    }
    size_t const order;
    std::vector<std::vector<double>> weight;
  };

  BnLayer const& bn_layer(size_t order) const { return *bn_layers_[order]; }

  ConvLayer const& conv_layer(size_t order) const {
    return *conv_layers_[order];
  }

  DenseLayer const& dense_layer(size_t order) const {
    return *dense_layers_[order];
  }

 private:
  bool LoadConvLayers() {
    Tick tick(__FN__);
    for (size_t i = 0; i < conv_layers_.size(); ++i) {
      conv_layers_[i].reset(new ConvLayer(kConvLayerInfos[i]));
      if (!LoadConvLayer(*conv_layers_[i])) {
        return false;
      }
    }
    return true;
  }

  bool LoadConvLayer(ConvLayer& layer) {
    std::string order_str = std::to_string(layer.order + 1);
    std::string conv_name = path_ + "/features_conv_" + order_str + "/" +
                            "conv" + order_str + ".txt";
    std::ifstream conv_file(conv_name);
    std::string bias_name =
        path_ + "/features_conv_" + order_str + "/" + "bias.txt";
    std::ifstream bias_file(bias_name);

    std::string conv_line;
    std::string bias_line;
    for (size_t k = 0; k < layer.bias.size(); ++k) {
      auto& bias = layer.bias[k];
      auto& coefs = layer.coefs[k];
      if (!std::getline(conv_file, conv_line)) {
        assert(false);
        return false;
      }
      if (conv_line != "conv" + order_str + "_ " + std::to_string(k)) {
        assert(false);
        return false;
      }

      for (size_t c = 0; c < coefs.size(); ++c) {
        auto& coef = coefs[c];

        for (size_t i = 0; i < 3; ++i) {
          if (!std::getline(conv_file, conv_line)) {
            assert(false);
            return false;
          }
          for (auto& s : conv_line) {
            if (s == '\t' || s == '[' || s == ']' || s == '\r' || s == '\n')
              s = ' ';
          }
          std::istringstream iss(conv_line);
          if (!(iss >> coef[i][0] >> coef[i][1] >> coef[i][2])) {
            assert(false);
            return false;
          }
        }
      }

      if (!std::getline(bias_file, bias_line)) {
        assert(false);
        return false;
      }
      std::istringstream iss(bias_line);
      if (!(iss >> bias)) {
        assert(false);
        return false;
      }
    }

    return true;
  }

  bool LoadBnLayers() {
    Tick tick(__FN__);
    for (size_t i = 0; i < bn_layers_.size(); ++i) {
      bn_layers_[i].reset(new BnLayer(kBnLayerInfos[i]));
      if (!LoadBnLayer(*bn_layers_[i])) {
        return false;
      }
    }
    return true;
  }

  bool LoadBnLayer(BnLayer& layer) {
    std::string order_str = std::to_string(layer.order + 1);

    std::string beta_name = path_ + "/features_bn_" + order_str + "/" + "bn" +
                            order_str + "_beta.txt";
    std::ifstream beta_file(beta_name);

    std::string miu_name = path_ + "/features_bn_" + order_str + "/" + "bn" +
                           order_str + "_miu.txt";
    std::ifstream miu_file(miu_name);

    std::string sigma_name = path_ + "/features_bn_" + order_str + "/" + "bn" +
                             order_str + "_sigma_process.txt";
    std::ifstream sigma_file(sigma_name);

    std::string gamma_name = path_ + "/features_bn_" + order_str + "/" + "bn" +
                             order_str + "_gamma.txt";
    std::ifstream gamma_file(gamma_name);

    auto read = [](std::ifstream& ifs, double& v) {
      std::string line;
      if (!std::getline(ifs, line)) return false;

      std::istringstream iss(line);
      iss >> v;
      return true;
    };

    for (size_t c = 0; c < layer.beta.size(); ++c) {
      if (!read(beta_file, layer.beta[c])) {
        assert(false);
        return false;
      }
      if (!read(miu_file, layer.mu[c])) {
        assert(false);
        return false;
      }
      double gamma, sigma;
      if (!read(gamma_file, gamma)) {
        assert(false);
        return false;
      }
      if (!read(sigma_file, sigma)) {
        assert(false);
        return false;
      }
      layer.alpha[c] = gamma * sigma;
    }
    return true;
  }

  bool LoadDenseLayers() {
    Tick tick(__FN__);
    for (size_t i = 0; i < dense_layers_.size(); ++i) {
      dense_layers_[i].reset(new DenseLayer(kDenseLayerInfos[i]));
      if (!LoadDenseLayer(*dense_layers_[i])) {
        return false;
      }
    }
    return true;
  }

  bool LoadDenseLayer(DenseLayer& layer) {
    std::string order_str = std::to_string(layer.order + 1);

    std::string bias_name =
        path_ + "/features_dense_" + order_str + "/" + "bias.txt";
    std::ifstream bias_file(bias_name);

    std::string weight_name =
        path_ + "/features_dense_" + order_str + "/" + "weights.txt";
    std::ifstream weight_file(weight_name);

    for (size_t c = 0; c < layer.weight.size(); ++c) {
      std::string line;
      if (!std::getline(bias_file, line)) {
        assert(false);
        return false;
      }

      std::istringstream iss(line);
      iss >> layer.weight[c].back();
    }

    for (size_t c = 0; c < layer.weight[0].size() - 1; ++c) {
      std::string line;
      if (!std::getline(weight_file, line)) {
        assert(false);
        return false;
      }

      for (auto& s : line) {
        if (s == '\t' || s == ',' || s == ';' || s == '\r' || s == '\n')
          s = ' ';
      }

      std::istringstream iss(line);
      for (size_t d = 0; d < layer.weight.size(); ++d) {
        iss >> layer.weight[d][c];
      }
    }
    return true;
  }

 private:
  std::string const path_;
  std::array<std::unique_ptr<BnLayer>, 14> bn_layers_;
  std::array<std::unique_ptr<ConvLayer>, 13> conv_layers_;
  std::array<std::unique_ptr<DenseLayer>, 2> dense_layers_;
};

struct Image {
  size_t const order;
  boost::multi_array<double, 3> pixels;

  Image(ImageInfo const& info)
      : order(info.order),
        pixels(boost::extents[info.channel_count][info.dimension]
                             [info.dimension]) {
  }

  size_t channel_count() const { return pixels.size(); }
  size_t dimension() const { return pixels[0].size(); }
  void dump() const {
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          std::cout << pixels[i][j][k] << ";";
        }
      }
    }
    std::cout << "\n\n";
  }
};

inline bool LoadTestImage(std::string const& path, Image& image) {
  assert(image.channel_count() == 3);
  std::string line;
  for (size_t c = 0; c < image.channel_count(); ++c) {
    auto name = path + "/test_image_" + std::to_string(c) + ".txt";
    std::ifstream file(name);
    for (size_t i = 0; i < image.dimension(); ++i) {
      if (!std::getline(file, line)) {
        assert(false);
        return false;
      }

      for (auto& s : line) {
        if (s == ',') s = ' ';
      }

      std::istringstream iss(line);
      for (size_t j = 0; j < image.dimension(); ++j) {
        if (!(iss >> image.pixels[c][i][j])) {
          assert(false);
          return false;
        }
      }
    }
  }
  return true;
}

inline void InferConv(Para::ConvLayer const& layer,
                      Image const& input_image, Image& output_image) {
  size_t const C = layer.channel_count();
  size_t const D = layer.dimension;
  size_t const K = layer.kernel_count();

  assert(input_image.dimension() == D);
  assert(input_image.channel_count() == C);
  assert(output_image.dimension() == D);
  assert(output_image.channel_count() == K);

  auto get_image = [&input_image](size_t h, size_t i, size_t j) -> double {
    auto d = input_image.dimension();
    if (i == 0 || j == 0 || i == (d + 1) || j == (d + 1)) return 0;
    return input_image.pixels[h][i - 1][j - 1];
  };

  auto DD = D * D;
  auto CDD = C * DD;
  // auto KDD = K * DD;
  auto KCDD = K * CDD;

  boost::multi_array<double, 2> b(boost::extents[CDD][9]);
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

  boost::multi_array<double, 2> c(boost::extents[KCDD][9]);
  for (size_t i = 0; i < KCDD; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      c[i][j] = b[i % CDD][j];
    }
  }

  boost::multi_array<double, 2> p(boost::extents[KCDD][9]);
  for (size_t i = 0; i < KCDD; ++i) {
    size_t coef_offset = i / DD;
    size_t coef_k = coef_offset / C;
    size_t coef_c = coef_offset % C;
    for (size_t j = 0; j < 9; ++j) {
      p[i][j] = layer.coefs[coef_k][coef_c][j / 3][j % 3];
    }
  }

  std::vector<double> x(KCDD);
  for (size_t i = 0; i < KCDD; ++i) {
    x[i] = std::inner_product(c[i].begin(), c[i].end(), p[i].begin(), 0.0);
  }

  auto& output = output_image.pixels;
  for (size_t i = 0; i < K; ++i) {
    for (size_t j = 0; j < D; ++j) {
      for (size_t k = 0; k < D; ++k) {
        for (size_t l = 0; l < C; ++l) {
          output[i][j][k] += x[i * CDD + l * DD + j * D + k];
        }
        output[i][j][k] += layer.bias[i];
      }
    }
  }

#ifdef _DEBUG
  std::vector<std::vector<std::vector<double>>> debug_output(K);
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

inline void InferReluBn(Para::BnLayer const& layer, Image const& input_image,
                        Image& output_image) {
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  for (size_t i = 0; i < input_data.size(); ++i) {
    for (size_t j = 0; j < input_data[0].size(); ++j) {
      for (size_t k = 0; k < input_data[0][0].size(); ++k) {
        output_data[i][j][k] =
            input_data[i][j][k] < 0 ? 0 : input_data[i][j][k];
        output_data[i][j][k] =
            layer.alpha[i] * (output_data[i][j][k] - layer.mu[i]) +
            layer.beta[i];
      }
    }
  }

#ifdef _DEBUG
  std::vector<std::vector<std::vector<double>>> debug_data(output_data.size());
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

inline void InferMaxPooling(Image const& input_image,
                            Image& output_image) {
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  for (size_t i = 0; i < output_data.size(); ++i) {
    for (size_t j = 0; j < output_data[0].size(); ++j) {
      for (size_t k = 0; k < output_data[0][0].size(); ++k) {
        std::array<double, 4> rect_data;
        rect_data[0] = input_data[i][j * 2][k * 2];
        rect_data[1] = input_data[i][j * 2][k * 2 + 1];
        rect_data[2] = input_data[i][j * 2 + 1][k * 2];
        rect_data[3] = input_data[i][j * 2 + 1][k * 2 + 1];
        output_data[i][j][k] =
            *std::max_element(rect_data.begin(), rect_data.end());
      }
    }
  }
#ifdef _DEBUG
  std::vector<std::vector<std::vector<double>>> debug_data(output_data.size());
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

inline void InferDense(Para::DenseLayer const& layer,
                       Image const& input_image, Image& output_image) {
  Tick tick(__FN__);
  auto const& input_data = input_image.pixels;
  auto& output_data = output_image.pixels;
  assert(input_data[0].size() == 1 && input_data[0][0].size() == 1);
  assert(output_data[0].size() == 1 && output_data[0][0].size() == 1);
  std::vector<double> input(input_data.size() + 1);
  for (size_t i = 0; i < input_data.size(); ++i) {
    input[i] = input_data[i][0][0];
  }
  input.back() = 1.0;

  for (size_t i = 0; i < output_data.size(); ++i) {    
    //for (size_t j = 0; j < input.size(); ++j) {
    //  std::cout << input[j] << ";";
    //  std::cout << layer.weight[i][j] << "\n";
    //}
    output_data[i][0][0] = std::inner_product(input.begin(), input.end(),
                                              layer.weight[i].begin(), 0.0);
    //std::cout << output_data[i][0][0] << "\n";
  }
}

inline void Test() {
  std::string path = "e:/code/crypto/pod_doc/vgg16_2";
  Para para(path + "/features");

  std::array<std::unique_ptr<Image>, 35> images;
  for (size_t i = 0; i < images.size(); ++i) {
    images[i].reset(new Image(kImageInfos[i]));
  }

  LoadTestImage(path + "/test_image", *images[0]);

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

  // -2.879224,-12.153086,6.915432,10.567657,4.284396,5.724507,-0.843318,
  // -6.088152,2.954222,-9.378846
  images[34]->dump();

  uint32_t max_value = 0;
  for (auto const& image : images) {
    auto const& pixels = image->pixels;
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          if (uint32_t(std::abs(pixels[i][j][k])) > max_value) {
            max_value = uint32_t(std::abs(pixels[i][j][k]));
          }
        }
      }
    }
  }
  std::cout << "max_value: " << max_value; // 176
}
}  // namespace dbl
}  // namespace vgg16
}  // namespace clink