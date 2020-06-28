#pragma once

#include "./para_dbl.h"
#include "ecc/ecc.h"
#include "misc/misc.h"
#include "utils/fst.h"

namespace clink {
namespace vgg16 {

struct Para {
  struct ConvLayer {
    ConvLayer() {}
    ConvLayer(ConvLayerInfo const& info)
        : order(info.order), dimension(info.dimension) {
      coefs.resize(info.kernel_count);
      for (auto& i : coefs) {
        i.resize(info.channel_count);
      }
      bias.resize(info.kernel_count);
    }

    ConvLayer(dbl::Para::ConvLayer const& dbl)
        : order(dbl.order), dimension(dbl.dimension) {
      namespace fp = circuit::fp;
      coefs.resize(dbl.coefs.size());
      for (auto& i : coefs) {
        i.resize(dbl.coefs[0].size());
      }
      bias.resize(dbl.bias.size());

      for (size_t i = 0; i < coefs.size(); ++i) {
        for (size_t j = 0; j < coefs[0].size(); ++j) {
          for (size_t k = 0; k < 3; ++k) {
            for (size_t l = 0; l < 3; ++l) {
              coefs[i][j][k][l] =
                  fp::DoubleToRational<8, 24>(dbl.coefs[i][j][k][l]);
            }
          }
        }
      }
      for (size_t i = 0; i < bias.size(); ++i) {
        bias[i] = fp::DoubleToRational<8, 24>(dbl.bias[i]);
      }
    }

    bool operator==(ConvLayer const& b) const {
      return order == b.order && dimension == b.dimension && coefs == b.coefs &&
             bias == b.bias;
    }

    bool operator!=(ConvLayer const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("vgg16.para.conv", ("o", order), ("d", dimension),
                         ("c", coefs), ("b", bias));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("vgg16.para.conv", ("o", order), ("d", dimension),
                         ("c", coefs), ("b", bias));
    }

    size_t order;
    size_t dimension;
    // size=K*C*3*3
    std::vector<std::vector<std::array<std::array<Fr, 3>, 3>>> coefs;
    std::vector<Fr> bias;  // size=K
    size_t kernel_count() const { return bias.size(); }
    size_t channel_count() const { return coefs[0].size(); }
  };

  struct BnLayer {
    BnLayer() {}
    BnLayer(BnLayerInfo const& info) : order(info.order) {
      mu.resize(info.channel_count);
      alpha.resize(info.channel_count);
      beta.resize(info.channel_count);
    }
    BnLayer(dbl::Para::BnLayer const& dbl) : order(dbl.order) {
      namespace fp = circuit::fp;
      mu.resize(dbl.mu.size());
      alpha.resize(dbl.alpha.size());
      beta.resize(dbl.beta.size());
      for (size_t i = 0; i < mu.size(); ++i) {
        mu[i] = fp::DoubleToRational<8, 24>(dbl.mu[i]);
        alpha[i] = fp::DoubleToRational<8, 24>(dbl.alpha[i]);
        beta[i] = fp::DoubleToRational<8, 24>(dbl.beta[i]);
      }
    }

    bool operator==(BnLayer const& b) const {
      return order == b.order && mu == b.mu && alpha == b.alpha &&
             beta == b.beta;
    }

    bool operator!=(BnLayer const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("vgg16.para.bn", ("o", order), ("m", mu), ("a", alpha),
                         ("b", beta));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("vgg16.para.bn", ("o", order), ("m", mu), ("a", alpha),
                         ("b", beta));
    }

    size_t order;
    std::vector<Fr> mu;
    std::vector<Fr> alpha;
    std::vector<Fr> beta;
  };

  struct DenseLayer {
    DenseLayer() {}
    DenseLayer(DenseLayerInfo const& info) : order(info.order) {
      weight.resize(info.type_count);
      for (auto& i : weight) {
        i.resize(info.input_count + 1);
      }
    }
    DenseLayer(dbl::Para::DenseLayer const& dbl) : order(dbl.order) {
      namespace fp = circuit::fp;
      weight.resize(dbl.weight.size());
      for (auto& i : weight) {
        i.resize(dbl.weight[0].size());
      }
      for (size_t i = 0; i < weight.size(); ++i) {
        for (size_t j = 0; j < weight[0].size(); ++j) {
          weight[i][j] = fp::DoubleToRational<8, 24>(dbl.weight[i][j]);
        }
      }
    }
    bool operator==(DenseLayer const& b) const {
      return order == b.order && weight == b.weight;
    }

    bool operator!=(DenseLayer const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("vgg16.para.dense", ("o", order), ("w", weight));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("vgg16.para.dense", ("o", order), ("w", weight));
    }
    size_t order;
    std::vector<std::vector<Fr>> weight;
  };

  bool operator==(Para const& b) const {
    for (size_t i = 0; i < conv_layers_.size(); ++i) {
      if (conv_layer(i) != b.conv_layer(i)) return false;
    }
    for (size_t i = 0; i < bn_layers_.size(); ++i) {
      if (bn_layer(i) != b.bn_layer(i)) return false;
    }
    for (size_t i = 0; i < dense_layers_.size(); ++i) {
      if (dense_layer(i) != b.dense_layer(i)) return false;
    }

    return true;
  }

  bool operator!=(Para const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.para", ("c", conv_layers_), ("b", bn_layers_),
                       ("d", dense_layers_));
  }

  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.para", ("c", conv_layers_), ("b", bn_layers_),
                       ("d", dense_layers_));
  }

  BnLayer const& bn_layer(size_t order) const { return bn_layers_[order]; }

  ConvLayer const& conv_layer(size_t order) const {
    return conv_layers_[order];
  }

  DenseLayer const& dense_layer(size_t order) const {
    return dense_layers_[order];
  }

  std::array<ConvLayer, 13> const& conv_layers() const { return conv_layers_; }

  std::array<BnLayer, 14> const& bn_layers() const { return bn_layers_; }

  std::array<DenseLayer, 2> const& dense_layers() const {
    return dense_layers_;
  }

  Para() {}

  Para(dbl::Para const& dbl_para) {
    Tick tick(__FN__);
    for (size_t i = 0; i < conv_layers_.size(); ++i) {
      conv_layers_[i] = ConvLayer(dbl_para.conv_layer(i));
    }
    for (size_t i = 0; i < bn_layers_.size(); ++i) {
      bn_layers_[i] = BnLayer(dbl_para.bn_layer(i));
    }
    for (size_t i = 0; i < dense_layers_.size(); ++i) {
      dense_layers_[i] = DenseLayer(dbl_para.dense_layer(i));
    }
  }

  Para(std::string const& file) {
    Tick tick(__FN__);
    if (!Load(file)) {
      throw std::invalid_argument("invalid para file: " + file);
    }
  }

  bool Load(std::string const& file) {
    try {
      yas::file_istream is(file.c_str());
      yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
      ia.serialize(*this);
      return true;
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
  }

  bool Save(std::string const& file) const {
    Tick tick(__FN__);
    try {
      boost::system::error_code dummy;
      fs::remove(file,dummy);
      yas::file_ostream os(file.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      oa.serialize(*this);
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }

#ifdef _DEBUG_CHECK
    try {
      Para check(file);
      if (check != *this) {
        std::cout << "oops\n";
        return false;
      }
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }
#endif

    return true;
  }

 private:
  std::array<ConvLayer, 13> conv_layers_;
  std::array<BnLayer, 14> bn_layers_;
  std::array<DenseLayer, 2> dense_layers_;
};

struct Image {
  size_t const order;
  boost::multi_array<Fr, 3> pixels; // TODO: change to normal vector  
  Image(ImageInfo const& info)
      : order(info.order),
        pixels(boost::extents[info.channel_count][info.dimension]
                             [info.dimension]) {
  }

  Image(dbl::Image const& dbl)
      : order(dbl.order),
        pixels(boost::extents[dbl.channel_count()][dbl.dimension()]
                             [dbl.dimension()]) {
    namespace fp = circuit::fp;
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          pixels[i][j][k] = fp::DoubleToRational<8, 24>(dbl.pixels[i][j][k]);
        }
      }
    }
  }

  std::vector<Fr> copy_vec() const {
    return std::vector<Fr>(pixels.data(),
                           pixels.data() + pixels.num_elements());
  }

  bool operator==(Image const& b) const {
    return order == b.order && pixels == b.pixels;
  }

  bool operator!=(Image const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    std::vector<Fr> data;    
    data.reserve(total_size());
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          data.push_back(pixels[i][j][k]);
        }
      }
    }
    ar& YAS_OBJECT_NVP("vgg16.image", ("p", data));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    std::vector<Fr> data;
    ar& YAS_OBJECT_NVP("vgg16.image", ("p", data));
    if (data.size() != total_size()) {
      throw std::runtime_error("oops");
    }
    auto d = dimension();
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          pixels[i][j][k] = data[i * d * d + j * d + k];
        }
      }
    }
  }

  size_t channel_count() const { return pixels.size(); }
  size_t dimension() const { return pixels[0].size(); }
  size_t total_size() const {
    return pixels.size() * pixels[0].size() * pixels[0][0].size();
  }
  Fr const& get(size_t offset) const {
    size_t D = dimension();
    size_t DD = D * D;
    size_t i = offset / DD;
    size_t ii = offset % DD;
    size_t j = ii / D;
    size_t k = ii % D;
    return pixels[i][j][k];
  }

  template<size_t D, size_t N>
  void dump() const {
    namespace fp = circuit::fp;   
    for (size_t i = 0; i < pixels.size(); ++i) {
      for (size_t j = 0; j < pixels[0].size(); ++j) {
        for (size_t k = 0; k < pixels[0][0].size(); ++k) {
          Fr const& data = pixels[i][j][k];
          double dbl_data = fp::RationalToDouble<D, N>(data);
          std::cout << dbl_data << ";";
        }
      }
    }
    std::cout << "\n\n";
  }


};
}  // namespace vgg16
}  // namespace clink