#pragma once

#include "./infer.h"

namespace clink::vgg16 {

struct ImageCommitmentPub {
  std::array<G1, 35> c;

  ImageCommitmentPub() {}
  ImageCommitmentPub(std::string const& file) {
    if (!Load(file)) {
      throw std::invalid_argument("invalid image commitment pub file: " + file);
    }
  }

  bool operator==(ImageCommitmentPub const& b) const { return c == b.c; }

  bool operator!=(ImageCommitmentPub const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.image.compub", ("c", c));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.image.compub", ("c", c));
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
    try {
      boost::system::error_code dummy;
      fs::remove(file);
      yas::file_ostream os(file.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      oa.serialize(*this);
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }

#ifdef _DEBUG_CHECK
    ImageCommitmentPub check;
    if (!check.Load(file)) {
      std::cout << "oops\n";
      return false;
    }
    if (check != *this) {
      std::cout << "oops\n";
      return false;
    }
#endif

    return true;
  }
};

struct ImageCommitmentSec {
  std::array<Fr, 35> r;

  ImageCommitmentSec() {}
  ImageCommitmentSec(std::string const& file) {
    if (!Load(file)) {
      throw std::invalid_argument("invalid image commitment sec file: " + file);
    }
  }

  bool operator==(ImageCommitmentSec const& b) const { return r == b.r; }

  bool operator!=(ImageCommitmentSec const& b) const { return !(*this == b); }

  template <typename Ar>
  void serialize(Ar& ar) const {
    ar& YAS_OBJECT_NVP("vgg16.image.comsec", ("r", r));
  }
  template <typename Ar>
  void serialize(Ar& ar) {
    ar& YAS_OBJECT_NVP("vgg16.image.comsec", ("r", r));
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
    try {
      boost::system::error_code dummy;
      fs::remove(file);
      yas::file_ostream os(file.c_str());
      yas::binary_oarchive<yas::file_ostream, YasBinF()> oa(os);
      oa.serialize(*this);
    } catch (std::exception& e) {
      std::cerr << e.what() << "\n";
      return false;
    }

#ifdef _DEBUG_CHECK
    ImageCommitmentSec check;
    if (!check.Load(file)) {
      std::cout << "oops\n";
      return false;
    }
    if (check != *this) {
      std::cout << "oops\n";
      return false;
    }
#endif

    return true;
  }
};

inline void ComputePerImageCommitment(
    std::array<std::unique_ptr<Image>, 35> const& images,
    ImageCommitmentPub& pub, ImageCommitmentSec& sec) {
  sec.r[0] = FrZero();  // the image[0] is public input
  FrRand(sec.r.data() + 1, sec.r.size() - 1);

  auto parallel_f = [&images, &pub, &sec](int64_t i) {
    auto const& image = *images[i];
    auto const& r = sec.r[i];
    auto get_x = [&image](int64_t j) -> Fr const& { return image.get(j); };
    pub.c[i] = pc::PcComputeCommitmentG(image.total_size(), get_x, r);
  };

  parallel::For(images.size(), parallel_f);
}

}  // namespace clink::vgg16