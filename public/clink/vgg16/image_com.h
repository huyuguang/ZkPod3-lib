#pragma once

#include "./infer.h"

namespace clink::vgg16 {

struct ImageCommitmentPub {
  std::array<G1, 35> c;

  ImageCommitmentPub() {}

  ImageCommitmentPub(std::string const& file) {
    if (!YasLoadBin(file, *this)) {
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
};

struct ImageCommitmentSec {
  std::array<Fr, 35> r;

  ImageCommitmentSec() {}

  ImageCommitmentSec(std::string const& file) {
    if (!YasLoadBin(file, *this)) {
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
};

inline void ComputePerImageCommitment(
    std::array<std::unique_ptr<Image>, 35> const& images,
    ImageCommitmentPub& pub, ImageCommitmentSec& sec) {
  sec.r[0] = FrZero();  // the image[0] is public input
  FrRand(sec.r.data() + 1, sec.r.size() - 1);

  auto parallel_f = [&images, &pub, &sec](int64_t i) {
    auto const& image = *images[i];
    auto const& r = sec.r[i];
    pub.c[i] = pc::ComputeCom(image.data, r);
  };

  parallel::For(images.size(), parallel_f);
}

}  // namespace clink::vgg16