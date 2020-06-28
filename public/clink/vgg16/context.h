#pragma once

#include "./image_com.h"
#include "./para_com.h"

namespace clink::vgg16 {
class ProveContext {
 public:
  ProveContext(std::string const& working_path) {
    auto pub_path = working_path + "/pub";
    auto sec_path = working_path + "/sec";
    para_.reset(new Para(sec_path + "/para"));
    for (size_t i = 0; i < images_.size(); ++i) {
      images_[i].reset(new Image(kImageInfos[i]));
      std::string file = sec_path + "/image_" + std::to_string(i);
      yas::file_istream is(file.c_str());
      yas::binary_iarchive<yas::file_istream, YasBinF()> ia(is);
      ia.serialize(*images_[i]);
      const_images_[i] = images_[i].get();
    }

    auxi_.reset(new AuxiPub(pub_path + "/auxi"));
    para_com_pub_.reset(new ParaCommitmentPub(pub_path + "/para_com_pub"));
    para_com_sec_.reset(new ParaCommitmentSec(sec_path + "/para_com_sec"));
    image_com_pub_.reset(new ImageCommitmentPub(pub_path + "/image_com_pub"));
    image_com_sec_.reset(new ImageCommitmentSec(sec_path + "/image_com_sec"));

#ifdef _DEBUG_CHECK

#endif
  }

  std::array<Image const*, 35> const& const_images() const {
    return const_images_;
  }

  AuxiPub const& auxi() const { return *auxi_; }

  Para const& para() const { return *para_; }

  ParaCommitmentPub const& para_com_pub() const { return *para_com_pub_; }

  ParaCommitmentSec const& para_com_sec() const { return *para_com_sec_; }

  ImageCommitmentPub const& image_com_pub() const { return *image_com_pub_; }

  ImageCommitmentSec const& image_com_sec() const { return *image_com_sec_; }

 private:
  std::array<std::unique_ptr<Image>, 35> images_;
  std::array<Image const*, 35> const_images_;
  std::unique_ptr<AuxiPub> auxi_;
  std::unique_ptr<Para> para_;
  std::unique_ptr<ParaCommitmentPub> para_com_pub_;
  std::unique_ptr<ParaCommitmentSec> para_com_sec_;
  std::unique_ptr<ImageCommitmentPub> image_com_pub_;
  std::unique_ptr<ImageCommitmentSec> image_com_sec_;
};

class VerifyContext {
 public:
  VerifyContext(std::string const& pub_path) {
    auxi_.reset(new AuxiPub(pub_path + "/auxi"));
    para_com_pub_.reset(new ParaCommitmentPub(pub_path + "/para_com_pub"));
    image_com_pub_.reset(new ImageCommitmentPub(pub_path + "/image_com_pub"));
  }

  AuxiPub const& auxi() const { return *auxi_; }

  ParaCommitmentPub const& para_com_pub() const { return *para_com_pub_; }

  ImageCommitmentPub const& image_com_pub() const { return *image_com_pub_; }

 private:
  std::unique_ptr<AuxiPub> auxi_;
  std::unique_ptr<ParaCommitmentPub> para_com_pub_;
  std::unique_ptr<ImageCommitmentPub> image_com_pub_;
};
}  // namespace clink::vgg16