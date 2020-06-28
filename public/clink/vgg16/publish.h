#pragma once

#include "./para_fr.h"
#include "./para_com.h"

namespace clink::vgg16 {

// --vgg16_publish "E:/code/crypto/pod_doc/vgg16_2 e:/vgg16_publish"

inline bool Publish(std::string const& para_path,
                    std::string const& working_path) {
  Tick tick(__FN__);
  if (para_path.empty() || !fs::is_directory(para_path)) {
    std::cerr << "Open publish_dir " << para_path << " failed\n";
    return false;
  }

  if (!fs::is_directory(working_path) &&
      !fs::create_directories(working_path)) {
    std::cerr << "Create " << working_path << " failed\n";
    return false;
  }

  std::string pub_path = working_path + "/pub";
  std::string sec_path = working_path + "/sec";
  if (!fs::is_directory(pub_path) &&
      !fs::create_directories(pub_path)) {
    std::cerr << "Create " << pub_path << " failed\n";
    return false;
  }
  if (!fs::is_directory(sec_path) &&
      !fs::create_directories(sec_path)) {
    std::cerr << "Create " << sec_path << " failed\n";
    return false;
  }

  try {
    dbl::Para dbl_para(para_path);
    std::unique_ptr<Para> para(new Para(dbl_para));
    if (!para->Save(sec_path + "/para")) {
      std::cerr << "save para failed\n";
      return false;
    }

    std::unique_ptr<AuxiPub> auxi(new AuxiPub);
    if (!auxi->Save(pub_path + "/auxi")) {
      std::cerr << "save auxi failed\n";
      return false;
    }

    ParaCommitmentPub com_pub;
    ParaCommitmentSec com_sec;
    ComputeParaCommitment(*para, *auxi, com_pub, com_sec);
    if (!com_pub.Save(pub_path + "/para_com_pub")) {
      std::cerr << "save com_pub failed\n";
      return false;
    }
    if (!com_sec.Save(sec_path + "/para_com_sec")) {
      std::cerr << "save com_sec failed\n";
      return false;
    }
    return true;
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return false;
  }
}

}  // namespace clink::vgg16