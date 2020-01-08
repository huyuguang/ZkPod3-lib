#pragma once

#include <stdint.h>

#include <string>

#include "misc/misc.h"
#include "scheme/bulletin_plain.h"
#include "utils/mkl_tree.h"
#include "scheme/public_misc.h"

namespace scheme::plain {

class BobData {
 public:
  BobData(Bulletin const& bulletin, std::string const& public_path);
  BobData(std::string const& bulletin_file, std::string const& public_path);
  Bulletin const& bulletin() const { return bulletin_; }
  std::vector<G1> sigmas() const { return sigmas_; }
  bool SaveDecryped(std::string const& file, std::vector<Range> const& demands,
                    std::vector<Fr> const& decrypted);

 private:
  void LoadData();
  bool NeedVerify();

 private:
  Bulletin bulletin_;
  std::string public_path_;

 private:
  std::vector<G1> sigmas_;
};

typedef std::shared_ptr<BobData> BobDataPtr;

}  // namespace scheme::plain