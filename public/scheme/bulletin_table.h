#pragma once

#include <stdint.h>

#include <string>
#include <vector>

#include "misc/misc.h"
#include "ecc/ecc.h"
#include "utils/mkl_tree.h"
#include "public.h"

namespace scheme {
namespace table {

struct Bulletin {
  uint64_t n;
  uint64_t s;
  h256_t sigma_mkl_root;
  h256_t vrf_meta_digest;
};

inline bool IsBulletinValid(Bulletin const& bulletin) {
  auto const& ecc_pub = GetEccPub();
  bool ret = bulletin.n && bulletin.s > 0 && bulletin.s <= ecc_pub.u1().size();
  if (!ret) {
    std::cerr << "invalid bulletin: n=" << bulletin.n << ", s=" << bulletin.s
              << "\n";
  }
  return ret;
}

inline bool SaveBulletin(std::string const& output, Bulletin const& bulletin) {
  if (!IsBulletinValid(bulletin)) {
    assert(false);
    return false;
  }

  try {
    pt::ptree tree;
    tree.put("mode", "table");
    tree.put("n", bulletin.n);
    tree.put("s", bulletin.s);
    tree.put("sigma_mkl_root", misc::HexToStr(bulletin.sigma_mkl_root));
    tree.put("vrf_meta_digest", misc::HexToStr(bulletin.vrf_meta_digest));
    pt::write_json(output, tree);
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool LoadBulletin(std::string const& input, Bulletin& bulletin) {
  try {
    pt::ptree tree;
    pt::read_json(input, tree);
    if (tree.get<std::string>("mode") != "table") throw std::exception();
    bulletin.n = tree.get<uint64_t>("n");
    bulletin.s = tree.get<uint64_t>("s");
    misc::HexStrToH256(tree.get<std::string>("sigma_mkl_root"),
                       bulletin.sigma_mkl_root);
    misc::HexStrToH256(tree.get<std::string>("vrf_meta_digest"),
                       bulletin.vrf_meta_digest);
    if (!IsBulletinValid(bulletin)) throw std::exception();
    return true;
  } catch (std::exception&) {
    std::cout << "bulletin invalid: " << input << "\n";
    assert(false);
    return false;
  }
}

}  // namespace table
}  // namespace scheme