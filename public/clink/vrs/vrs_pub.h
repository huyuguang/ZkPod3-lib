#pragma once

#include "../details.h"
#include "./vrs_scheme.h"

namespace clink {

template <typename Scheme>
struct VrsPub {
  struct Item {
    int64_t begin;
    int64_t end;
    int64_t n() const { return end - begin; }
  };

  static std::vector<Item> SplitLargeTask(int64_t n) {
    auto const kMaxUnitPerZkp = Scheme::kMaxUnitPerZkp;
    std::vector<Item> items((n + kMaxUnitPerZkp - 1) / kMaxUnitPerZkp);
    for (int64_t i = 0; i < (int64_t)items.size(); ++i) {
      auto& item = items[i];
      item.begin = i * kMaxUnitPerZkp;
      item.end = item.begin + kMaxUnitPerZkp;
      if (item.end > n) {
        item.end = n;
      }
    }
    return items;
  }

  static void GeneratePlain(Fr* out, h256_t const& plain_seed,
                            int64_t position) {
    CryptoPP::Keccak_256 hash;
    h256_t plain;
    hash.Update(plain_seed.data(), plain_seed.size());
    auto offset_big = boost::endian::native_to_big(position);
    hash.Update((uint8_t const*)&offset_big, sizeof(offset_big));
    hash.Final(plain.data());
    *out = H256ToFr(plain);
  }

  static Fr GeneratePlain(h256_t const& plain_seed, int64_t position) {
    Fr ret;
    GeneratePlain(&ret, plain_seed, position);
    return ret;
  }

  static bool VerifySecret(G1 const& h, G1 const& g, G1 const& key_com,
                           Fr const& r, Fr const& key) {
    return key_com == h * r + g * key;
  }
};
}  // namespace clink