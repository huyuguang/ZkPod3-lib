#pragma once

#include "./context.h"
#include "./image_com.h"
#include "./safevec.h"
#include "hyrax/hyrax.h"

namespace clink::vgg16 {
struct AdaptProveItem {
  std::string order_tag;
  std::vector<std::vector<Fr>> x;
  std::vector<std::vector<Fr>> a;
  std::vector<G1> cx;
  std::vector<Fr> rx;
  void Init(size_t count, std::string const& tag) {
    x.resize(count);
    a.resize(count);
    cx.resize(count);
    rx.resize(count);
    order_tag = tag;
  }
  bool CheckFormat() const {
    if (x.size() != a.size()) return false;
    if (cx.size() != a.size()) return false;
    if (rx.size() != a.size()) return false;
    for (size_t i = 0; i < x.size(); ++i) {
      if (x[i].size() != a[i].size()) return false;
    }
    return true;
  }
  bool CheckData() const {
    Fr sum = FrZero();
    for (size_t i = 0; i < x.size(); ++i) {
      sum += InnerProduct(x[i], a[i]);
    }
    if (sum != FrZero()) return false;
    bool all_success = false;
    auto parallel_f = [this](int64_t i) {
      return cx[i] == pc::ComputeCom(x[i], rx[i]);
    };
    parallel::For(&all_success, x.size(), parallel_f);
    return all_success;
  }
};

struct AdaptVerifyItem {
  std::string order_tag;
  std::vector<std::vector<Fr>> a;
  std::vector<G1> cx;
  bool CheckFormat() const { return a.size() == cx.size(); }
  void Init(size_t count, std::string const& tag) {
    order_tag = tag;
    a.resize(count);
    cx.resize(count);
  }
};

using AdaptProveItemMan = SafeVec<AdaptProveItem>;
using AdaptVerifyItemMan = SafeVec<AdaptVerifyItem>;
using ParallelVoidTaskMan = SafeVec<parallel::VoidTask>;
using ParallelBoolTaskMan = SafeVec<parallel::BoolTask>;

inline h256_t AdaptItemDigest(std::vector<G1> const& cx,
                              std::vector<std::vector<Fr>> const& a) {
  h256_t digest;
  CryptoPP::Keccak_256 hash;
  for (auto const& i : cx) {
    HashUpdate(hash, i);
  }
  for (auto const& i : a) {
    for (auto const& j : i) {
      HashUpdate(hash, j);
    }
  }
  hash.Final(digest.data());
  return digest;
}

template<typename Item>
void AdaptUpdateSeed(h256_t& seed, std::vector<Item> const& items) {
  std::vector<h256_t> digests(items.size());
  auto parallel_f = [&digests,&items](int64_t i) {
    digests[i] = AdaptItemDigest(items[i].cx, items[i].a);
  };
  parallel::For(items.size(), parallel_f);

  CryptoPP::Keccak_256 hash;
  HashUpdate(hash, seed);
  for (auto const& i : digests) {
    HashUpdate(hash, i);
  }
  hash.Final(seed.data());
}

inline void AdaptComputeFst(h256_t const& seed, std::vector<Fr>& e) {
  std::string salt = "vgg16 adapt " + std::to_string(e.size());
  ComputeFst(seed, salt, e);
}


inline std::vector<size_t> AdaptGetCombinedItemOrder(
    std::vector<std::vector<Fr>> const& combined_a) {
  std::vector<size_t> order(combined_a.size());
  for (size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }

  std::sort(order.begin(), order.end(), [&combined_a](size_t a, size_t b) {
    auto len_a = combined_a[a].size();
    auto len_b = combined_a[b].size();
    if (len_a == len_b) return a > b;
    return len_a > len_b;
  });

  return order;
}

template<typename T>
void AdaptPermuteCombinedItems(std::vector<size_t> const& order,
                               std::vector<T>& combined_v) {
  if (order.size() != combined_v.size()) {
    std::cerr << __FN__ << " oops";
    throw std::runtime_error("oops");
  }

  std::vector<T> v2(combined_v.size());
  for (size_t i = 0; i < order.size();++i) {
    v2[i] = std::move(combined_v[order[i]]);
  }
  combined_v.swap(v2);
}

inline void AdaptProve(h256_t seed, AdaptProveItemMan& item_man,
                       hyrax::A4::Proof& proof) {
  Tick tick(__FN__);
  std::vector<AdaptProveItem> items;
  item_man.take(items);
  if (items.empty()) return;

  std::sort(items.begin(), items.end(),
            [](AdaptProveItem const& a, AdaptProveItem const& b) {
              return a.order_tag < b.order_tag;
            });

#ifdef _DEBUG_CHECK
  for (auto const& i : items) {
    //std::cout << i.order_tag << "\n";

    if (!i.CheckFormat()) {
      std::string errmsg = i.order_tag + " oops";
      std::cout << errmsg << "\n";
      throw std::runtime_error(errmsg);
    }
    if (!i.CheckData()) {
      std::string errmsg = i.order_tag + " oops";
      std::cout << errmsg << "\n";
      throw std::runtime_error(errmsg);
    }
  }
#endif    
  
  for (auto const& i : items) {
    std::cout << __FN__ << " " << i.order_tag << "," << i.a.size() << "*"
              << i.a[0].size() << "\n";
  }

  AdaptUpdateSeed(seed, items);

  std::vector<Fr> e(items.size());
  AdaptComputeFst(seed, e);
  std::cout << __FN__ << "  " << e[0] << "\n";

  auto pf = [&items,&e](int64_t i) {
    for (auto& j : items[i].a) {
      j *= e[i];
    }
  };
  parallel::For(items.size(), pf);

  // combine
  size_t vector_count = 0;
  for (auto const& i : items) {
    vector_count += i.x.size();
  }

  std::vector<std::vector<Fr>> combined_x(vector_count);
  std::vector<std::vector<Fr>> combined_a(vector_count);
  std::vector<G1> combined_cx(vector_count);
  std::vector<Fr> combined_rx(vector_count);

  size_t cursor = 0;
  for (auto& item : items) {
    for (size_t i = 0; i < item.x.size(); ++i) {
      combined_x[cursor] = std::move(item.x[i]);
      combined_a[cursor] = std::move(item.a[i]);
      combined_cx[cursor] = std::move(item.cx[i]);
      combined_rx[cursor] = std::move(item.rx[i]);
      ++cursor;
    }
  }

  // TODO: move to hyrax::A4
  auto order = AdaptGetCombinedItemOrder(combined_a);
  AdaptPermuteCombinedItems(order, combined_a);
  AdaptPermuteCombinedItems(order, combined_x);
  AdaptPermuteCombinedItems(order, combined_cx);
  AdaptPermuteCombinedItems(order, combined_rx);

  hyrax::A4::ProveInput input(std::move(combined_x), std::move(combined_a),
                              FrZero(), pc::kGetRefG1, pc::PcU());
  hyrax::A4::CommitmentPub com_pub;
  com_pub.cx = std::move(combined_cx);
  com_pub.cz = G1Zero();

  hyrax::A4::CommitmentSec com_sec;
  com_sec.r = std::move(combined_rx);
  com_sec.t = FrZero();

  hyrax::A4::AlignData(input, com_pub, com_sec);

  hyrax::A4::Prove(proof, seed, std::move(input), std::move(com_pub), 
                   std::move(com_sec));
}

inline bool AdaptVerify(h256_t seed, AdaptVerifyItemMan& item_man,
                        hyrax::A4::Proof const& proof) {
  Tick tick(__FN__);

  std::vector<AdaptVerifyItem> items;
  item_man.take(items);
  if (items.empty()) return true;

  std::sort(items.begin(), items.end(),
            [](AdaptVerifyItem const& a, AdaptVerifyItem const& b) {
              return a.order_tag < b.order_tag;
            });

  for (auto const& i : items) {    
    if (!i.CheckFormat()) {
      std::cout << i.order_tag << " format error\n";
      throw std::runtime_error("oops");
    }
  }

  AdaptUpdateSeed(seed, items);

  std::vector<Fr> e(items.size());
  AdaptComputeFst(seed, e);
  std::cout << __FN__ << " " << e[0] << "\n";

  auto pf = [&items,&e](int64_t i) {
    for (auto& j : items[i].a) {
      j *= e[i];
    }
  };
  parallel::For(items.size(), pf);

  // combine
  size_t vector_count = 0;
  for (auto const& i : items) {
    vector_count += i.a.size();
  }

  std::vector<std::vector<Fr>> combined_a(vector_count);
  std::vector<G1> combined_cx(vector_count);

  size_t cursor = 0;
  for (auto& item : items) {
    for (size_t i = 0; i < item.a.size(); ++i) {
      combined_a[cursor] = std::move(item.a[i]);
      combined_cx[cursor] = std::move(item.cx[i]);
      ++cursor;
    }
  }

  auto order = AdaptGetCombinedItemOrder(combined_a);
  AdaptPermuteCombinedItems(order, combined_a);
  AdaptPermuteCombinedItems(order, combined_cx);
  
  hyrax::A4::CommitmentPub com_pub;
  com_pub.cx = std::move(combined_cx);
  com_pub.cz = G1Zero();
  com_pub.Align();
  hyrax::A4::VerifyInput input(com_pub, pc::kGetRefG1, std::move(combined_a),
                               pc::PcU());
  bool success = hyrax::A4::Verify(proof, seed, input);
  std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}

}  // namespace clink::vgg16
