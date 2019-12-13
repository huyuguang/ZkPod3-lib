#pragma once

#include <stdint.h>

#include <string>
#include <vector>

#include "basic_types.h"
#include "bp.h"
#include "chain.h"
#include "ecc.h"
#include "ecc_pub.h"
#include "misc.h"
#include "mkl_tree.h"
#include "mpz.h"
#include "multiexp.h"
#include "public.h"
#include "vrf.h"

namespace scheme {

enum Mode { kPlain, kTable };

enum Action {
  kVrfQuery,
  kOtVrfQuery,
  kVrfPod,
  kOtVrfPod,
  kComplaintPod,
  kOtComplaintPod,
  kAtomicSwapPod,
  kAtomicSwapPodVc,
};

namespace details {

inline h256_t KToH256(G1 const& g) {
  assert(g.isNormalized());

  uint8_t x[kFpBinSize] = {0};
  if (g.z == 1) g.x.serialize(x, kFpBinSize);

  uint8_t y[kFpBinSize] = {0};
  if (g.z == 1) g.y.serialize(y, kFpBinSize);

  h256_t digest;
  CryptoPP::Keccak_256 hash;
  hash.Update(x, sizeof(x));
  hash.Update(y, sizeof(y));
  hash.Final(digest.data());
  return digest;
}

}  // namespace details

inline void LoadMij(uint8_t const* data_start, uint8_t const* data_end,
                    uint64_t i, uint64_t j, uint64_t s, Fr& mij) {
  auto offset = i * s + j;
  uint8_t const* p = data_start + offset * 31;
  uint8_t const* q = p + 31;
  if (p >= data_end) {
    mij = FrZero();
  } else {
    if (q > data_end) q = data_end;
    mij = BinToFr31(p, q);
  }
}

inline bool CopyData(std::string const& src, std::string const& dst) {
  try {
    io::mapped_file_params src_params;
    src_params.path = src;
    src_params.flags = io::mapped_file_base::readonly;
    io::mapped_file_source src_view(src_params);

    io::mapped_file_params dst_params;
    dst_params.path = dst;
    dst_params.flags = io::mapped_file_base::readwrite;
    dst_params.new_file_size = src_view.size();
    io::mapped_file dst_view(dst_params);

    memcpy(dst_view.data(), src_view.data(), src_view.size());

    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool SaveMkl(std::string const& output,
                    std::vector<h256_t> const& mkl_tree) {
  Tick _tick_(__FUNCTION__);
  constexpr size_t kItemSize = 32;  // h256_t
  try {
    io::mapped_file_params params;
    params.path = output;
    params.flags = io::mapped_file_base::readwrite;
    params.new_file_size = mkl_tree.size() * kItemSize;
    io::mapped_file view(params);
    uint8_t* start = (uint8_t*)view.data();
    for (size_t i = 0; i < mkl_tree.size(); ++i) {
      h256_t const& h = mkl_tree[i];
      memcpy(start + i * kItemSize, h.data(), kItemSize);
    }
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool LoadMkl(std::string const& input, uint64_t n,
                    std::vector<h256_t>& mkl_tree) {
  constexpr size_t kItemSize = 32;  // h256_t
  try {
    io::mapped_file_params params;
    params.path = input;
    params.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view(params);
    auto tree_size = mkl::GetTreeSize(n);
    if (view.size() != tree_size * kItemSize) {
      assert(false);
      return false;
    }
    mkl_tree.resize(tree_size);
    auto start = (uint8_t*)view.data();
    for (size_t i = 0; i < tree_size; ++i) {
      h256_t& h = mkl_tree[i];
      memcpy(h.data(), start + i * h.size(), h.size());
    }

    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool SaveSigma(std::string const& output, std::vector<G1> const& sigma) {
  Tick _tick_(__FUNCTION__);
  try {
    io::mapped_file_params params;
    params.path = output;
    params.flags = io::mapped_file_base::readwrite;
    params.new_file_size = sigma.size() * kG1CompBinSize;
    io::mapped_file view(params);
    uint8_t* start = (uint8_t*)view.data();
    for (size_t i = 0; i < sigma.size(); ++i) {
      G1ToBin(sigma[i], start + i * kG1CompBinSize);
    }
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool LoadSigma(std::string const& input, uint64_t n, h256_t const* root,
                      std::vector<G1>& sigmas) {
  try {
    io::mapped_file_params params;
    params.path = input;
    params.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view(params);
    if (view.size() != n * kG1CompBinSize) return false;
    auto start = (uint8_t*)view.data();

    if (root) {
      auto get_sigma = [start, n](uint64_t i) -> h256_t {
        assert(i < n);
        (void)n;
        h256_t h;
        memcpy(h.data(), start + i * kG1CompBinSize, kG1CompBinSize);
        return h;
      };
      if (*root != mkl::CalcRoot(std::move(get_sigma), n)) {
        assert(false);
        return false;
      }
    }

    sigmas.resize(n);

    auto parallel_f = [start, &sigmas](int64_t i) mutable {
      sigmas[i] = BinToG1(start + i * kG1CompBinSize);
    };
    parallel::For((int64_t)n, parallel_f);
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool SaveMatrix(std::string const& output, std::vector<Fr> const& m) {
  Tick _tick_(__FUNCTION__);
  try {
    io::mapped_file_params params;
    params.path = output;
    params.flags = io::mapped_file_base::readwrite;
    params.new_file_size = m.size() * kFrBinSize;
    io::mapped_file view(params);
    uint8_t* start = (uint8_t*)view.data();
    for (size_t i = 0; i < m.size(); ++i) {
      FrToBin(m[i], start + i * kFrBinSize);
    }
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool LoadMatrix(std::string const& input, uint64_t ns,
                       std::vector<Fr>& m) {
  try {
    io::mapped_file_params params;
    params.path = input;
    params.flags = io::mapped_file_base::readonly;
    io::mapped_file_source view(params);
    if (view.size() != kFrBinSize * ns) return false;

    auto start = (uint8_t*)view.data();
    m.resize(ns);
    for (uint64_t i = 0; i < m.size(); ++i) {
      if (!BinToFr32(start + i * kFrBinSize, &m[i])) {
        assert(false);
        return false;
      }
    }
    return true;
  } catch (std::exception&) {
    return false;
  }
}

inline std::vector<G1> CalcSigma(std::vector<Fr> const& m, uint64_t n,
                                 uint64_t s) {
  assert(m.size() == n * s);

  auto const& ecc_pub = GetEccPub();

  auto const& u1 = ecc_pub.u1();
  std::vector<G1> sigmas(n);
  auto parallel_f = [s, &m, &u1, &sigmas](int64_t i) {
    G1& sigma = sigmas[i];
    Fr const* mi0 = &m[i * s];
    sigma = MultiExpBdlo12(u1.data(), mi0, s);
  };
  parallel::For((int64_t)n, parallel_f);
  return sigmas;
}

inline std::vector<h256_t> BuildSigmaMklTree(std::vector<G1> const& sigmas) {
  Tick _tick_(__FUNCTION__);
  auto get_sigma = [&sigmas](uint64_t i) -> h256_t {
    return G1ToBin(sigmas[i]);
  };
  return mkl::BuildTree(sigmas.size(), get_sigma);
}

inline bool GetBulletinMode(std::string const& file, Mode& mode) {
  try {
    pt::ptree tree;
    pt::read_json(file, tree);
    auto str = tree.get<std::string>("mode");
    if (str == "table") {
      mode = Mode::kTable;
    } else if (str == "plain") {
      mode = Mode::kPlain;
    } else {
      return false;
    }
    return true;
  } catch (std::exception&) {
    assert(false);
    return false;
  }
}

inline bool IsElementUnique(std::vector<Fr> const v) {
  std::vector<Fr const*> pv(v.size());
  for (size_t i = 0; i < v.size(); ++i) pv[i] = &v[i];

  auto compare = [](Fr const* a, Fr const* b) {
    return a->getMpz() < b->getMpz();
  };

  std::sort(pv.begin(), pv.end(), compare);

  return std::adjacent_find(pv.begin(), pv.end(), compare) == pv.end();
}

// since we need to verify the mkl path in contract, we use plain G1
inline h256_t CalcRootOfK(std::vector<G1> const& k) {
  Tick _tick_(__FUNCTION__);
  auto get_k = [&k](uint64_t i) -> h256_t {
    assert(i < k.size());
    return details::KToH256(k[i]);
  };
  return mkl::CalcRoot(std::move(get_k), k.size());
}

// since we need to verify the mkl path in contract, we use plain G1
inline h256_t CalcPathOfK(std::vector<G1> const& k, uint64_t ij,
                          std::vector<h256_t>& path) {
  auto root = mkl::CalcPath(
      [&k](uint64_t i) -> h256_t {
        assert(i < k.size());
        return details::KToH256(k[i]);
      },
      k.size(), ij, &path);
  return root;
}

inline h256_t CalcRangesDigest(std::vector<Range> const& r) {
  h256_t digest;
  CryptoPP::Keccak_256 hash;
  for (auto& i : r) {
    auto a = boost::endian::native_to_big(i.start);
    auto b = boost::endian::native_to_big(i.count);
    hash.Update((uint8_t*)&a, sizeof(a));
    hash.Update((uint8_t*)&b, sizeof(b));
  }
  hash.Final(digest.data());
  return digest;
}

inline h256_t CalcFrDataDigest(std::vector<Fr> const& m) {
  h256_t digest;
  CryptoPP::Keccak_256 hash;
  for (auto& i : m) {
    h256_t bin = FrToBin(i);
    hash.Update(bin.data(), bin.size());
  }
  hash.Final(digest.data());
  return digest;
}

inline h256_t CalcG1DataDigest(std::vector<G1> const& d) {
  h256_t digest;
  CryptoPP::Keccak_256 hash;
  for (auto& i : d) {
    h256_t bin = G1ToBin(i);
    hash.Update(bin.data(), bin.size());
  }
  hash.Final(digest.data());
  return digest;
}

// since we need to verify the mkl path in contract, we use plain G1
inline bool VerifyPathOfK(G1 const& ki, uint64_t i, uint64_t n,
                          h256_t const& root, std::vector<h256_t> const& path) {
  h256_t k_bin = details::KToH256(ki);
  return mkl::VerifyPath(i, k_bin, n, root, path);
}

inline void BuildK(std::vector<Fr> const& v, std::vector<G1>& k, uint64_t s) {
  Tick _tick_(__FUNCTION__);

  assert(v.size() % s == 0);

  int64_t n = (int64_t)(v.size() / s);
  k.resize(n);

  auto parallel_f = [&v, &k, s](int64_t i) {
    Fr const* vi0 = &v[i * s];
    k[i] = MultiExpU1(s, [vi0](uint64_t j) -> Fr const& { return vi0[j]; });
    k[i].normalize();
  };
  parallel::For(n, parallel_f);
}

inline h256_t CalcSeed2(std::vector<h256_t> const& h) {
  h256_t digest;
  CryptoPP::Keccak_256 hash;
  for (auto& i : h) {
    hash.Update(i.data(), i.size());
  }
  hash.Final(digest.data());
  return digest;
}

inline bool CheckDemandPhantoms(uint64_t n, std::vector<Range> const& demands,
                                std::vector<Range> const& phantoms) {
  if (demands.empty() || phantoms.empty()) return false;

  for (auto const& demand : demands) {
    if (!demand.count || demand.start >= n || demand.count > n ||
        (demand.start + demand.count) > n)
      return false;
  }

  for (auto const& phantom : phantoms) {
    if (!phantom.count || phantom.start >= n || phantom.count > n ||
        (phantom.start + phantom.count) > n)
      return false;
  }

  for (size_t i = 1; i < demands.size(); ++i) {
    if (demands[i].start <= demands[i - 1].start + demands[i - 1].count)
      return false;
  }

  for (size_t i = 1; i < phantoms.size(); ++i) {
    if (phantoms[i].start <= phantoms[i - 1].start + phantoms[i - 1].count)
      return false;
  }

  for (size_t i = 0; i < demands.size(); ++i) {
    auto d_start = demands[i].start;
    auto d_end = d_start + demands[i].count;
    bool find = false;
    for (size_t j = 0; j < phantoms.size(); ++j) {
      auto p_start = phantoms[j].start;
      auto p_end = p_start + phantoms[j].count;
      if (p_start <= d_start && p_end >= d_end) {
        find = true;
        break;
      }
    }
    if (!find) return false;
  }
  return true;
}

inline uint64_t GetRangesOffsetByIndexOfM(std::vector<Range> const& ranges,
                                          uint64_t index) {
  uint64_t offset = 0;
  for (auto const& range : ranges) {
    if (index >= range.start && index < (range.start + range.count)) {
      offset += index - range.start;
      return offset;
    }
    offset += range.count;
  }
  assert(false);
  throw std::runtime_error(__FUNCTION__);
}

inline bool CheckPhantoms(uint64_t n, std::vector<Range> const& phantoms) {
  for (auto const& phantom : phantoms) {
    if (!phantom.count || phantom.start >= n || phantom.count > n ||
        (phantom.start + phantom.count) > n)
      return false;
  }

  for (size_t i = 1; i < phantoms.size(); ++i) {
    if (phantoms[i].start <= phantoms[i - 1].start + phantoms[i - 1].count)
      return false;
  }
  return true;
}

inline bool CheckDemands(uint64_t n, std::vector<Range> const& demands) {
  if (demands.empty()) return false;

  for (auto const& demand : demands) {
    if (!demand.count || demand.start >= n || demand.count > n ||
        (demand.start + demand.count) > n) {
      return false;
    }
  }

  for (size_t i = 1; i < demands.size(); ++i) {
    if (demands[i].start <= demands[i - 1].start + demands[i - 1].count)
      return false;
  }
  return true;
}

}  // namespace scheme

namespace std {

inline std::istream& operator>>(std::istream& in, scheme::Mode& t) {
  std::string token;
  in >> token;
  if (token == "plain") {
    t = scheme::Mode::kPlain;
  } else if (token == "table") {
    t = scheme::Mode::kTable;
  } else {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, scheme::Mode const& t) {
  if (t == scheme::Mode::kPlain) {
    os << "plain";
  } else if (t == scheme::Mode::kTable) {
    os << "table";
  } else {
    os.setstate(std::ios_base::failbit);
  }
  return os;
}

inline std::istream& operator>>(std::istream& in, scheme::Action& t) {
  std::string token;
  in >> token;
  if (token == "vrf_query") {
    t = scheme::Action::kVrfQuery;
  } else if (token == "ot_vrf_query") {
    t = scheme::Action::kOtVrfQuery;
  } else if (token == "vrf_pod") {
    t = scheme::Action::kVrfPod;
  } else if (token == "ot_vrf_pod") {
    t = scheme::Action::kOtVrfPod;
  } else if (token == "complaint_pod") {
    t = scheme::Action::kComplaintPod;
  } else if (token == "ot_complaint_pod") {
    t = scheme::Action::kOtComplaintPod;
  } else if (token == "atomic_swap_pod") {
    t = scheme::Action::kAtomicSwapPod;
  } else if (token == "atomic_swap_pod_vc") {
    t = scheme::Action::kAtomicSwapPodVc;
  } else {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, scheme::Action const& t) {
  if (t == scheme::Action::kVrfQuery) {
    os << "vrf_query";
  } else if (t == scheme::Action::kOtVrfQuery) {
    os << "ot_vrf_query";
  } else if (t == scheme::Action::kVrfPod) {
    os << "vrf_pod";
  } else if (t == scheme::Action::kOtVrfPod) {
    os << "ot_vrf_pod";
  } else if (t == scheme::Action::kComplaintPod) {
    os << "complaint_pod";
  } else if (t == scheme::Action::kOtComplaintPod) {
    os << "ot_complaint_pod";
  } else if (t == scheme::Action::kAtomicSwapPod) {
    os << "atomic_swap_pod";
  } else if (t == scheme::Action::kAtomicSwapPodVc) {
    os << "atomic_swap_pod_vc";
  } else {
    os.setstate(std::ios_base::failbit);
  }
  return os;
}

inline std::istream& operator>>(std::istream& in, Range& t) {
  try {
    std::string token;
    in >> token;
    t = Range::from_string(token);  // throw
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, Range const& t) {
  std::string s = Range::to_string(t);
  os << s;
  return os;
}

}  // namespace std