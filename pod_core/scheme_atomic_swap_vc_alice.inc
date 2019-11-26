#include "chain.h"
#include "misc.h"
#include "scheme_atomic_swap_vc_misc.h"
#include "scheme_misc.h"
#include "tick.h"

namespace scheme::atomic_swap_vc {

template <typename AliceData>
Alice<AliceData>::Alice(std::shared_ptr<AliceData> a, h256_t const& self_id,
                        h256_t const& peer_id)
    : a_(a),
      self_id_(self_id),
      peer_id_(peer_id),
      n_(a_->bulletin().n),
      s_(a_->bulletin().s) {
  alice_nonce_ = misc::RandH256();
}

template <typename AliceData>
Alice<AliceData>::~Alice() {
  if (!vrs_cache_file_.empty()) {
    if (!got_peer_receipt_) {
      vrs::ReturnCacheFile(vrs_cache_file_);
    } else {
      vrs::ExhaustCacheFile(vrs_cache_file_);
    }
  }
}

template <typename AliceData>
void Alice<AliceData>::BuildMapping() {
  Tick _tick_(__FUNCTION__);
  mappings_.resize(demands_count_);
  size_t index = 0;
  for (auto const& p : demands_) {
    for (size_t i = p.start; i < (p.start + p.count); ++i) {
      mappings_[index++].global_index = i;
    }
  }
}

template <typename AliceData>
bool Alice<AliceData>::OnRequest(Request request, Response& response) {
  Tick _tick_(__FUNCTION__);

  if (!CheckDemands(n_, request.demands)) {
    assert(false);
    return false;
  }

  demands_count_ =
      std::accumulate(request.demands.begin(), request.demands.end(), 0ULL,
                      [](uint64_t t, Range const& i) { return t + i.count; });

  demands_ = std::move(request.demands);
  bob_nonce_ = request.bob_nonce;
  response.alice_nonce = alice_nonce_;
  BuildMapping();

  char const* data_dir = std::getenv("options:data_dir");
  std::string cache_dir = data_dir ? data_dir : ".";
  cache_dir += "/vrs_cache";
  auto vrs_count = (demands_count_ + 1) * s_;

  vrs::CacheUPtr cache;
  vrs_cache_file_ = vrs::SelectCacheFile(cache_dir, vrs_count);
  if (!vrs_cache_file_.empty()) {
    cache.reset(new vrs::Cache());
    if (!vrs::LoadCache(vrs_cache_file_, *cache, false)) {
      std::cerr << "LoadCache failed\n";
      cache.reset();
      vrs_cache_file_.clear();
    } else {
      vrs::UpgradeCache(*cache, vrs_count);
    }
  }

  std::vector<std::vector<G1>> cached_var_coms;
  std::vector<std::vector<Fr>> cached_var_coms_r;
  if (cache) {
    response.vrs_plain_seed = cache->seed;
    seed0_ = cache->key;
    seed0_com_r_ = cache->key_com_r;
    cached_var_coms = std::move(cache->var_coms);
    cached_var_coms_r = std::move(cache->var_coms_r);
  } else {
    response.vrs_plain_seed = misc::RandH256();
    seed0_ = FrRand();
    seed0_com_r_ = FrRand();
  }

  std::vector<Fr> plain((demands_count_ + 1) * s_);
  auto parallel_f = [&plain, &response](uint64_t i) {
    plain[i] = vrs::GeneratePlain(response.vrs_plain_seed, i);
  };
  parallel::For(plain.size(), parallel_f);

  v_.resize(plain.size());
  auto parallel_f_v = [this, &plain](uint64_t i) {
    v_[i] = vrs::Mimc5Enc(plain[i], seed0_);
  };
  parallel::For(plain.size(), parallel_f_v);

  if (evil_) {
    // NOTE: use rand() for test
    uint64_t evil_i = rand() % demands_count_;
    uint64_t evil_j = s_ - 1;  // last col
    v_[evil_i * s_ + evil_j] = FrRand();
    std::cout << "evil: " << evil_i << "," << evil_j << "\n";
  }

  BuildK(v_, response.k, s_);

  std::vector<h256_t> seed2_h{{peer_id_, self_id_, alice_nonce_, bob_nonce_,
                               CalcRangesDigest(demands_),
                               CalcRootOfK(response.k)}};
  seed2_ = CalcSeed2(seed2_h);

  ChainKeccak256(seed2_, demands_count_, w_);
  w_.push_back(FrOne());

  // compute mij' = vij + wi * mij
  auto const& m = a_->m();
  response.m.resize(demands_count_ * s_);

  auto parallel_f2 = [this, &response, &m](uint64_t i) mutable {
    auto const& map = mappings_[i];
    auto is = i * s_;
    auto m_is = map.global_index * s_;
    for (uint64_t j = 0; j < s_; ++j) {
      auto ij = is + j;
      auto m_ij = m_is + j;
      response.m[ij] = v_[ij] + w_[i] * m[m_ij];
    }
  };
  parallel::For(mappings_.size(), parallel_f2);

  response.vw.resize(s_);
  auto parallel_f3 = [this, &response](int64_t j) mutable {
    response.vw[j] = FrZero();
    for (size_t i = 0; i <= demands_count_; ++i) {
      response.vw[j] += v_[i * s_ + j] * w_[i];
    }
  };
  parallel::For(s_, parallel_f3);

  // for (size_t j = 0; j < s_; ++j) {
  //  response.vw[j] = FrZero();
  //  for (size_t i = 0; i <= demands_count_; ++i) {
  //    response.vw[j] += v_[i * s_ + j] * w_[i];
  //  }
  //}

  // use seed2_ as rom_seed
  vrs::PublicInput public_input(plain.size(),
                                [&plain](int64_t i) { return plain[i]; });

  vw_com_r_ = FrRand();
  vrs::SecretInput secret_input(seed0_, seed0_com_r_, vw_com_r_);
  vrs::LargeProverLowRam prover(public_input, secret_input,
                                std::move(cached_var_coms),
                                std::move(cached_var_coms_r));
  auto get_w = [this](int64_t i) { return w_[i / s_]; };
  vrs::ProveOutput vrs_output;
  prover.Prove(seed2_, get_w, response.vrs_proofs, vrs_output);
  response.vw_com_r = vw_com_r_;

  receipt_.h = vrs_output.h;
  receipt_.g = vrs_output.g;
  receipt_.seed0_com = vrs_output.key_com;

#ifdef _DEBUG
  if (!evil_) {
    assert(v_ == prover.v());
    assert(parallel::Accumulate(response.vw.begin(), response.vw.end(),
                                FrZero()) == prover.vw());
  }
#endif

  return true;
}

template <typename AliceData>
bool Alice<AliceData>::OnReceipt(Receipt const& receipt, Secret& secret) {
  if (receipt != receipt_) {
    assert(false);
    return false;
  }

  secret.seed0 = seed0_;
  secret.seed0_r = seed0_com_r_;

  got_peer_receipt_ = true;

  return true;
}

}  // namespace scheme::atomic_swap_vc
