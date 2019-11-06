#include "misc.h"
#include "tick.h"

namespace scheme::complaint {

template <typename BobData>
Bob<BobData>::Bob(std::shared_ptr<BobData> b, h256_t const& self_id,
                  h256_t const& peer_id, std::vector<Range> demands)
    : b_(b),
      self_id_(self_id),
      peer_id_(peer_id),
      n_(b_->bulletin().n),
      s_(b_->bulletin().s),
      demands_(std::move(demands)) {

  if (!CheckDemands(n_, demands_)) {
    throw std::invalid_argument("invalid demands");
  }

  demands_count_ =
      std::accumulate(demands_.begin(), demands_.end(), 0ULL,
                      [](uint64_t a, Range const& b) { return a + b.count; });

  BuildMapping();

  bob_nonce_ = misc::RandH256();
}

template <typename BobData>
void Bob<BobData>::BuildMapping() {
  Tick _tick_(__FUNCTION__);
  mappings_.resize(demands_count_);
  size_t index = 0;
  for (auto const& d : demands_) {
    for (size_t i = d.start; i < (d.start + d.count); ++i) {
      auto& map = mappings_[index];
      map.global_index = i;
      ++index;
    }
  }
}

template <typename BobData>
void Bob<BobData>::GetRequest(Request& request) {
  request.bob_nonce = bob_nonce_;
  request.demands = demands_;
}

template <typename BobData>
bool Bob<BobData>::OnResponse(Response response, Receipt& receipt) {
  Tick _tick_(__FUNCTION__);

  if (response.k.size() != (demands_count_ + 1)) {
    assert(false);
    return false;
  }
  if (response.m.size() != demands_count_ * s_) {
    assert(false);
    return false;
  }

  k_ = std::move(response.k);

  k_mkl_root_ = CalcRootOfK(k_);

  alice_nonce_ = response.alice_nonce;

  std::vector<h256_t> seed2_h{{self_id_, peer_id_, alice_nonce_, bob_nonce_,
                               CalcRangesDigest(demands_), k_mkl_root_}};
  seed2_ = CalcSeed2(seed2_h);

  ChainKeccak256(seed2_, demands_count_, w_);
  w_.push_back(FrOne());

  encrypted_m_ = std::move(response.m);

  if (!CheckEncryptedM()) {
    assert(false);
    return false;
  }

  receipt.count = demands_count_;
  receipt.k_mkl_root = k_mkl_root_;
  receipt.seed2 = seed2_;

  receipt_ = receipt;
  return true;
}

template <typename BobData>
bool Bob<BobData>::CheckEncryptedM() {
  Tick _tick_(__FUNCTION__);

  auto const& sigmas = b_->sigmas();

  // uint64_t phantom_offset = phantom_.start - demand_.start;
  std::vector<int64_t> rets(mappings_.size());
//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
//  for (int64_t i = 0; i < (int64_t)mappings_.size(); ++i) {
//    auto const& mapping = mappings_[i];
//    G1 const& sigma = sigmas[mapping.global_index];
//    G1 left = sigma * w_[i] + k_[i];
//    G1 right = MultiExpU1(s_, [i, this](uint64_t j) -> Fr const& {
//      return encrypted_m_[i * s_ + j];
//    });
//
//    rets[i] = left == right;
//  }
  auto parallel_f = [this, &rets, &sigmas](int64_t i) mutable {
    auto const& mapping = mappings_[i];
    G1 const& sigma = sigmas[mapping.global_index];
    G1 left = sigma * w_[i] + k_[i];
    G1 right = MultiExpU1(s_, [i, this](uint64_t j) -> Fr const& {
      return encrypted_m_[i * s_ + j];
    });

    rets[i] = left == right;
  };
  parallel::For((int64_t)mappings_.size(), parallel_f, "CheckEncryptedM");

  if (std::any_of(rets.begin(), rets.end(), [](int64_t r) { return !r; })) {
    assert(false);
    std::cerr << "ASSERT: " << __FUNCTION__ << ": " << __LINE__ << "\n";
    return false;
  }
  return true;
}

template <typename BobData>
bool Bob<BobData>::OnSecret(Secret const& secret) {
  Tick _tick_(__FUNCTION__);

  // compute v
  std::vector<Fr> v;
  ChainKeccak256(secret.seed0, (demands_count_ + 1) * s_, v);

  if (!CheckK(v)) {
    assert(claim_i_ >= 0);
    return false;
  } else {
#ifdef _DEBUG
    Claim claim;
    BuildClaim(rand() % k_.size(), claim);
    assert(!VerifyClaim(s_, receipt_, secret, claim));
#endif

    DecryptM(v);
    return true;
  }
}

template <typename BobData>
bool Bob<BobData>::GenerateClaim(Claim& claim) {
  if (claim_i_ == -1) {
    assert(false);
    return false;
  }
  BuildClaim(claim_i_, claim);
  return true;
}

template <typename BobData>
bool Bob<BobData>::CheckK(std::vector<Fr> const& v) {
  Tick _tick_(__FUNCTION__);

  // compute k
  std::vector<G1> k;
  BuildK(v, k, s_);

  assert(k.size() == k_.size());

  // compare k
  for (int64_t i = 0; i < (int64_t)k.size(); ++i) {
    if (k[i] != k_[i]) {
      claim_i_ = i;
      return false;
    }
  }
  return true;
}

template <typename BobData>
void Bob<BobData>::BuildClaim(uint64_t i, Claim& claim) {
  Tick _tick_(__FUNCTION__);
  claim.i = i;
  claim.ki = k_[i];
  auto root = CalcPathOfK(k_, i, claim.mkl_path);
  if (root != k_mkl_root_) {
    assert(false);
    throw std::runtime_error("oops, mkl root mismatch");
  }
}

template <typename BobData>
void Bob<BobData>::DecryptM(std::vector<Fr> const& v) {
  Tick _tick_(__FUNCTION__);

  std::vector<Fr> inv_w = w_;
  FrInv(inv_w.data(), inv_w.size());

//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
//  for (int64_t i = 0; i < (int64_t)mappings_.size(); ++i) {
//    auto is = i * s_;
//    for (uint64_t j = 0; j < s_; ++j) {
//      auto ij = is + j;
//      encrypted_m_[ij] = (encrypted_m_[ij] - v[ij]) * inv_w[i];
//    }
//  }
  auto parallel_f = [this, &v, &inv_w](int64_t i) mutable {
    auto is = i * s_;
    for (uint64_t j = 0; j < s_; ++j) {
      auto ij = is + j;
      encrypted_m_[ij] = (encrypted_m_[ij] - v[ij]) * inv_w[i];
    }
  };
  parallel::For((int64_t)mappings_.size(), parallel_f, "DecryptM");

  decrypted_m_ = std::move(encrypted_m_);
}

template <typename BobData>
bool Bob<BobData>::SaveDecrypted(std::string const& file) {
  Tick _tick_(__FUNCTION__);
  return b_->SaveDecryped(file, demands_, decrypted_m_);
}
}  // namespace scheme::complaint
