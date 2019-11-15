#include "misc.h"
#include "scheme_misc.h"
#include "scheme_ot_complaint_notary.h"
#include "tick.h"

namespace scheme::ot_complaint {

template <typename BobData>
Bob<BobData>::Bob(std::shared_ptr<BobData> b, h256_t const& self_id,
                  h256_t const& peer_id, std::vector<Range> demands,
                  std::vector<Range> phantoms)
    : b_(b),
      self_id_(self_id),
      peer_id_(peer_id),
      n_(b_->bulletin().n),
      s_(b_->bulletin().s),
      demands_(std::move(demands)),
      phantoms_(std::move(phantoms)) {
  if (!CheckDemandPhantoms(n_, demands_, phantoms_)) {
    throw std::invalid_argument("invalid demands or phantoms");
  }
  for (auto const& i : demands_) demands_count_ += i.count;
  for (auto const& i : phantoms_) phantoms_count_ += i.count;
  BuildMapping();

  bob_nonce_ = misc::RandH256();

  ot_self_pk_ = G1Rand();
  ot_beta_ = FrRand();
  ot_rand_a_ = FrRand();
  ot_rand_b_ = FrRand();
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
      map.phantom_offset = GetRangesOffsetByIndexOfM(phantoms_, i);
      ++index;
    }
  }
}

template <typename BobData>
void Bob<BobData>::GetNegoReqeust(NegoBRequest& request) {
  request.t = ot_self_pk_;
}

template <typename BobData>
bool Bob<BobData>::OnNegoRequest(NegoARequest const& request,
                                 NegoAResponse& response) {
  ot_peer_pk_ = request.s;
  response.s_exp_beta = ot_peer_pk_ * ot_beta_;
  return true;
}

template <typename BobData>
bool Bob<BobData>::OnNegoResponse(NegoBResponse const& response) {
  ot_sk_ = response.t_exp_alpha * ot_beta_;
  return true;
}

template <typename BobData>
void Bob<BobData>::GetRequest(Request& request) {
  request.bob_nonce = bob_nonce_;
  request.phantoms = phantoms_;

  request.ot_vi.reserve(demands_count_);
  for (auto const& i : demands_) {
    for (size_t j = i.start; j < i.start + i.count; ++j) {
      auto fr = MapToFr(j);
      request.ot_vi.push_back(ot_sk_ * (ot_rand_b_ * fr));
    }
  }
  request.ot_v = ot_self_pk_ * (ot_rand_a_ * ot_rand_b_);
  request_ = request;
}

template <typename BobData>
bool Bob<BobData>::OnResponse(Response response, Receipt& receipt) {
  Tick _tick_(__FUNCTION__);
  if (response.k.size() != (phantoms_count_ + 1)) {
    assert(false);
    return false;
  }
  if (response.ot_ui.size() != demands_count_) {
    assert(false);
    return false;
  }
  if (response.m.size() != phantoms_count_ * s_) {
    assert(false);
    return false;
  }

  alice_nonce_ = response.alice_nonce;
  k_ = std::move(response.k);
  ot_ui_ = std::move(response.ot_ui);
  k_mkl_root_ = CalcRootOfK(k_);

  std::vector<h256_t> seed2_h{
      {self_id_, peer_id_, alice_nonce_, bob_nonce_,
       CalcRangesDigest(phantoms_), CalcG1DataDigest(request_.ot_vi),
       G1ToBin(request_.ot_v), CalcG1DataDigest(ot_ui_), k_mkl_root_}};
  seed2_ = CalcSeed2(seed2_h);
  // std::cout << "bob seed: " << misc::HexToStr(seed2_) << "\n";

  ChainKeccak256(seed2_, phantoms_count_, w_);
  w_.push_back(FrOne());

  encrypted_m_.resize(demands_count_ * s_);

//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
//  for (int64_t i = 0; i < (int64_t)ot_ui_.size(); ++i) {
//    Fp12 e;
//    G1 ui_exp_a = ot_ui_[i] * ot_rand_a_;
//    mcl::bn256::pairing(e, ui_exp_a, ot_peer_pk_);
//    uint8_t buf[32 * 12];
//    auto ret_len = e.serialize(buf, sizeof(buf));
//    if (ret_len != sizeof(buf)) {
//      assert(false);
//      throw std::runtime_error("oops");
//    }
//    Fr fr_e = MapToFr(buf, sizeof(buf));
//
//    auto phantom_offset = mappings_[i].phantom_offset;
//    for (size_t j = 0; j < s_; ++j) {
//      encrypted_m_[i * s_ + j] = response.m[phantom_offset * s_ + j] - fr_e;
//    }
//  }
  auto parallel_f = [this, &response](int64_t i) mutable {
    Fp12 e;
    G1 ui_exp_a = ot_ui_[i] * ot_rand_a_;
    mcl::bn256::pairing(e, ui_exp_a, ot_peer_pk_);
    uint8_t buf[32 * 12];
    auto ret_len = e.serialize(buf, sizeof(buf));
    if (ret_len != sizeof(buf)) {
      assert(false);
      throw std::runtime_error("oops");
    }
    Fr fr_e = MapToFr(buf, sizeof(buf));

    auto phantom_offset = mappings_[i].phantom_offset;
    for (size_t j = 0; j < s_; ++j) {
      encrypted_m_[i * s_ + j] = response.m[phantom_offset * s_ + j] - fr_e;
    }
  };
  parallel::For((int64_t)ot_ui_.size(), parallel_f);

  if (!CheckEncryptedM()) {
    assert(false);
    return false;
  }

  receipt.count = phantoms_count_;
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
//    G1 left = sigma * w_[mapping.phantom_offset] + k_[mapping.phantom_offset];
//    G1 right = MultiExpU1(s_, [i, this](uint64_t j) -> Fr const& {
//      return encrypted_m_[i * s_ + j];
//    });
//
//    rets[i] = left == right;
//  }
  auto parallel_f = [this, &rets, &sigmas](int64_t i) mutable {
    auto const& mapping = mappings_[i];
    G1 const& sigma = sigmas[mapping.global_index];
    G1 left = sigma * w_[mapping.phantom_offset] + k_[mapping.phantom_offset];
    G1 right = MultiExpU1(s_, [i, this](uint64_t j) -> Fr const& {
      return encrypted_m_[i * s_ + j];
    });

    rets[i] = left == right;
  };
  parallel::For((int64_t)mappings_.size(), parallel_f);

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
  ChainKeccak256(secret.seed0, (phantoms_count_ + 1) * s_, v);

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
//#ifdef MULTICORE
//#pragma omp parallel for
//#endif
//  for (int64_t i = 0; i < (int64_t)mappings_.size(); ++i) {
//    auto const& mapping = mappings_[i];
//    Fr inv_w = FrInv(w_[mapping.phantom_offset]);
//    for (uint64_t j = 0; j < s_; ++j) {
//      encrypted_m_[i * s_ + j] =
//          (encrypted_m_[i * s_ + j] - v[mapping.phantom_offset * s_ + j]) *
//          inv_w;
//    }
//  }
  auto parallel_f = [this, &v](int64_t i) mutable {
    auto const& mapping = mappings_[i];
    Fr inv_w = FrInv(w_[mapping.phantom_offset]);
    for (uint64_t j = 0; j < s_; ++j) {
      encrypted_m_[i * s_ + j] =
          (encrypted_m_[i * s_ + j] - v[mapping.phantom_offset * s_ + j]) *
          inv_w;
    }
  };
  parallel::For((int64_t)mappings_.size(), parallel_f);

  decrypted_m_ = std::move(encrypted_m_);
}

template <typename BobData>
bool Bob<BobData>::SaveDecrypted(std::string const& file) {
  Tick _tick_(__FUNCTION__);
  return b_->SaveDecryped(file, demands_, decrypted_m_);
}
}  // namespace scheme::ot_complaint
