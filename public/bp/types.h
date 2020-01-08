#pragma once

#include "./details.h"

namespace bp {

class Challenge {
 public:
  Challenge(uint64_t count, uint64_t seed, bool zero_u) : count_(count) {
    auto seed_big = boost::endian::native_to_big(seed);
    Init((uint8_t const*)&seed_big, sizeof(seed_big), zero_u);
  }
  Challenge(uint64_t count, h256_t const& seed, bool zero_u) : count_(count) {
    Init(seed.data(), seed.size(), zero_u);
  }
  Challenge(uint64_t count, uint8_t const* seed, uint64_t size, bool zero_u)
      : count_(count) {
    Init(seed, size, zero_u);
  }

  uint64_t count() const { return count_; }
  std::vector<Fr> const& x() const { return x_; }
  std::vector<Fr> const& x_inverse() const { return x_inverse_; }
  std::vector<Fr> const& x_square() const { return x_square_; }
  std::vector<Fr> const& x_square_inverse() const { return x_square_inverse_; }
  G1 const& u() const { return u_; }
  bool zero_u() const { return u_ == G1Zero(); }

 private:
  void Init(uint8_t const* seed, uint64_t size, bool zero_u) {
    Tick tick(__FUNCTION__);
    using namespace details;
    auto x_count = PackXCount(count_);
    x_.resize(x_count);
    x_inverse_.resize(x_count);
    x_square_.resize(x_count);
    x_square_inverse_.resize(x_count);

    std::string str_seed((char const*)seed, size);
    for (uint64_t i = 0; i < x_count; ++i) {
      auto str_x = str_seed + "_" + std::to_string(i) + "_x";
      x_[i] = StrHashToFr(str_x);
      Fr::inv(x_inverse_[i], x_[i]);
      Fr::sqr(x_square_[i], x_[i]);
      Fr::inv(x_square_inverse_[i], x_square_[i]);
    }

    if (zero_u) {
      u_ = G1Zero();
    } else {
      auto str_u = str_seed + "_u";
      u_ = MapToG1(str_u);
    }
  }

 private:
  uint64_t const count_;
  std::vector<Fr> x_;
  std::vector<Fr> x_inverse_;
  std::vector<Fr> x_square_;
  std::vector<Fr> x_square_inverse_;
  G1 u_;
};

struct P1Committment {
  G1 q(G1 const& u) const { return p + u * c; }
  G1 p;
  Fr c;

  bool operator==(P1Committment const& v) const { return p == v.p && c == v.c; }
  bool operator!=(P1Committment const& v) const { return !((*this) == v); }

  constexpr size_t GetBufSize() const { return kG1FlatBinSize + kFrBinSize; }

  size_t serialize(void* buf, uint64_t max_size) const {
    if (max_size < GetBufSize()) return 0;
    auto pbuf = (uint8_t*)buf;
    G1ToFlatBin(p, pbuf);
    FrToBin(c, pbuf + kG1FlatBinSize);
    return GetBufSize();
  }

  size_t deserialize(void const* buf, uint64_t buf_size) {
    if (buf_size < GetBufSize()) return 0;
    auto pbuf = (uint8_t const*)buf;
    if (!FlatBinToG1(pbuf, &p)) return 0;
    if (!BinToFr32(pbuf + kG1FlatBinSize, &c)) return 0;
    return GetBufSize();
  }
};

struct P2Proof {
  G1 q;
  std::vector<G1> left;  // size() == log(count)
  std::vector<G1> right;
  Fr a;
  Fr b;

  bool operator==(P2Proof const& v) const {
    return q == v.q && a == v.a && b == v.b && left == v.left &&
           right == v.right;
  }
  bool operator!=(P2Proof const& v) const { return !((*this) == v); }

  size_t GetBufSize() const {
    assert(left.size() == right.size());
    return kG1FlatBinSize + kFrBinSize + kFrBinSize + sizeof(uint64_t) +
           2 * kG1FlatBinSize * left.size();
  }

  size_t serialize(void* buf, uint64_t max_size) const {
    assert(left.size() == right.size());

    size_t buf_size = GetBufSize();

    if (max_size < buf_size) return 0;

    auto pbuf = (uint8_t*)buf;

    G1ToFlatBin(q, pbuf);
    pbuf += kG1FlatBinSize;

    FrToBin(a, pbuf);
    pbuf += kFrBinSize;

    FrToBin(b, pbuf);
    pbuf += kFrBinSize;

    uint64_t count = left.size();
    count = boost::endian::native_to_big(count);
    memcpy(pbuf, &count, sizeof(count));
    pbuf += sizeof(count);

    for (auto const& i : left) {
      G1ToFlatBin(i, pbuf);
      pbuf += kG1FlatBinSize;
    }
    for (auto const& i : right) {
      G1ToFlatBin(i, pbuf);
      pbuf += kG1FlatBinSize;
    }

    assert(buf_size == (uint64_t)(pbuf - (uint8_t*)buf));

    return buf_size;
  }

  size_t deserialize(void const* buf, uint64_t buf_size) {
    if (buf_size <=
        (kG1FlatBinSize + kFrBinSize + kFrBinSize + sizeof(uint64_t)))
      return 0;

    auto pbuf = (uint8_t const*)buf;

    if (!FlatBinToG1(pbuf, &q)) return 0;
    pbuf += kG1FlatBinSize;

    if (!BinToFr32(pbuf, &a)) return 0;
    pbuf += kFrBinSize;

    if (!BinToFr32(pbuf, &b)) return 0;
    pbuf += kFrBinSize;

    uint64_t count;
    memcpy(&count, pbuf, sizeof(count));
    count = boost::endian::big_to_native(count);
    pbuf += sizeof(count);

    size_t data_size = kG1FlatBinSize + kFrBinSize + kFrBinSize +
                       sizeof(uint64_t) + 2 * kG1FlatBinSize * count;
    if (buf_size < data_size) return 0;

    left.resize(count);
    for (auto& i : left) {
      if (!FlatBinToG1(pbuf, &i)) return 0;
      pbuf += kG1FlatBinSize;
    }

    right.resize(count);
    for (auto& i : right) {
      if (!FlatBinToG1(pbuf, &i)) return 0;
      pbuf += kG1FlatBinSize;
    }

    return buf_size;
  }
};

struct P1Proof {
  P1Committment committment;
  P2Proof p2_proof;

  bool operator==(P1Proof const& v) const {
    return committment == v.committment && p2_proof == v.p2_proof;
  }
  bool operator!=(P1Proof const& v) const { return !((*this) == v); }

  size_t GetBufSize() const {
    return committment.GetBufSize() + p2_proof.GetBufSize();
  }

  size_t serialize(void* buf, uint64_t max_size) const {
    auto pbuf = (uint8_t*)buf;
    auto size1 = committment.serialize(pbuf, max_size);
    if (!size1) return 0;
    auto size2 = p2_proof.serialize(pbuf + size1, max_size - size1);
    if (!size2) return 0;
    auto size = size1 + size2;
    assert(size == GetBufSize());
    return size;
  }

  size_t deserialize(void const* buf, uint64_t buf_size) {
    auto pbuf = (uint8_t const*)buf;
    auto size1 = committment.deserialize(pbuf, buf_size);
    if (!size1) return 0;
    auto size2 = p2_proof.deserialize(pbuf + size1, buf_size - size1);
    if (!size2) return 0;
    auto size = size1 + size2;
    assert(size == GetBufSize());
    return size;
  }
};
}  // namespace bp