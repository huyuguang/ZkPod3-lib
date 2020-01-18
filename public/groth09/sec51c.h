#pragma once

#include "groth09/details.h"
#include "utils/fst.h"
#include "hyrax/a3.h"

// t: public vector<Fr>, size = n
// X, Y: secret vector<Fr>, size = n
// z: secret Fr
// open: com(gx, X), com(gy, Y), com(gz, z)
// prove: z = <X, Y o t>
// proof size: 2log(n) Fr and 4 G1
// prove cost: 2*mulexp(n)
// verify cost: mulexp(n)
namespace groth09::sec51c {

class ProverInput {
 private:
  std::vector<Fr> const& x_;   // size = n
  std::vector<Fr> const& y_;   // size = n
  std::vector<Fr> const& t_;   // size = n
  std::vector<Fr> const& yt_;  // yt = y o t
  Fr z_;                       // z = <x, y o t>

 public:
  std::vector<Fr> const& t() const { return t_; }
  Fr const& z() const { return z_; }
  std::vector<Fr> const& x() const { return x_; }
  std::vector<Fr> const& y() const { return y_; }
  std::vector<Fr> const& yt() const { return yt_; }
  Fr const& x(size_t i) const { return x_[i]; }
  Fr const& y(size_t i) const { return y_[i]; }
  Fr const& yt(size_t i) const { return yt_[i]; }
  int64_t n() const { return x().size(); }

 public:
  ProverInput(std::vector<Fr> const& x, std::vector<Fr> const& y,
              std::vector<Fr> const& t, std::vector<Fr> const& yt, Fr const& z)
      : x_(x), y_(y), t_(t), yt_(yt), z_(z) {
    assert(!x.empty());
    assert(x.size() == y.size());
    assert(x.size() == t.size());
    assert(x.size() == yt.size());
    assert(yt == details::HadamardProduct(y, t));
    assert(z == InnerProduct(x, yt));
  }
};

struct CommitmentPub {
  CommitmentPub() {}
  CommitmentPub(G1 const& a, G1 const& b, G1 const& c) : a(a), b(b), c(c) {}
  G1 a;
  G1 b;
  G1 c;
};

struct CommitmentSec {
  CommitmentSec() {}
  CommitmentSec(Fr const& r, Fr const& s, Fr const& t) : r(r), s(s), t(t) {}
  Fr r;
  Fr s;
  Fr t;
};

struct RomProof {
  G1 com_v;
  hyrax::a3::RomProof proof_a3;
};
}  // namespace groth09::sec51c
