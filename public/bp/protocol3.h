#pragma once

#include "./protocol1.h"

// protocol 3
// a, b: secret vector<Fr>
// c = <a,b>
// open com(a,b): P=h*alpha + g1*a + g2*b
// open com(c): Q=h*beta + u * c
// prove: c = <a,b>

namespace bp {

struct Protocol3Proof {};

//Protocol3Proof Protocol3Prove(h256_t seed, std::vector<G1>&& g,
//                              std::vector<G1>&& h, std::vector<Fr>&& a,
//                              std::vector<Fr>&& b, G1 const& p, G1 const& u,
//                              Fr const& c, G1 const& q) {}
}