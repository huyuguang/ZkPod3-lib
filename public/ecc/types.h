#pragma once

#include "public.h"

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4191)
#endif

#include <mcl/bn256.hpp>
#include <mcl/window_method.hpp>

#ifdef _WIN32
#pragma warning(pop)
#endif

using mcl::bn256::Fp;
using mcl::bn256::Fp12;
using mcl::bn256::Fp2;
using mcl::bn256::Fp6;
using mcl::bn256::Fr;
using mcl::bn256::G1;
using mcl::bn256::G2;
typedef mcl::fp::WindowMethod<G1> G1WM;
typedef mcl::fp::WindowMethod<G2> G2WM;

typedef std::shared_ptr<G1> G1Ptr;
typedef std::shared_ptr<Fr> FrPtr;

enum {
  kFrBinSize = 32,
  kFpBinSize = 32,
  kG1FlatBinSize = 1 + kFpBinSize * 2,  // x,y and 0
  kFp2BinSize = 64,
  kG2FlatBinSize = 1 + kFp2BinSize * 2,
  kG1CompBinSize = 32,
  kG2CompBinSize = 64,
};
