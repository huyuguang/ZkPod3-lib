#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>
#include <iostream>
//#include "../pod_core/zkp_key.h"
#include "ecc.h"
#include "ecc_pub.h"
#include "pds_pub.h"
#include "groth09/groth09.h"
#include "hyrax/hyrax.h"
#include "misc.h"
#include "public.h"
#include "tick.h"
#include "groth09/test.h"
#include "vrs/test.h"

int main(int /*argc*/, char** /*argv*/) {
  InitEcc();

  std::string ecc_pub_file = "ecc_pub.bin";
  if (!OpenOrCreateEccPub(ecc_pub_file)) {
    std::cerr << "Open or create ecc pub file " << ecc_pub_file << " failed\n";
    return -1;
  }
  std::string ecc_pds_file = "pds_pub.bin";
  if (!OpenOrCreatePdsPub(ecc_pds_file)) {
    std::cerr << "Open or create pds pub file " << ecc_pds_file << " failed\n";
    return -1;
  }

  //groth09::Test();
  //vrs::Test();
  //vrs::TestLarge();
  //vrs::TestCache();

  hyrax::a1::TestRom();
  return 0;
}
