#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>

#include <iostream>

#include "ecc/ecc.h"
#include "groth09/groth09.h"
#include "groth09/test.h"
#include "hyrax/hyrax.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "pc_utils/test.h"
#include "pod/test.h"
#include "public.h"
#include "vrs/test.h"
#include "cmd/cmd.h"


int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

#ifdef _DEBUG
  int thread_num = 1;  // disable parallel
#else
  int thread_num = 0;
#endif

#ifdef USE_TBB
  int tbb_thread_num =
      thread_num ? (int)thread_num : tbb::task_scheduler_init::automatic;
  tbb::task_scheduler_init init(tbb_thread_num);
  std::cout << "use tbb, thread_num " << tbb_thread_num << "\n";
#else
  setenv("options:thread_num", std::to_string(thread_num).c_str(), true);
#endif

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

  // TestMultiexp();

  // groth09::Test();
  // vrs::Test();
  // vrs::TestLarge();
  // vrs::TestCache();
  // hyrax::a1::TestRom();
  // pod::Test();
  //pc_utils::Test();
  //cmd::substr_query::Test();
  cmd::match_query::Test();
  return 0;
}
