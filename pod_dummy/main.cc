#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>

#include <iostream>

#include "cmd/test.h"
#include "ecc/ecc.h"
#include "groth09/test.h"
#include "hyrax/test.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "pc_utils/test.h"
#include "pod/test.h"
#include "public.h"
#include "vrs/sha256c_gadget.h"
#include "vrs/test.h"

bool InitAll(std::string const& data_dir) {
  InitEcc();

  // auto ecc_pub_file = data_dir + "/" + "ecc_pub.bin";
  // if (!OpenOrCreateEccPub(ecc_pub_file)) {
  //  std::cerr << "Open or create ecc pub file " << ecc_pub_file << "
  //  failed\n"; return false;
  //}

  auto ecc_pds_file = data_dir + "/" + "pds_pub.bin";
  if (!OpenOrCreatePdsPub(ecc_pds_file)) {
    std::cerr << "Open or create pds pub file " << ecc_pds_file << " failed\n";
    return false;
  }

  return true;
}

enum VrsSchemeType { kMimic5 = 0, kSha256c = 1 };

struct ParamIntPair {
  bool valid() const { return x && y; }
  int64_t x = 0;
  int64_t y = 0;
};

inline std::istream& operator>>(std::istream& in, ParamIntPair& t) {
  try {
    char sperator;
    in >> t.x;
    in >> sperator;
    in >> t.y;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}
//
// inline std::ostream& operator<<(std::ostream& os, ParamIntPair const& t) {
//  char sperator = '*';
//  os << t.x;
//  os << sperator;
//  os << t.y;
//  return os;
//}

int main(int argc, char** argv) {
  setlocale(LC_ALL, "");
  std::string data_dir;
  uint32_t thread_num;
  int vrs_scheme;
  int64_t hyrax_a1_n = 0;
  int64_t hyrax_a2_n = 0;
  int64_t hyrax_a3_n = 0;
  int64_t gro09_51a_n = 0;
  int64_t gro09_51b_n = 0;
  ParamIntPair gro09_52a;
  ParamIntPair gro09_52b;
  ParamIntPair gro09_53a;
  ParamIntPair gro09_53b;
  ParamIntPair gro09_43b;
  int64_t vrs_basic_n = 0;
  int64_t vrs_large_n = 0;
  ParamIntPair pod;
  ParamIntPair matrix;
  ParamIntPair equal_ip;
  bool overlap = false;
  bool divide = false;
  bool match = false;
  bool substr = false;
  bool pack = false;
  bool substrpack = false;
  bool matchpack = false;
  bool match_query = false;
  bool substr_query = false;

  try {
    po::options_description options("command line options");
    options.add_options()("help,h", "Use -h or --help to list all arguments")(
        "data_dir,d", po::value<std::string>(&data_dir)->default_value("."),
        "Provide the data dir")(
        "thread_num", po::value<uint32_t>(&thread_num)->default_value(0),
        "Provide the number of the parallel thread, 1: disable, 0: default.")(
        "vrs_scheme", po::value<int>(&vrs_scheme)->default_value(kMimic5),
        "Provide the scheme type, 0: mimc5, 1:sha256c")(
        "hyrax_a1", po::value<int64_t>(&hyrax_a1_n), "")(
        "hyrax_a2", po::value<int64_t>(&hyrax_a2_n), "")(
        "hyrax_a3", po::value<int64_t>(&hyrax_a3_n), "")(
        "gro09_51a", po::value<int64_t>(&gro09_51a_n), "")(
        "gro09_51b", po::value<int64_t>(&gro09_51b_n), "")(
        "gro09_52a", po::value<ParamIntPair>(&gro09_52a), "m*n, ex: 10*20")(
        "gro09_52b", po::value<ParamIntPair>(&gro09_52b), "m*n, ex: 10*20")(
        "gro09_53a", po::value<ParamIntPair>(&gro09_53a), "m*n, ex: 10*20")(
        "gro09_53b", po::value<ParamIntPair>(&gro09_53b), "m*n, ex: 10*20")(
        "gro09_43b", po::value<ParamIntPair>(&gro09_43b), "m*n, ex: 10*20")(
        "vrs_basic", po::value<int64_t>(&vrs_basic_n), "")(
        "vrs_large", po::value<int64_t>(&vrs_large_n), "")(
        "pod", po::value<ParamIntPair>(&pod), "n*s, ex: 100*1023")(
        "matrix", po::value<ParamIntPair>(&matrix), "n*s, ex:10*20")(
        "equal_ip", po::value<ParamIntPair>(&equal_ip), "xn,yn, ex:10,20")(
        "overlap", "")(
        "divide", "")(
        "match", "")(
        "substr", "")(
        "pack", "")(
        "substrpack", "")(
        "matchpack", "")(
        "match_query", "")(
        "substr_query", "");

    boost::program_options::variables_map vmap;

    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, options), vmap);
    boost::program_options::notify(vmap);

    if (vmap.count("help")) {
      std::cout << options << std::endl;
      return -1;
    }

    if (vrs_scheme != kMimic5 && vrs_scheme != kSha256c) {
      std::cout << options << std::endl;
      return -1;
    }

    if (vmap.count("overlap")) {
      overlap = true;
    }

    if (vmap.count("divide")) {
      divide = true;
    }

    if (vmap.count("match")) {
      match = true;
    }

    if (vmap.count("substr")) {
      substr = true;
    }

    if (vmap.count("pack")) {
      pack = true;
    }

    if (vmap.count("substrpack")) {
      substrpack = true;
    }

    if (vmap.count("matchpack")) {
      matchpack = true;
    }

    if (vmap.count("match_query")) {
      match_query = true;
    }

    if (vmap.count("substr_query")) {
      substr_query = true;
    }
  } catch (std::exception& e) {
    std::cout << "Unknown parameters.\n"
              << e.what() << "\n"
              << "-h or --help to list all arguments.\n";
    return -1;
  }

  int tbb_thread_num =
      thread_num ? (int)thread_num : tbb::task_scheduler_init::automatic;
  tbb::task_scheduler_init init(tbb_thread_num);

  parallel::CheckAllocationHook();

  if (!InitAll(data_dir)) {
    std::cerr << "Init failed\n";
    return -1;
  }

  std::map<std::string, bool> rets;

  if (hyrax_a1_n) {
    rets["hyrax::a1"] = hyrax::a1::TestRom();
  }
  if (hyrax_a2_n) {
    rets["hyrax::a2"] = hyrax::a2::TestRom(hyrax_a2_n);
  }
  if (hyrax_a3_n) {
    rets["hyrax::a3"] = hyrax::a3::TestRom(hyrax_a3_n);
  }
  if (gro09_51a_n) {
    rets["groth09::sec51a"] = groth09::sec51a::TestRom(gro09_51a_n);
  }
  if (gro09_51b_n) {
    rets["groth09::sec51b"] = groth09::sec51b::TestRom(gro09_51b_n);
  }
  if (gro09_52a.valid()) {
    rets["groth09::sec52a"] =
        groth09::sec52a::TestRom(gro09_52a.x, gro09_52a.y);
  }
  if (gro09_52b.valid()) {
    rets["groth09::sec52b"] =
        groth09::sec52b::TestRom(gro09_52b.x, gro09_52b.y);
  }
  if (gro09_53a.valid()) {
    rets["groth09::sec53a"] =
        groth09::sec52a::TestRom(gro09_53a.x, gro09_53a.y);
  }
  if (gro09_53b.valid()) {
    rets["groth09::sec53b"] =
        groth09::sec52b::TestRom(gro09_53b.x, gro09_53b.y);
  }
  if (gro09_43b.valid()) {
    rets["groth09::sec43b"] =
        groth09::sec43b::TestRom(gro09_43b.x, gro09_43b.y);
  }
  if (matrix.valid()) {
    rets["pc_utils::matrix"] = pc_utils::matrix::Test(matrix.x, matrix.y);
  }
  if (equal_ip.valid()) {
    rets["pc_utils::equal_ip"] = pc_utils::equal_ip::Test(equal_ip.x, equal_ip.y);
  }
  if (overlap) {
    rets["pc_utils::overlap"] = pc_utils::overlap::Test();
  }
  if (divide) {
    rets["pc_utils::divide"] = pc_utils::divide::Test();
  }
  if (match) {
    rets["pc_utils::match"] = pc_utils::match::Test();
  }
  if (substr) {
    rets["pc_utils::substr"] = pc_utils::substr::Test();
  }
  if (pack) {
    rets["pc_utils::pack"] = pc_utils::pack::Test();
  }
  if (substrpack) {
    rets["pc_utils::substrpack"] = pc_utils::substrpack::Test();
  }
  if (matchpack) {
    rets["pc_utils::matchpack"] = pc_utils::matchpack::Test();
  }
  if (match_query) {
    rets["cmd::match_query"] = cmd::match_query::Test();
  }
  if (substr_query) {
    rets["cmd::substr_query"] = cmd::substr_query::Test();
  }

  if (vrs_scheme == VrsSchemeType::kMimic5) {
    if (vrs_basic_n) {
      rets["vrs::basic"] = vrs::TestBasic<vrs::Mimc5Scheme>(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::large"] = vrs::TestLarge<vrs::Mimc5Scheme>(vrs_basic_n);
    }
    if (pod.valid()) {
      rets["pod"] = pod::TestBasic<vrs::Mimc5Scheme>(pod.x, pod.y);
    }
  } else {
    if (vrs_basic_n) {
      rets["vrs::basic"] = vrs::TestBasic<vrs::Sha256cScheme>(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::large"] = vrs::TestLarge<vrs::Sha256cScheme>(vrs_basic_n);
    }
    if (pod.valid()) {
      rets["pod"] = pod::TestBasic<vrs::Sha256cScheme>(pod.x, pod.y);
    }
  }

  std::cout << "\n=============================\n";
  std::cout << rets.size() << " test finished\n";
  bool all_success = true;
  for (auto const& i : rets) {
    std::cout << i.first << ": " << (i.second ? "success" : "failed") << "\n";
    if (!i.second) all_success = false;
  }

  return all_success ? 0 : 1;
}
