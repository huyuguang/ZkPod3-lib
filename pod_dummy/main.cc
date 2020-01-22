#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>

#include <iostream>

#include "ecc/ecc.h"
#include "groth09/test.h"
#include "hyrax/test.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "pc_utils/test.h"
#include "pod/test.h"
#include "public.h"
#include "vrs/test.h"
#include "cmd/test.h"
#include "vrs/sha256c_gadget.h"

bool InitAll(std::string const& data_dir) {
  InitEcc();

  //auto ecc_pub_file = data_dir + "/" + "ecc_pub.bin";
  //if (!OpenOrCreateEccPub(ecc_pub_file)) {
  //  std::cerr << "Open or create ecc pub file " << ecc_pub_file << " failed\n";
  //  return false;
  //}

  auto ecc_pds_file = data_dir + "/" + "pds_pub.bin";
  if (!OpenOrCreatePdsPub(ecc_pds_file)) {
    std::cerr << "Open or create pds pub file " << ecc_pds_file << " failed\n";
    return false;
  }

  return true;
}

enum VrsSchemeType { kMimic5 = 0, kSha256c = 1 };

struct Gro09Param {
  int64_t m = 0;
  int64_t n = 0;
};
inline std::istream& operator>>(std::istream& in, Gro09Param& t) {
  try {
    char sperator;
    in >> t.m;
    in >> sperator;
    in >> t.n;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, Gro09Param const& t) {
  char sperator = '*';
  os << t.m;
  os << sperator;
  os << t.n;
  return os;
}

struct PodParam {
  int64_t n = 0;
  int64_t s = 0;
};

inline std::istream& operator>>(std::istream& in, PodParam& t) {
  try {
    char sperator;
    in >> t.n;
    in >> sperator;
    in >> t.s;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

inline std::ostream& operator<<(std::ostream& os, PodParam const& t) {
  char sperator = '*';
  os << t.n;
  os << sperator;
  os << t.s;
  return os;
}


int main(int argc, char** argv) {
  setlocale(LC_ALL, "");
  std::string data_dir;
  uint32_t thread_num;
  int vrs_scheme;
  int64_t hyrax_a1_n;
  int64_t hyrax_a2_n;
  int64_t hyrax_a3_n;
  int64_t gro09_51a_n;
  int64_t gro09_51b_n;
  Gro09Param gro09_52a;
  Gro09Param gro09_52b;
  Gro09Param gro09_53a;
  Gro09Param gro09_53b;
  Gro09Param gro09_43b;
  int64_t vrs_basic_n;
  int64_t vrs_large_n;
  PodParam pod;

  try {
    po::options_description options("command line options");
    options.add_options()("help,h", "Use -h or --help to list all arguments")(
        "data_dir,d", po::value<std::string>(&data_dir)->default_value("."),
        "Provide the data dir")(
        "thread_num", po::value<uint32_t>(&thread_num)->default_value(0),
        "Provide the number of the parallel thread, 1: disable, 0: default.")(
        "vrs_scheme", po::value<int>(&vrs_scheme)->default_value(kMimic5),
        "Provide the scheme type, 0: mimc5, 1:sha256c")("hyrax_a1", "")(
        "hyrax_a1", po::value<int64_t>(&hyrax_a1_n)->default_value(0), "")(
        "hyrax_a2", po::value<int64_t>(&hyrax_a2_n)->default_value(0), "")(
        "hyrax_a3", po::value<int64_t>(&hyrax_a3_n)->default_value(0), "")(
        "gro09_51a", po::value<int64_t>(&gro09_51a_n)->default_value(0), "")(
        "gro09_51b", po::value<int64_t>(&gro09_51b_n)->default_value(0), "")(
        "gro09_52a", po::value<Gro09Param>(&gro09_52a)->multitoken(), "m*n, ex: 10*20")(
        "gro09_52b", po::value<Gro09Param>(&gro09_52b)->multitoken(), "m*n, ex: 10*20")(
        "gro09_53a", po::value<Gro09Param>(&gro09_53a)->multitoken(), "m*n, ex: 10*20")(
        "gro09_53b", po::value<Gro09Param>(&gro09_53b)->multitoken(), "m*n, ex: 10*20")(
        "gro09_43b", po::value<Gro09Param>(&gro09_43b)->multitoken(), "m*n, ex: 10*20")(
        "vrs_basic", po::value<int64_t>(&vrs_basic_n)->default_value(0), "")(
        "vrs_large", po::value<int64_t>(&vrs_large_n)->default_value(0), "")(
        "pod", po::value<PodParam>(&pod),"n*s, ex:100*1023");

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
  if (gro09_52a.m) {
    rets["groth09::sec52a"] = groth09::sec52a::TestRom(gro09_52a.m, gro09_52a.n);
  }
  if (gro09_52b.m) {
    rets["groth09::sec52b"] = groth09::sec52b::TestRom(gro09_52b.m, gro09_52b.n);
  }
  if (gro09_53a.m) {
    rets["groth09::sec53a"] = groth09::sec52a::TestRom(gro09_53a.m, gro09_53a.n);
  }
  if (gro09_53b.m) {
    rets["groth09::sec53b"] = groth09::sec52b::TestRom(gro09_53b.m, gro09_53b.n);
  }
  if (gro09_43b.m) {
    rets["groth09::sec43b"] = groth09::sec43b::TestRom(gro09_43b.m, gro09_43b.n);
  }

  if (vrs_scheme == VrsSchemeType::kMimic5) {
    if (vrs_basic_n) {
      rets["vrs::basic"] = vrs::TestBasic<vrs::Mimc5Scheme>(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::large"] = vrs::TestLarge<vrs::Mimc5Scheme>(vrs_basic_n);
    }
    if (pod.n) {
      rets["pod"] = pod::TestBasic<vrs::Mimc5Scheme>(pod.n, pod.s);
    }
  } else {
    if (vrs_basic_n) {
      rets["vrs::basic"] = vrs::TestBasic<vrs::Sha256cScheme>(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::large"] = vrs::TestLarge<vrs::Sha256cScheme>(vrs_basic_n);
    }
    if (pod.n) {
      rets["pod"] = pod::TestBasic<vrs::Sha256cScheme>(pod.n, pod.s);
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
