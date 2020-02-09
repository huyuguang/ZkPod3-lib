#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>

#include <iostream>

#include "bp/bp.h"
#include "cmd/cmd.h"
#include "ecc/ecc.h"
#include "groth09/groth09.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "clink/clink.h"
#include "public.h"

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

struct ParamIntStr {
  bool valid() const { return n && !str.empty(); }
  int64_t n = 0;
  std::string str;
};

inline std::istream& operator>>(std::istream& in, ParamIntStr& t) {
  try {
    char sperator;
    in >> t.n;
    in >> sperator;
    in >> t.str;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

struct Param2IntStr {
  bool valid() const { return n && s && !str.empty(); }
  int64_t n = 0;
  int64_t s = 0;
  std::string str;
};

inline std::istream& operator>>(std::istream& in, Param2IntStr& t) {
  try {
    char sperator;
    in >> t.n;
    in >> sperator;
    in >> t.s;
    in >> sperator;
    in >> t.str;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

int main(int argc, char** argv) {
  setlocale(LC_ALL, "");
  std::string data_dir;
  uint32_t thread_num;
  int vrs_scheme;
  bool hyrax_a1 = false;
  int64_t hyrax_a2_n = 0;
  int64_t hyrax_a3_n = 0;
  int64_t gro09_51a_n = 0;
  int64_t gro09_51b_n = 0;
  int64_t gro09_51c_n = 0;
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
  ParamIntStr substrpack;
  ParamIntStr matchpack;
  Param2IntStr match_query;
  Param2IntStr substr_query;
  int64_t bp_p1_n = 0;
  int64_t bp_p2_n = 0;
  int64_t bp_p3_n = 0;

  try {
    po::options_description options("command line options");
    options.add_options()("help,h", "Use -h or --help to list all arguments")(
        "data_dir,d", po::value<std::string>(&data_dir)->default_value("."),
        "Provide the data dir")(
        "thread_num", po::value<uint32_t>(&thread_num)->default_value(0),
        "Provide the number of the parallel thread, 1: disable, 0: default.")(
        "vrs_scheme", po::value<int>(&vrs_scheme)->default_value(kMimic5),
        "Provide the scheme type, 0: mimc5, 1:sha256c")("hyrax_a1", "")(
        "hyrax_a2", po::value<int64_t>(&hyrax_a2_n), "")(
        "hyrax_a3", po::value<int64_t>(&hyrax_a3_n), "")(
        "gro09_51a", po::value<int64_t>(&gro09_51a_n), "")(
        "gro09_51b", po::value<int64_t>(&gro09_51b_n), "")(
        "gro09_51c", po::value<int64_t>(&gro09_51c_n), "")(
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
        "overlap", "")("divide", "")("match", "")("substr", "")("pack", "")(
        "substrpack", po::value<ParamIntStr>(&substrpack), "n,key, ex: 10,abc")(
        "matchpack", po::value<ParamIntStr>(&substrpack), "n,key, ex: 10,abc")(
        "match_query", po::value<Param2IntStr>(&substr_query),
        "n,s,key, ex: 24,1000,abc")(
        "substr_query", po::value<Param2IntStr>(&substr_query),
        "n,s,key, ex: 24,1000,abc")("bp_p1", po::value<int64_t>(&bp_p1_n), "")(
        "bp_p2", po::value<int64_t>(&bp_p2_n), "")(
        "bp_p3", po::value<int64_t>(&bp_p3_n), "");

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

    if (vmap.count("hyrax_a1")) {
      hyrax_a1 = true;
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

  if (hyrax_a1) {
    rets["hyrax::A1"] = hyrax::A1::Test();
  }
  if (hyrax_a2_n) {
    rets["hyrax::A2"] = hyrax::A2::Test(hyrax_a2_n);
  }
  if (hyrax_a3_n) {
    rets["hyrax::A3"] = hyrax::A3::Test(hyrax_a3_n);
  }
  if (gro09_51a_n) {
    rets["groth09::sec51a"] = groth09::Sec51a::Test(gro09_51a_n);
  }
  if (gro09_51b_n) {
    rets["groth09::sec51b"] = groth09::Sec51b::Test(gro09_51b_n);
  }
  if (gro09_51c_n) {
    rets["groth09::sec51c"] = groth09::Sec51c::Test(gro09_51c_n);
  }
  if (gro09_52a.valid()) {
    rets["groth09::sec52a"] = groth09::Sec52a::Test(gro09_52a.x, gro09_52a.y);
  }
  if (gro09_52b.valid()) {
    rets["groth09::sec52b"] = groth09::Sec52b::Test(gro09_52b.x, gro09_52b.y);
  }
  if (gro09_53a.valid()) {
    rets["groth09::sec53a"] = groth09::Sec53a::Test(gro09_53a.x, gro09_53a.y);
  }
  if (gro09_53b.valid()) {
    rets["groth09::sec53b_1"] =
        groth09::Sec53b<groth09::Sec51b>::Test(gro09_53b.x, gro09_53b.y);
    rets["groth09::sec53b_2"] =
        groth09::Sec53b<groth09::Sec51c>::Test(gro09_53b.x, gro09_53b.y);
  }
  if (gro09_43b.valid()) {
    rets["groth09::sec43b_1"] =
        groth09::Sec43b<groth09::Sec53b<groth09::Sec51b>, hyrax::A2>::Test(
            gro09_43b.x, gro09_43b.y);
    rets["groth09::sec43b_2"] =
        groth09::Sec43b<groth09::Sec53b<groth09::Sec51c>, hyrax::A3>::Test(
            gro09_43b.x, gro09_43b.y);
  }
  if (matrix.valid()) {
    rets["clink::matrix_a2"] =
        clink::Matrix<hyrax::A2>::Test(matrix.x, matrix.y);
    rets["clink::matrix_a3"] =
        clink::Matrix<hyrax::A3>::Test(matrix.x, matrix.y);
  }
  if (equal_ip.valid()) {
    rets["clink::equal_ip_1"] =
        clink::EqualIp<hyrax::A2>::Test(equal_ip.x, equal_ip.y);
    rets["clink::equal_ip_2"] =
        clink::EqualIp<hyrax::A3>::Test(equal_ip.x, equal_ip.y);
  }
  if (overlap) {
    rets["clink::overlap_1"] = clink::Overlap<hyrax::A2>::Test();
    rets["clink::overlap_2"] = clink::Overlap<hyrax::A3>::Test();
  }
  if (divide) {
    rets["clink::divide_1"] = clink::Divide<hyrax::A2>::Test();
    rets["clink::divide_2"] = clink::Divide<hyrax::A3>::Test();
  }
  if (match) {
    rets["clink::match_1"] =
        clink::Match<groth09::OrdinaryPolicy>::Test();
    rets["clink::match_2"] =
        clink::Match<groth09::SuccinctPolicy>::Test();
  }
  if (substr) {
    rets["clink::substr_1"] =
        clink::Substr<groth09::OrdinaryPolicy>::Test();
    rets["clink::substr_2"] =
        clink::Substr<groth09::SuccinctPolicy>::Test();
  }
  if (pack) {
    rets["clink::pack_1"] = clink::Pack<hyrax::A2>::Test();
    rets["clink::pack_2"] = clink::Pack<hyrax::A3>::Test();
  }
  if (substrpack.valid()) {
    rets["clink::substrpack_1"] =
        clink::SubstrPack<groth09::OrdinaryPolicy>::Test(substrpack.n,
                                                            substrpack.str);
    rets["clink::substrpack_2"] =
        clink::SubstrPack<groth09::SuccinctPolicy>::Test(substrpack.n,
                                                            substrpack.str);
  }
  if (matchpack.valid()) {
    rets["clink::matchpack_1"] =
        clink::MatchPack<groth09::OrdinaryPolicy>::Test(matchpack.n,
                                                           matchpack.str);
    rets["clink::matchpack_2"] =
        clink::MatchPack<groth09::SuccinctPolicy>::Test(matchpack.n,
                                                           matchpack.str);
  }
  if (match_query.valid()) {
    rets["cmd::match_query_1"] = cmd::MatchQuery<groth09::OrdinaryPolicy>::Test(
        match_query.n, match_query.s, match_query.str);
    rets["cmd::match_query_2"] = cmd::MatchQuery<groth09::SuccinctPolicy>::Test(
        match_query.n, match_query.s, match_query.str);
  }
  if (substr_query.valid()) {
    rets["cmd::substr_query_1"] =
        cmd::SubstrQuery<groth09::OrdinaryPolicy>::Test(
            substr_query.n, substr_query.s, substr_query.str);
    rets["cmd::substr_query_2"] =
        cmd::SubstrQuery<groth09::SuccinctPolicy>::Test(
            substr_query.n, substr_query.s, substr_query.str);
  }

  if (vrs_scheme == VrsSchemeType::kMimic5) {
    if (vrs_basic_n) {
      rets["vrs::basic_1"] =
          clink::VrsBasic<clink::VrsMimc5Scheme,
                             groth09::OrdinaryPolicy>::Test(vrs_basic_n);
      rets["vrs::basic_2"] =
          clink::VrsBasic<clink::VrsMimc5Scheme,
                             groth09::SuccinctPolicy>::Test(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::basic_1"] =
          clink::VrsLarge<clink::VrsMimc5Scheme,
                             groth09::OrdinaryPolicy>::Test(vrs_large_n);
      rets["vrs::basic_2"] =
          clink::VrsLarge<clink::VrsMimc5Scheme,
                             groth09::SuccinctPolicy>::Test(vrs_large_n);
    }
    if (pod.valid()) {
      rets["pod_1"] =
          clink::Pod<clink::VrsMimc5Scheme,
                        groth09::OrdinaryPolicy>::Test(pod.x, pod.y);
      rets["pod_2"] =
          clink::Pod<clink::VrsMimc5Scheme,
                        groth09::SuccinctPolicy>::Test(pod.x, pod.y);
    }
  } else {
    if (vrs_basic_n) {
      rets["vrs::basic_1"] =
          clink::VrsBasic<clink::VrsSha256cScheme,
                             groth09::OrdinaryPolicy>::Test(vrs_basic_n);
      rets["vrs::basic_2"] =
          clink::VrsBasic<clink::VrsSha256cScheme,
                             groth09::SuccinctPolicy>::Test(vrs_basic_n);
    }
    if (vrs_large_n) {
      rets["vrs::basic_1"] =
          clink::VrsLarge<clink::VrsSha256cScheme,
                             groth09::OrdinaryPolicy>::Test(vrs_large_n);
      rets["vrs::basic_2"] =
          clink::VrsLarge<clink::VrsSha256cScheme,
                             groth09::SuccinctPolicy>::Test(vrs_large_n);
    }

    if (pod.valid()) {
      rets["pod_1"] =
          clink::Pod<clink::VrsSha256cScheme,
                        groth09::OrdinaryPolicy>::Test(pod.x, pod.y);
      rets["pod_2"] =
          clink::Pod<clink::VrsSha256cScheme,
                        groth09::SuccinctPolicy>::Test(pod.x, pod.y);
    }
  }

  if (bp_p1_n) {
    rets["bp::p1"] = bp::p1::Test(bp_p1_n);
  }

  if (bp_p2_n) {
    rets["bp::p2"] = bp::p2::Test(bp_p2_n);
  }

  if (bp_p3_n) {
    rets["bp::p3"] = bp::p3::Test(bp_p3_n);
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
