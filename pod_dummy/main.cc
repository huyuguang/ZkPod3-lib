#include <assert.h>
#include <cryptopp/blake2.h>
#include <cryptopp/osrng.h>
#include <cryptopp/randpool.h>

#include <iostream>

#include "bp/bp.h"
#include "circuit/test_all.h"
#include "clink/clink.h"
#include "cmd/cmd.h"
#include "debug/flags.h"
#include "ecc/ecc.h"
#include "groth09/groth09.h"
//#include "iop/iop.h"
#include "log/tick.h"
#include "misc/misc.h"
#include "public.h"

bool InitAll(std::string const& data_dir) {
  InitEcc();

  // auto ecc_pub_file = data_dir + "/" + "ecc_pub.bin";
  // if (!OpenOrCreateEccPub(ecc_pub_file)) {
  //  std::cerr << "Open or create ecc pub file " << ecc_pub_file << "
  //  failed\n"; return false;
  //}

  auto ecc_pds_file = data_dir + "/" + "pds_pub.bin";
  if (!pc::OpenOrCreatePdsPub(ecc_pds_file)) {
    std::cerr << "Open or create pds pub file " << ecc_pds_file << " failed\n";
    return false;
  }

  return true;
}

enum VrsSchemeType { kMimic5 = 0, kSha256c = 1, kPoseidon = 2 };
enum PolicyType { kOrdinary = 0, kSuccinct = 1 };

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

struct Param2Str {
  bool valid() const { return !s1.empty() && !s2.empty(); }
  std::string s1;
  std::string s2;
};

inline std::istream& operator>>(std::istream& in, Param2Str& t) {
  try {
    char sperator;
    in >> t.s1;
    in >> sperator;
    in >> t.s2;
  } catch (std::exception&) {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}

struct Param3Str {
  bool valid() const { return !s1.empty() && !s2.empty() && !s3.empty(); }
  std::string s1;
  std::string s2;
  std::string s3;
};

inline std::istream& operator>>(std::istream& in, Param3Str& t) {
  try {
    char sperator;
    in >> t.s1;
    in >> sperator;
    in >> t.s2;
    in >> sperator;
    in >> t.s3;
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

std::unique_ptr<tbb::task_scheduler_init> tbb_init;

int main(int argc, char** argv) {
  setlocale(LC_ALL, "");
  std::string data_dir;
  uint32_t thread_num;
  int vrs_scheme;
  int policy;
  bool hyrax_a1 = false;
  int64_t hyrax_a2_n = 0;
  int64_t hyrax_a3_n = 0;
  ParamIntPair hyrax_a4;
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
  bool circuit = false;
  bool opening = false;
  bool equality = false;
  bool equality2 = false;
  bool vcp_mnist = false;
  bool iop = false;
  ParamIntPair r1cs;
  int64_t pack_n = 0;
  ParamIntStr substrpack;
  ParamIntStr matchpack;
  Param2IntStr match_query;
  Param2IntStr substr_query;
  int64_t bp_p1_n = 0;
  int64_t bp_p2_n = 0;
  int64_t bp_p31_n = 0;
  int64_t vrs_cache_n = 0;
  int64_t pc_commitment_n = 0;
  int64_t multiexp_n = 0;
  int64_t mcl_n = 0;
  Param2Str vgg16_publish;
  Param2Str vgg16_infer;
  Param2Str vgg16_prove;
  bool vgg16_test = false;

  try {
    po::options_description options("command line options");
    options.add_options()("help,h", "Use -h or --help to list all arguments")(
        "data_dir,d", po::value<std::string>(&data_dir)->default_value("."),
        "Provide the data dir")(
        "thread_num", po::value<uint32_t>(&thread_num)->default_value(0),
        "Provide the number of the parallel thread, 1: disable, 0: default.")(
        "vrs_scheme", po::value<int>(&vrs_scheme)->default_value(kMimic5),
        "Provide the scheme type, 0: mimc5, 1:sha256c, 2:poseidon")(
        "policy", po::value<int>(&policy)->default_value(kOrdinary),
        "Provide the policy type, 0: ordinary(fast with large proof size), "
        "1:succinct(slow with small proof size)")(
        "vrs_cache", po::value<int64_t>(&vrs_cache_n), "")("hyrax_a1", "")(
        "hyrax_a2", po::value<int64_t>(&hyrax_a2_n), "")(
        "hyrax_a3", po::value<int64_t>(&hyrax_a3_n), "")(
        "hyrax_a4", po::value<ParamIntPair>(&hyrax_a4), "m*n, ex: 10*20")(
        "gro09_51a", po::value<int64_t>(&gro09_51a_n), "")(
        "gro09_51b", po::value<int64_t>(&gro09_51b_n), "")(
        "gro09_51c", po::value<int64_t>(&gro09_51c_n), "")(
        "gro09_52a", po::value<ParamIntPair>(&gro09_52a), "m*n, ex: 10*20")(
        "gro09_52b", po::value<ParamIntPair>(&gro09_52b), "m*n, ex: 10*20")(
        "gro09_53a", po::value<ParamIntPair>(&gro09_53a), "m*n, ex: 10*20")(
        "gro09_53b", po::value<ParamIntPair>(&gro09_53b), "m*n, ex: 10*20")(
        "gro09_43b", po::value<ParamIntPair>(&gro09_43b), "m*n, ex: 10*20")(
        "r1cs", po::value<ParamIntPair>(&r1cs), "m*n, ex: 10*20")(
        "vrs_basic", po::value<int64_t>(&vrs_basic_n), "")(
        "vrs_large", po::value<int64_t>(&vrs_large_n), "")(
        "pod", po::value<ParamIntPair>(&pod), "n*s, ex: 100*1023")(
        "matrix", po::value<ParamIntPair>(&matrix), "n*s, ex:10*20")(
        "equal_ip", po::value<ParamIntPair>(&equal_ip), "xn,yn, ex:10,20")(
        "overlap", "")("divide", "")("match", "")("substr", "")("circuit", "")(
        "pack", po::value<int64_t>(&pack_n), "n, ex: 31000")(
        "substrpack", po::value<ParamIntStr>(&substrpack), "n,key, ex: 10,abc")(
        "matchpack", po::value<ParamIntStr>(&substrpack), "n,key, ex: 10,abc")(
        "match_query", po::value<Param2IntStr>(&substr_query),
        "n,s,key, ex: 24,1000,abc")(
        "substr_query", po::value<Param2IntStr>(&substr_query),
        "n,s,key, ex: 24,1000,abc")("bp_p1", po::value<int64_t>(&bp_p1_n), "")(
        "bp_p2", po::value<int64_t>(&bp_p2_n), "")(
        "bp_p31", po::value<int64_t>(&bp_p31_n), "")(
        "pc_commitment", po::value<int64_t>(&pc_commitment_n), "")(
        "multiexp", po::value<int64_t>(&multiexp_n), "")(
        "disable_vrs_cache", "")("mcl", po::value<int64_t>(&mcl_n), "")(
        "opening", "")("equality", "")("equality2", "")("vcp_mnist", "")("iop", "")(
        "vgg16_publish", po::value<Param2Str>(&vgg16_publish),
        "\"para_path working_path\", ex: /vgg16/features /temp/vgg16")(
        "vgg16_infer", po::value<Param2Str>(&vgg16_infer),
        "\"test_image_path working_path\"")(
        "vgg16_prove", po::value<Param2Str>(&vgg16_prove),
        "test_image_path working_path")("vgg16_test", "");

    boost::program_options::variables_map vmap;

    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, options), vmap);
    boost::program_options::notify(vmap);

    if (vmap.count("help")) {
      std::cout << options << std::endl;
      return -1;
    }

    if (vrs_scheme != kMimic5 && vrs_scheme != kSha256c &&
        vrs_scheme != kPoseidon) {
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

    if (vmap.count("circuit")) {
      circuit = true;
    }

    if (vmap.count("disable_vrs_cache")) {
      debug::flags::disable_vrs_cache = true;
    }

    if (vmap.count("opening")) {
      opening = true;
    }

    if (vmap.count("equality")) {
      equality = true;
    }

    if (vmap.count("equality2")) {
      equality2 = true;
    }

    if (vmap.count("vcp_mnist")) {
      vcp_mnist = true;
    }

    if (vmap.count("iop")) {
      iop = true;
    }

    if (vmap.count("vgg16_test")) {
      vgg16_test = true;
    }
  } catch (std::exception& e) {
    std::cout << "Unknown parameters.\n"
              << e.what() << "\n"
              << "-h or --help to list all arguments.\n";
    return -1;
  }

  tbb_init = parallel::InitTbb((int)thread_num);

  if (!InitAll(data_dir)) {
    std::cerr << "Init failed\n";
    return -1;
  }

  std::map<std::string, bool> rets;

  if (vrs_cache_n) {
    if (vrs_scheme == VrsSchemeType::kMimic5) {
      using VrsCache = clink::VrsCache<clink::VrsMimc5Scheme>;
      rets["vrs_cache"] = VrsCache::CreateAndSave(data_dir, vrs_cache_n);
    } else if (vrs_scheme == VrsSchemeType::kSha256c) {
      using VrsCache = clink::VrsCache<clink::VrsSha256cScheme>;
      rets["vrs_cache"] = VrsCache::CreateAndSave(data_dir, vrs_cache_n);
    } else if (vrs_scheme == VrsSchemeType::kPoseidon) {
      using VrsCache = clink::VrsCache<clink::VrsPoseidonScheme>;
      rets["vrs_cache"] = VrsCache::CreateAndSave(data_dir, vrs_cache_n);
    } else {
      assert(false);
      rets["vrs_cache"] = false;
    }
  }

  if (mcl_n) {
    rets["mcl"] = TestMcl(mcl_n);
  }

  if (multiexp_n) {
    rets["multiexp"] = TestMultiexp(multiexp_n);
  }

  if (pc_commitment_n) {
    rets["pc_commiment"] = pc::TestPcCommitment(pc_commitment_n);
  }

  if (bp_p1_n) {
    rets["bp::p1"] = bp::p1::Test(bp_p1_n);
  }

  if (bp_p2_n) {
    rets["bp::p2"] = bp::p2::Test(bp_p2_n);
  }

  if (bp_p31_n) {
    rets["bp::p31"] = bp::p31::Test(bp_p31_n);
  }

  if (hyrax_a1) {
    rets["hyrax::A1"] = hyrax::A1::Test();
  }

  if (hyrax_a2_n) {
    rets["hyrax::A2"] = hyrax::A2::Test(hyrax_a2_n);
  }

  if (hyrax_a3_n) {
    rets["hyrax::A3"] = hyrax::A3::Test(hyrax_a3_n);
  }

  if (hyrax_a4.valid()) {
    rets["hyrax::A4"] = hyrax::A4::Test(hyrax_a4.x,hyrax_a4.y);
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
    if (policy == PolicyType::kOrdinary) {
      rets["groth09::sec53b(ordinary)"] =
          groth09::Sec53b<groth09::Sec51b>::Test(gro09_53b.x, gro09_53b.y);
    } else {
      rets["groth09::sec53b(succinct)"] =
          groth09::Sec53b<groth09::Sec51c>::Test(gro09_53b.x, gro09_53b.y);
    }
  }

  if (gro09_43b.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["groth09::sec43b(ordinary)"] =
          groth09::Sec43b<groth09::Sec53b<groth09::Sec51b>, hyrax::A2>::Test(
              gro09_43b.x, gro09_43b.y);
    } else {
      rets["groth09::sec43b(succinct)"] =
          groth09::Sec43b<groth09::Sec53b<groth09::Sec51c>, hyrax::A3>::Test(
              gro09_43b.x, gro09_43b.y);
    }
  }

  if (r1cs.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::r1cs(ordinary)"] =
          clink::ParallelR1cs<groth09::OrdinaryPolicy>::Test(r1cs.x, r1cs.y);
    } else {
      rets["clink::r1cs(succinct)"] =
          clink::ParallelR1cs<groth09::SuccinctPolicy>::Test(r1cs.x, r1cs.y);
    }
  }

  if (matrix.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::matrix_a2"] =
          clink::Matrix<hyrax::A2>::Test(matrix.x, matrix.y);
    } else {
      rets["clink::matrix_a3"] =
          clink::Matrix<hyrax::A3>::Test(matrix.x, matrix.y);
    }
  }

  if (equal_ip.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::equal_ip(ordinary)"] =
          clink::EqualIp<hyrax::A2>::Test(equal_ip.x, equal_ip.y);
    } else {
      rets["clink::equal_ip(succinct)"] =
          clink::EqualIp<hyrax::A3>::Test(equal_ip.x, equal_ip.y);
    }
  }

  if (overlap) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::overlap(ordinary)"] = clink::Overlap<hyrax::A2>::Test();
    } else {
      rets["clink::overlap(succinct)"] = clink::Overlap<hyrax::A3>::Test();
    }
  }

  if (divide) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::divide(ordinary)"] = clink::Divide<hyrax::A2>::Test();
    } else {
      rets["clink::divide(succinct)"] = clink::Divide<hyrax::A3>::Test();
    }
  }

  if (match) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::match(ordinary)"] =
          clink::Match<groth09::OrdinaryPolicy>::Test();
    } else {
      rets["clink::match(succinct)"] =
          clink::Match<groth09::SuccinctPolicy>::Test();
    }
  }

  if (substr) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::substr(ordinary)"] =
          clink::Substr<groth09::OrdinaryPolicy>::Test();
    } else {
      rets["clink::substr(succinct)"] =
          clink::Substr<groth09::SuccinctPolicy>::Test();
    }
  }

  if (circuit) {
    rets["circuit"] = circuit::Test();
  }

  if (pack_n) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::pack(ordinary)"] = clink::Pack<hyrax::A2>::Test(pack_n);
    } else {
      rets["clink::pack(succinct)"] = clink::Pack<hyrax::A3>::Test(pack_n);
    }
  }

  if (substrpack.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::substrpack(ordinary)"] =
          clink::SubstrPack<groth09::OrdinaryPolicy>::Test(substrpack.n,
                                                           substrpack.str);
    } else {
      rets["clink::substrpack(succinct)"] =
          clink::SubstrPack<groth09::SuccinctPolicy>::Test(substrpack.n,
                                                           substrpack.str);
    }
  }

  if (matchpack.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["clink::matchpack(ordinary)"] =
          clink::MatchPack<groth09::OrdinaryPolicy>::Test(matchpack.n,
                                                          matchpack.str);
    } else {
      rets["clink::matchpack(succinct)"] =
          clink::MatchPack<groth09::SuccinctPolicy>::Test(matchpack.n,
                                                          matchpack.str);
    }
  }

  if (match_query.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["cmd::match_query(ordinary)"] =
          cmd::MatchQuery<groth09::OrdinaryPolicy>::Test(
              match_query.n, match_query.s, match_query.str, data_dir);
    } else {
      rets["cmd::match_query(succinct)"] =
          cmd::MatchQuery<groth09::SuccinctPolicy>::Test(
              match_query.n, match_query.s, match_query.str, data_dir);
    }
  }

  if (substr_query.valid()) {
    if (policy == PolicyType::kOrdinary) {
      rets["cmd::substr_query(ordinary)"] =
          cmd::SubstrQuery<groth09::OrdinaryPolicy>::Test(
              substr_query.n, substr_query.s, substr_query.str, data_dir);
    } else {
      rets["cmd::substr_query(succinct)"] =
          cmd::SubstrQuery<groth09::SuccinctPolicy>::Test(
              substr_query.n, substr_query.s, substr_query.str, data_dir);
    }
  }

  if (vrs_basic_n) {
    if (vrs_scheme == VrsSchemeType::kMimic5) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(mimc5+ordinary)"] =
            clink::VrsBasic<clink::VrsMimc5Scheme,
                            groth09::OrdinaryPolicy>::Test(vrs_basic_n);
      } else {
        rets["vrs::basic(mimc5+succinct)"] =
            clink::VrsBasic<clink::VrsMimc5Scheme,
                            groth09::SuccinctPolicy>::Test(vrs_basic_n);
      }
    } else if (vrs_scheme == VrsSchemeType::kSha256c) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(sha256c+ordinary)"] =
            clink::VrsBasic<clink::VrsSha256cScheme,
                            groth09::OrdinaryPolicy>::Test(vrs_basic_n);
      } else {
        rets["vrs::basic(sha256c+succinct)"] =
            clink::VrsBasic<clink::VrsSha256cScheme,
                            groth09::SuccinctPolicy>::Test(vrs_basic_n);
      }
    } else if (vrs_scheme == VrsSchemeType::kPoseidon) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(poseidon+ordinary)"] =
            clink::VrsBasic<clink::VrsPoseidonScheme,
                            groth09::OrdinaryPolicy>::Test(vrs_basic_n);
      } else {
        rets["vrs::basic(poseidon+succinct)"] =
            clink::VrsBasic<clink::VrsPoseidonScheme,
                            groth09::SuccinctPolicy>::Test(vrs_basic_n);
      }
    } else {
      assert(false);
    }
  }

  if (vrs_large_n) {
    if (vrs_scheme == VrsSchemeType::kMimic5) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(mimc5+ordinary)"] =
            clink::VrsLarge<clink::VrsMimc5Scheme,
                            groth09::OrdinaryPolicy>::Test(vrs_large_n);
      } else {
        rets["vrs::basic(mimc5+succinct)"] =
            clink::VrsLarge<clink::VrsMimc5Scheme,
                            groth09::SuccinctPolicy>::Test(vrs_large_n);
      }
    } else if (vrs_scheme == VrsSchemeType::kSha256c) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(sha256c+ordinary)"] =
            clink::VrsLarge<clink::VrsSha256cScheme,
                            groth09::OrdinaryPolicy>::Test(vrs_large_n);
      } else {
        rets["vrs::basic(sha256c+succinct)"] =
            clink::VrsLarge<clink::VrsSha256cScheme,
                            groth09::SuccinctPolicy>::Test(vrs_large_n);
      }
    } else if (vrs_scheme == VrsSchemeType::kPoseidon) {
      if (policy == PolicyType::kOrdinary) {
        rets["vrs::basic(poseidon+ordinary)"] =
            clink::VrsLarge<clink::VrsPoseidonScheme,
                            groth09::OrdinaryPolicy>::Test(vrs_large_n);
      } else {
        rets["vrs::basic(poseidon+succinct)"] =
            clink::VrsLarge<clink::VrsPoseidonScheme,
                            groth09::SuccinctPolicy>::Test(vrs_large_n);
      }
    } else {
      assert(false);
    }
  }

  if (pod.valid()) {
    if (vrs_scheme == VrsSchemeType::kMimic5) {
      if (policy == PolicyType::kOrdinary) {
        rets["pod(mimc5+ordinary)"] =
            clink::Pod<clink::VrsMimc5Scheme, groth09::OrdinaryPolicy>::Test(
                pod.x, pod.y, data_dir);
      } else {
        rets["pod(mimc5+succinct)"] =
            clink::Pod<clink::VrsMimc5Scheme, groth09::SuccinctPolicy>::Test(
                pod.x, pod.y, data_dir);
      }
    } else if (vrs_scheme == VrsSchemeType::kSha256c) {
      if (policy == PolicyType::kOrdinary) {
        rets["pod(sha256c+ordinary)"] =
            clink::Pod<clink::VrsSha256cScheme, groth09::OrdinaryPolicy>::Test(
                pod.x, pod.y, data_dir);
      } else {
        rets["pod(sha256c+succinct)"] =
            clink::Pod<clink::VrsSha256cScheme, groth09::SuccinctPolicy>::Test(
                pod.x, pod.y, data_dir);
      }
    } else if (vrs_scheme == VrsSchemeType::kPoseidon) {
      if (policy == PolicyType::kOrdinary) {
        rets["pod(poseidon+ordinary)"] =
            clink::Pod<clink::VrsPoseidonScheme, groth09::OrdinaryPolicy>::Test(
                pod.x, pod.y, data_dir);
      } else {
        rets["pod(poseidon+succinct)"] =
            clink::Pod<clink::VrsPoseidonScheme, groth09::SuccinctPolicy>::Test(
                pod.x, pod.y, data_dir);
      }
    } else {
      assert(false);
    }
  }

  if (opening) {
    rets["opending"] = clink::Opening::Test();
  }

  if (equality) {
    rets["equality"] = clink::Equality::Test();
  }

  if (equality2) {
    rets["equality2"] = clink::Equality2::Test();
  }

  if (vcp_mnist) {
    if (policy == PolicyType::kOrdinary) {
      rets["vcp_mnist(ordinary)"] =
          clink::Mnist<groth09::OrdinaryPolicy>::Test();
    } else {
      rets["vcp_mnist(succinct)"] =
          clink::Mnist<groth09::SuccinctPolicy>::Test();
    }
  }

  if (iop) {
    std::cerr << "not support\n";
    //rets["iop"] = iop::Test();
  }

  if (vgg16_test) {
    rets["vgg16_test"] = clink::vgg16::Test();
  }

  if (vgg16_publish.valid()) {
    rets["vgg16_publish"] =
        clink::vgg16::Publish(vgg16_publish.s1, vgg16_publish.s2);
  }

  if (vgg16_infer.valid()) {
    //clink::vgg16::dbl::Test();
    rets["vgg16_infer"] = clink::vgg16::TestInfer(
        vgg16_infer.s1, vgg16_infer.s2);
  }

  if (vgg16_prove.valid()) {
    rets["vgg16_prove"] = clink::vgg16::TestProve(
        vgg16_prove.s1, vgg16_prove.s2);
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
