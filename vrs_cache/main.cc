#include "ecc.h"
#include "public.h"
#include "vrs/vrs.h"

bool InitAll(std::string const& data_dir) {
  InitEcc();

  auto ecc_pub_file = data_dir + "/" + "ecc_pub.bin";
  if (!OpenOrCreateEccPub(ecc_pub_file)) {
    std::cerr << "Open or create ecc pub file " << ecc_pub_file << " failed\n";
    return false;
  }

  auto ecc_pds_file = data_dir + "/" + "pds_pub.bin";
  if (!OpenOrCreatePdsPub(ecc_pds_file)) {
    std::cerr << "Open or create pds pub file " << ecc_pub_file << " failed\n";
    return false;
  }

  return true;
}

int main(int argc, char** argv) {
  setlocale(LC_ALL, "");
  std::string data_dir;
  uint64_t count;
  uint32_t thread_num = 0;

  try {
    po::options_description options("command line options");
    options.add_options()("help,h", "Use -h or --help to list all arguments")(
        "data_dir,d", po::value<std::string>(&data_dir)->default_value("."),
        "Provide the data dir")(
        "count,c", po::value<uint64_t>(&count)->default_value(2),
        "Provide the count, must >1, should be (n+1)*s or multiple 32k")(
        "thread_num", po::value<uint32_t>(&thread_num)->default_value(0),
        "Provide the number of the parallel thread, 1: disable, 0: default.");

    boost::program_options::variables_map vmap;

    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, options), vmap);
    boost::program_options::notify(vmap);

    if (vmap.count("help")) {
      std::cout << options << std::endl;
      return -1;
    }

    if (count <= 1) {
      std::cout << options << std::endl;
      return -1;
    }
  } catch (std::exception& e) {
    std::cout << "Unknown parameters.\n"
              << e.what() << "\n"
              << "-h or --help to list all arguments.\n";
    return -1;
  }

#ifdef USE_TBB
  int tbb_thread_num =
      thread_num ? (int)thread_num : tbb::task_scheduler_init::automatic;
  tbb::task_scheduler_init init(tbb_thread_num);
  std::cout << "use tbb\n";
#else
  setenv("options:thread_num", std::to_string(thread_num).c_str(), true);
#endif

  if (!InitAll(data_dir)) {
    std::cerr << "Init failed\n";
    return -1;
  }

  std::string cache_dir = data_dir + "/vrs_cache";
  fs::create_directories(cache_dir);
  if (!fs::is_directory(cache_dir)) {
    std::cerr << "create directory failed: " << cache_dir << "\n";
    return -1;
  }

  auto cache = vrs::CreateCache(count);

  std::string cache_file;
  bool ret = vrs::SaveCache(cache_dir, cache, cache_file);

  if (ret) {
    std::cout << "Success: cache_file: " << cache_file << "\n";
  } else {
    assert(false);
    std::cout << "Save Failed: " << cache_file << "\n";
  }

  return ret ? 0 : -1;
}
