#pragma once

#include <libff/algebra/fields/field_utils.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libfqfft/evaluation_domain/evaluation_domain.hpp>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>
#include <libfqfft/polynomial_arithmetic/basis_change.hpp>
#include <libfqfft/polynomial_arithmetic/naive_evaluate.hpp>
#include <libfqfft/polynomial_arithmetic/xgcd.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include "ecc/ecc.h"

namespace snark {
namespace details {
inline void DisableLibffLog() {
  libff::inhibit_profiling_info = true;
  libff::inhibit_profiling_counters = true;
}

struct imembuf : public std::streambuf {
  imembuf(char const* array, size_t len) {
    char* p(const_cast<char*>(array));
    this->setg(p, p, p + len);
  }
  size_t get_count() { return gptr() - eback(); }
};

struct imemstream : virtual imembuf, std::istream {
  imemstream(char* array, size_t len)
      : imembuf(array, len), std::istream(static_cast<std::streambuf*>(this)) {}
};

struct omembuf : public std::streambuf {
  omembuf(char* array, size_t len) { this->setp(array, array + len); }
  size_t get_count() { return pptr() - pbase(); }
};

struct omemstream : virtual omembuf, std::ostream {
  omemstream(char* array, size_t len)
      : omembuf(array, len), std::ostream(static_cast<std::streambuf*>(this)) {}
};
}  // namespace details

typedef libsnark::r1cs_gg_ppzksnark_proof<
    libsnark::default_r1cs_gg_ppzksnark_pp>
    ZkProof;

typedef libff::Fr<libsnark::default_r1cs_gg_ppzksnark_pp> ZkFr;

typedef libsnark::r1cs_gg_ppzksnark_proving_key<
    libsnark::default_r1cs_gg_ppzksnark_pp>
    ZkPk;
typedef std::shared_ptr<ZkPk> ZkPkPtr;

typedef libsnark::r1cs_gg_ppzksnark_verification_key<
    libsnark::default_r1cs_gg_ppzksnark_pp>
    ZkVk;
typedef std::shared_ptr<ZkVk> ZkVkPtr;

inline void InitZkp(bool disable_log) {
  libsnark::default_r1cs_gg_ppzksnark_pp::init_public_params();
  if (disable_log) {
    details::DisableLibffLog();
  }
}

inline ZkFr ConvertToZkFr(Fr const& mcl_fr) {
  mpz_class m = mcl_fr.getMpz();
  return ZkFr(libff::bigint<ZkFr::num_limbs>(m.get_mpz_t()));
}

inline std::vector<ZkFr> ConvertToZkFr(std::vector<Fr> const& mcl_frs) {
  std::vector<ZkFr> zk_frs(mcl_frs.size());
  for (size_t i = 0; i < zk_frs.size(); ++i) {
    zk_frs[i] = ConvertToZkFr(mcl_frs[i]);
  }
  return zk_frs;
}

inline std::vector<ZkFr> ConvertToZkFr(std::vector<uint64_t> const& o) {
  std::vector<ZkFr> zk_frs(o.size());
  for (size_t i = 0; i < zk_frs.size(); ++i) {
    zk_frs[i] = ConvertToZkFr(Fr(o[i]));
  }
  return zk_frs;
}

// TODO: convert from zkfr
inline ZkPkPtr LoadZkPk(std::string const& file) {
  try {
    ZkPkPtr ret(new ZkPk());
    std::ifstream ifs;
    ifs.open(file, std::ifstream::in | std::ifstream::binary);
    ifs >> (*ret);
    return ret;
  } catch (std::exception& ex) {
    std::cerr << "Exception: " << ex.what() << "\n";
    return ZkPkPtr();
  }
}

inline ZkVkPtr LoadZkVk(std::string const& file) {
  try {
    ZkVkPtr ret(new ZkVk());
    std::ifstream ifs;
    ifs.open(file, std::ifstream::in | std::ifstream::binary);
    ifs >> (*ret);
    return ret;
  } catch (std::exception& ex) {
    std::cerr << "Exception: " << ex.what() << "\n";
    return ZkVkPtr();
  }
}

inline void ZkProofToBin(ZkProof const& proof, std::vector<uint8_t>& bin) {
  bin.resize(1024);
  details::omemstream out((char*)bin.data(), bin.size());
  out << proof;
  bin.resize(out.get_count());
}

inline void ZkProofFromBin(ZkProof& proof, std::vector<uint8_t> const& bin) {
  details::imemstream in((char*)bin.data(), bin.size());
  in >> proof;
  assert(in.get_count() == bin.size());
}
}  // namespace snark
