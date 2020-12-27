#include <depends/libff/libff/common/utils.cpp>
#include <libsnark/common/data_structures/integer_permutation.cpp>
#include <libsnark/common/routing_algorithms/as_waksman_routing_algorithm.cpp>

#ifdef _WIN32  // for linux we link the libmcl.a
#include <depends/mcl/src/fp.cpp>
#endif