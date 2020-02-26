#include <depends/libff/libff/common/utils.cpp>

#ifdef _WIN32 // for linux we link the libmcl.a
#include <depends/mcl/src/fp.cpp>
#endif