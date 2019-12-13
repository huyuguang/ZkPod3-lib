#pragma once

#include <cryptopp/osrng.h>

// NOTE: NonblockingRng should enough for linux & windows
// thread_local CryptoPP::AutoSeededRandomPool rng;
inline thread_local CryptoPP::NonblockingRng tls_rng;
