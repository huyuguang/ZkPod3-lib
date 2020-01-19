#pragma once

#include <stdint.h>

#define SHA256_DIGESTSIZE 32

uint32_t inline ReadBE32(const uint8_t* ptr) {
	return ((uint32_t) *((ptr)+3)) | ((uint32_t) *((ptr)+2) << 8) |
		((uint32_t) *((ptr)+1) << 16) | ((uint32_t) *((ptr)+0) << 24);
}

void inline WriteBE32(uint8_t* ptr, uint32_t x) {
	*((ptr)+3) = (uint8_t)((x));
	*((ptr)+2) = (uint8_t)((x) >> 8);
	*((ptr)+1) = (uint8_t)((x) >> 16);
	*((ptr)+0) = (uint8_t)((x) >> 24);
}

void Sha256Compress(const uint8_t data[64], uint8_t hash[32]);

void Sha256Compress2(const uint32_t data[16], uint32_t hash[8]);