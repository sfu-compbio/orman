/// 786

#include "common.h"


uint64_t _zaman() {
	struct timeval t;
	gettimeofday(&t, 0);
	return (t.tv_sec * 1000000ll + t.tv_usec);
}

uint64_t zaman_start() {
	static uint64_t z = 0;
	if (!z) z = _zaman();
	return (_zaman() - z) / 1000000;
}

uint64_t zaman_last() {
	static uint64_t z = 0;
	uint64_t _z = _zaman() - z;
	z = _zaman();
	return _z / 1000000;
}
