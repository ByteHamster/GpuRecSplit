#pragma once

#include <immintrin.h>
#include <vectorclass.h>

#if INSTRSET >= 9
#define SIMDRS_512_BIT
constexpr bool NO_AVX = false;

#ifdef __AVX512VPOPCNTDQ__
#define SIMDRS_512_BIT_POPCNT
#endif

#elif INSTRSET == 8
constexpr bool NO_AVX = false;
#else
constexpr bool NO_AVX = true;
#pragma message("SIMDRecSplit was compiled without AVX512 and AVX2 support => suboptimal performance")
#endif

#ifdef SIMDRS_512_BIT
using FullVecUi = Vec16ui;
using FullVecUq = Vec8uq;
using FullVecIb = Vec16ib;
using FullVecQ = Vec8q;
using FullVecC = Vec64c;
constexpr uint32_t FULL_VEC_32_COUNT = 16;
#else
using FullVecUi = Vec8ui;
using FullVecUq = Vec4uq;
using FullVecIb = Vec8ib;
using FullVecQ = Vec4q;
using FullVecC = Vec32c;
constexpr uint32_t FULL_VEC_32_COUNT = 8;
#endif
constexpr uint32_t FULL_VEC_64_COUNT = FULL_VEC_32_COUNT / 2;

FullVecUq shift(FullVecUq lhs, FullVecUq rhs) {
#ifdef SIMDRS_512_BIT
	return _mm512_sllv_epi64(lhs, rhs);
#else
	return _mm256_sllv_epi64(lhs, rhs);
#endif
}
