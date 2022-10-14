#pragma once

#include "xoroshiro128pp.hpp"
#include "recsplitCorrectness.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>

#ifndef STATS
#define STATS
#endif

#if defined(SIMD)
#include <function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE, bez::util::AllocType AT = bez::util::AllocType::MALLOC>
using RecSplit = bez::function::SIMDRecSplit<LEAF_SIZE, AT>;
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, bez::util::AllocType AT = bez::util::AllocType::MALLOC>
using RecSplit = bez::function::GPURecSplit<LEAF_SIZE, AT>;
#else
#include <function/RecSplit.hpp>
#endif

#define SAMPLES (11)

#ifndef LEAF
#define LEAF 5
#endif

using namespace std;
using namespace bez::function;

template<class RS>
void benchmark(RS &rs, const uint64_t n) {
	printf("Benchmarking...\n");

	uint64_t sample[SAMPLES];
	uint64_t h = 0;

	for (int k = SAMPLES; k-- != 0;) {
		s[0] = 0x5603141978c51071;
		s[1] = 0x3bbddc01ebdf4b72;
		auto begin = chrono::high_resolution_clock::now();
		for (uint64_t i = 0; i < n; i++) h ^= rs(hash128_t(next(), next() ^ h));
		auto end = chrono::high_resolution_clock::now();
		const uint64_t elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
		sample[k] = elapsed;
		printf("Elapsed: %.3fs; %.3f ns/key\n", elapsed * 1E-9, elapsed / (double)n);
	}

	const volatile uint64_t unused = h;
	(void)unused;
	sort(sample, sample + SAMPLES);
	printf("\nMedian: %.3fs; %.3f ns/key\n", sample[SAMPLES / 2] * 1E-9, sample[SAMPLES / 2] / (double)n);
}

int build(int argc, char **argv) {
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <n> <bucket size> <mphf>\n", argv[0]);
		return 1;
	}

	const uint64_t n = strtoll(argv[1], NULL, 0);
	const size_t bucket_size = strtoll(argv[2], NULL, 0);
	std::vector<hash128_t> keys;
	for (uint64_t i = 0; i < n; i++) keys.push_back(hash128_t(next(), next()));

	printf("Building...\n");

#ifdef ALLOC_TYPE
#define ALLOC_TYPE_APPEND ,ALLOC_TYPE
#else
#define ALLOC_TYPE_APPEND
#endif

#if defined(SIMD)
	int num_threads = std::thread::hardware_concurrency();
	num_threads = num_threads == 0 ? 1 : num_threads;
	auto begin = chrono::high_resolution_clock::now();
	RecSplit<LEAF ALLOC_TYPE_APPEND> rs(keys, bucket_size, num_threads);
#else
	auto begin = chrono::high_resolution_clock::now();
	RecSplit<LEAF ALLOC_TYPE_APPEND> rs(keys, bucket_size);
#endif

	auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count();
	printf("Construction time: %.3f s, %.0f ns/key\n", elapsed * 1E-9, elapsed / (double)n);

#ifdef DUMP
	fstream fs;
	fs.exceptions(fstream::failbit | fstream::badbit);
	fs.open(argv[3], fstream::out | fstream::binary | fstream::trunc);
	fs << rs;
	fs.close();
#endif

	testCorrectness(rs, keys);
    benchmark(rs, n);

	return 0;
}
