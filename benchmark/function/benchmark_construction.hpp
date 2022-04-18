#pragma once

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <test/xoroshiro128pp.hpp>
#include <src/util/csv_printer.hpp>

#if defined(SIMD)
#include <src/function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE, sux::util::AllocType AT, bool USE_BIJECTIONS_ROTATE>
using RecSplit = sux::function::SIMDRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE>;
const std::string FILE_NAME = "simdBenchmark";
#elif defined(GPU)
#include <src/function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, sux::util::AllocType AT, bool USE_BIJECTIONS_ROTATE>
using RecSplit = sux::function::GPURecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE>;
const std::string FILE_NAME = "gpuBenchmark";
#else
#include <src/function/RecSplit.hpp>
const std::string FILE_NAME = "stdBenchmark";
#endif

using namespace sux::function;

static constexpr size_t sizes[] = { 1'000'000, 1'000'000'000 };
static constexpr size_t bucket_sizes[] = { 5, 50, 500, 2000 };
static constexpr int MIN_TEST_LEAF_SIZE = 5;
static constexpr int MAX_TEST_LEAF_SIZE = 17;
static constexpr int ITERATIONS = 5;
static constexpr int QUERIES = 1'000'000;

template<class RS>
double benchmarkQueries(RS &rs) {
	uint64_t h = 0;
	auto begin = chrono::high_resolution_clock::now();
	for (int i = 0; i < QUERIES; i++) h ^= rs(hash128_t(next(), next() ^ h));
	auto end = chrono::high_resolution_clock::now();
	const uint64_t elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
	const volatile uint64_t unused = h;
	(void)unused;
	return (double)elapsed / QUERIES;
}

template<int FROM_LEAF, int TO_LEAF, bool USE_BIJECTIONS_ROTATE>
void construct(std::ofstream &out) {
	std::vector<hash128_t> keys;
	for (const size_t bucket_size : bucket_sizes) {
		for (size_t size : sizes) {
			if (FROM_LEAF > 10 && size > 1'000'000)
				continue;
			for (int iteration = 0; iteration < ITERATIONS; ++iteration) {
				for (uint64_t i = 0; i < size; i++) keys.push_back(hash128_t(next(), next()));

				auto begin = chrono::high_resolution_clock::now();
				RecSplit<FROM_LEAF, sux::util::AllocType::MALLOC, USE_BIJECTIONS_ROTATE> rs(keys, bucket_size);
				auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count();
				csvPrint(out, FROM_LEAF, bucket_size, size, iteration, elapsed / (double)size, rs.getBitsPerKey(), benchmarkQueries(rs));
				keys.clear();
			}
		}
	}
	keys.shrink_to_fit(); // Free memory before recursive call

    if constexpr (FROM_LEAF < TO_LEAF)
        construct<FROM_LEAF + 1, TO_LEAF, USE_BIJECTIONS_ROTATE>(out);
}

void constructAll() {
	uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	std::ofstream out(FILE_NAME + "NoRot.csv");
	csvPrint(out, "leafSize", "bucketSize", "size", "iteration", "timingPerElement[ns]", "bitsPerKey", "queryTime[ns]");
    construct<MIN_TEST_LEAF_SIZE, MAX_TEST_LEAF_SIZE, false>(out);
	out.close();

	s[0] = s0; // Reset seeds to ensure both variants receive the same elements for construction and queries
	s[1] = s1;
	std::ofstream outRot(FILE_NAME + "Rot.csv");
	csvPrint(outRot, "leafSize", "bucketSize", "size", "iteration", "timingPerElement[ns]", "bitsPerKey", "queryTime[ns]");
	construct<MIN_TEST_LEAF_SIZE, MAX_TEST_LEAF_SIZE, true>(outRot);
	outRot.close();
}
