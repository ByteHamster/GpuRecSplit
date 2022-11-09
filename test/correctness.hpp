#pragma once

#include <XorShift64.h>
#include "recsplitCorrectness.hpp"
#include <algorithm>
#include <cstdio>
#include <iostream>
#if defined(SIMD)
#include <function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE, bez::util::AllocType AT = bez::util::AllocType::MALLOC>
using TestRecSplit = bez::function::SIMDRecSplit<LEAF_SIZE, AT>;
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, bez::util::AllocType AT = bez::util::AllocType::MALLOC>
using TestRecSplit = bez::function::GPURecSplit<LEAF_SIZE, AT>;
#else
#include <function/RecSplit.hpp>
template<size_t LEAF_SIZE, bez::util::AllocType AT = bez::util::AllocType::MALLOC>
using TestRecSplit = bez::function::RecSplit<LEAF_SIZE, AT>;
#endif

using namespace std;

static constexpr size_t sizes[] = { 1, 10, 100, 123, 10000, 100000, 1000000 };// , 10000000, 12345678 };
static constexpr size_t bucket_sizes[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 75, 100, 231, 483, 1009, 1300, 2000 };
static constexpr int MIN_TEST_LEAF_SIZE = 2;
static constexpr int MAX_TEST_LEAF_SIZE = 24;

template<int FROM_LEAF, int TO_LEAF>
bool test() {
	bool correct = true;
	std::vector<bez::function::hash128_t> keys;
	for (const size_t bucket_size : bucket_sizes) {
		for (size_t size : sizes) {
			if (FROM_LEAF > 20)
				size = min((unsigned long long)size, 2000ULL);
			else if (FROM_LEAF > 16)
				size = min((unsigned long long)size, 1'000'000ULL);

            util::XorShift64 prng(0x5603141978c51071);
			for (uint64_t i = 0; i < size; i++) keys.push_back(bez::function::hash128_t(prng(), prng()));

#if (defined(SIMD) || defined(GPU))
			int num_threads = std::thread::hardware_concurrency();
			num_threads = num_threads == 0 ? 1 : num_threads;
			TestRecSplit<FROM_LEAF> rs(keys, bucket_size, num_threads);
#else
			TestRecSplit<FROM_LEAF> rs(keys, bucket_size);
#endif

			std::cout << "l = " << FROM_LEAF << ", b = " << bucket_size << ", n = " << size << ": ";
			correct &= testCorrectness(rs, keys);
			keys.clear();
		}
	}
	keys.shrink_to_fit(); // Free memory before recursive call
	if (correct)
		printf("All tests for leaf size %d are correct!\n", FROM_LEAF);
	else
		printf("There were errors for leaf size %d!\n", FROM_LEAF);

	if constexpr (FROM_LEAF < TO_LEAF)
		correct &= test<FROM_LEAF + 1, TO_LEAF>();
	return correct;
}

void testAll() {
	if (test<MIN_TEST_LEAF_SIZE, MAX_TEST_LEAF_SIZE>())
		printf("\nAll tests are correct!\n");
	else
		printf("\nThere were errors!\n");
}
