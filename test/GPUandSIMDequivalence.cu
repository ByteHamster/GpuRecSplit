#pragma once

#include "xoroshiro128pp.hpp"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <function/SIMDRecSplit.hpp>
#include <function/GPURecSplit.cuh>

using namespace std;
using namespace bez::function;

static constexpr size_t sizes[] = { 1, 10, 100, 123, 10000, 100000, 1000000, 10000000, 12345678 };
static constexpr size_t bucket_sizes[] = { 2, 5, 10, 25, 75, 100, 483, 1300, 2000 };
static constexpr int MIN_TEST_LEAF_SIZE = 2;
static constexpr int MAX_TEST_LEAF_SIZE = 24;

template<class RS1, class RS2>
bool testEquivalence(RS1 &rs1, RS2 &rs2, const std::vector<hash128_t> &keys) {
	for (const auto &key : keys) {
		std::size_t result1 = rs1(key);
		std::size_t result2 = rs2(key);
		if (result1 != result2) {
			std::cout << "There are differences for key (" << key.first << "," << key.second << "): " << result1 << " != " << result2 << std::endl;
			return false;
		}
	}
	std::cout << "Both RecSplit instances are equivalent!\n";
	return true;
}

template<int FROM_LEAF, int TO_LEAF>
bool test() {
	bool equivalent = true;
	std::vector<hash128_t> keys;
	for (const size_t bucket_size : bucket_sizes) {
		for (size_t size : sizes) {
			if (FROM_LEAF > 20)
				size = min((unsigned long long)size, 2000ULL);
			else if (FROM_LEAF > 16)
				size = min((unsigned long long)size, 1'000'000ULL);

			for (uint64_t i = 0; i < size; i++) keys.push_back(hash128_t(next(), next()));

			GPURecSplit<FROM_LEAF> gpurs(keys, bucket_size);
			SIMDRecSplit<FROM_LEAF> simdrs(keys, bucket_size);
			std::cout << "l = " << FROM_LEAF << ", b = " << bucket_size << ", n = " << size << ": ";
			equivalent &= testEquivalence(gpurs, simdrs, keys);
			keys.clear();
		}
	}
	keys.shrink_to_fit(); // Free memory before recursive call
	if (equivalent)
		printf("All tests for leaf size %d are correct!\n", FROM_LEAF);
	else
		printf("There were errors for leaf size %d!\n", FROM_LEAF);

	if constexpr (FROM_LEAF < TO_LEAF)
		equivalent &= test<FROM_LEAF + 1, TO_LEAF>();
	return equivalent;
}

int main() {
	if (test<MIN_TEST_LEAF_SIZE, MAX_TEST_LEAF_SIZE>())
		printf("\nAll tests are correct!\n");
	else
		printf("\nThere were errors!\n");
}
