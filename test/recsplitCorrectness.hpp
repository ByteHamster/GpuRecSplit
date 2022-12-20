#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <util/Hash128.h>

template<class RS>
bool testCorrectness(RS &rs, const std::vector<bez::function::hash128_t> &keys) {
	constexpr std::size_t bits = sizeof(std::uint64_t) * 8;
	std::vector<std::uint64_t> results((keys.size() + bits - 1) / bits);
	for (const auto &key : keys) {
		std::size_t result = rs(key);
		std::size_t int_num = result / bits;
		results[int_num] |= 1ULL << (result - int_num * bits);
	}

	for (std::size_t i = 0; i < results.size() - 1; ++i) {
		if (results[i] != std::uint64_t(-1)) {
			std::cout << "Missing result between " << (i * bits) << " (inclusive) and " << ((i + 1) * bits)
				<< " (exclusive): " << results[i] << "\n";
			return false;
		}
	}
	std::size_t m = keys.size() % bits;
	if (m > 0 && results[results.size() - 1] != (1ULL << m) - 1) {
		std::cout << "Missing result between " << ((results.size() - 1) * bits) << " (inclusive) and "
			<< keys.size() << " (exclusive): " << results[results.size() - 1] << "\n";
		return false;
	}
	std::printf("Result is a correct MPHF!\n");
	return true;
}