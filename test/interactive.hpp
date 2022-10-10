#define MORESTATS
#include "xoroshiro128pp.hpp"
#include <fstream>
#include <cstdlib>

#if defined(SIMD)
#include <function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
using RecSplit = sux::function::SIMDRecSplit<LEAF_SIZE, AT>;
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
using RecSplit = sux::function::GPURecSplit<LEAF_SIZE, AT>;
#else
#include <function/RecSplit.hpp>
#endif

using namespace sux::function;

void interactive(int argc, char** argv) {
	constexpr int LEAF_SIZE = 16;
	constexpr int BUCKET_SIZE = 100;
	RecSplit<LEAF_SIZE> recSplit;
	std::vector<std::string> strings;
	if (argc > 1) {
		int n = atoi(argv[1]);
		std::vector<hash128_t> keys;
		keys.reserve(n);
		for (uint64_t i = 0; i < n; i++) keys.push_back(hash128_t(next(), next()));
		recSplit = RecSplit<LEAF_SIZE>(keys, BUCKET_SIZE);
	} else {
		std::vector<hash128_t> keys;
		for(std::string key; getline(std::cin, key) && key != "";) {
			strings.push_back(key);
			keys.push_back(first_hash(key.c_str(), key.size()));
		}
		recSplit = RecSplit<LEAF_SIZE>(keys, BUCKET_SIZE);
	}
	for (std::string s; std::cin >> s;) {
		if (s == "/cout") {
			std::ofstream myfile;
			myfile.open("coutDump.data");

			myfile << recSplit;
			myfile.close();
		} else if (s == "/print") {
			std::vector<std::pair<size_t, std::string>> pairs;
			pairs.reserve(strings.size());
			for (const std::string& s : strings) pairs.emplace_back(recSplit(s), s);
			std::sort(pairs.begin(), pairs.end());
			for (const auto& p : pairs) std::cout << (p.second) << ": " << (p.first) << std::endl;
		} else {
			std::cout << recSplit(s) << std::endl;
		}
	}
}
