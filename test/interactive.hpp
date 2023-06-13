#define MORESTATS
#include <XorShift64.h>
#include <fstream>
#include <cstdlib>

#if defined(SIMD)
#include <function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
using RecSplit = bez::function::SIMDRecSplit<LEAF_SIZE, AT>;
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
using RecSplit = bez::function::GPURecSplit<LEAF_SIZE, AT>;
#else
#include <function/RecSplitRotate.hpp>
template<size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC>
using RecSplit = bez::function::recsplit_rotate::RecSplit<LEAF_SIZE, AT>;
#endif

void interactive(int argc, char** argv) {
	constexpr int LEAF_SIZE = 16;
	constexpr int BUCKET_SIZE = 100;
	std::vector<std::string> strings;
	std::vector<bez::function::hash128_t> keys;
	if (argc > 1) {
		int n = atoi(argv[1]);
		keys.reserve(n);
        util::XorShift64 prng(0x5603141978c51071);
		for (uint64_t i = 0; i < n; i++) keys.push_back(bez::function::hash128_t(prng(), prng()));
	} else {
		for (std::string key; getline(std::cin, key) && key != "";) {
			strings.push_back(key);
			keys.push_back(bez::function::first_hash(key.c_str(), key.size()));
		}
	}
	int num_threads = std::thread::hardware_concurrency();
	num_threads = num_threads == 0 ? 1 : num_threads;
	auto recSplit = RecSplit<LEAF_SIZE>(keys, BUCKET_SIZE, num_threads);
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
