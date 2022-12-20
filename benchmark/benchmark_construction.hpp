#pragma once

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <XorShift64.h>
#include <tlx/cmdline_parser.hpp>

#if defined(SIMD)
#include <function/SIMDRecSplit.hpp>
template<size_t LEAF_SIZE>
using RecSplitRotate = bez::function::SIMDRecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, true>;
template<size_t LEAF_SIZE>
using RecSplit = bez::function::SIMDRecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, false>;
std::string name = "SIMDRecSplit";
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE>
using RecSplitRotate = bez::function::GPURecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, true>;
template<size_t LEAF_SIZE>
using RecSplit = bez::function::GPURecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, false>;
std::string name = "GPURecSplit";
#else
#define SHOCKHASH_ENABLED
#include <shockhash/ShockHash.h>
#include <function/RecSplit.hpp>
template<size_t LEAF_SIZE>
using RecSplitRotate = bez::function::RecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, true>;
template<size_t LEAF_SIZE>
using RecSplit = bez::function::RecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, false>;
std::string name = "CpuRecSplit";
#endif

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory");

size_t numObjects = 1e6;
size_t numQueries = 1e6;
std::string leafMethod = "bruteforce";
size_t leafSize = 8;
size_t bucketSize = 1000;
size_t numThreads = 1;

template<typename RecSplit, typename hash128_t>
void construct() {
    auto time = std::chrono::system_clock::now();
    long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    std::cout<<"Generating input data (Seed: "<<seed<<")"<<std::endl;
    util::XorShift64 prng(seed);
	std::vector<hash128_t> keys;
    for (size_t i = 0; i < numObjects; i++) {
        keys.push_back(hash128_t(prng(), prng()));
    }

    std::cout<<"Constructing"<<std::endl;
    auto beginConstruction = std::chrono::high_resolution_clock::now();
    RecSplit rs(keys, bucketSize
                        #if (defined(SIMD) || defined(GPU))
                            , numThreads
                        #endif
                        );
    unsigned long constructionDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginConstruction).count();

    std::cout<<"Querying"<<std::endl;
    uint64_t h = 0;
    auto beginQueries = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numQueries; i++) {
        h ^= rs(hash128_t(prng(), prng() ^ h));
    }
    auto queryDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginQueries).count();
    DO_NOT_OPTIMIZE(h);

    std::cout << "RESULT"
              << " name=" << name
              << " l=" << leafSize
              << " b=" << bucketSize
              << " leafMethod=" << leafMethod
              << " N=" << numObjects
              << " numQueries=" << numQueries
              << " threads=" << numThreads
              << " queryTimeMilliseconds=" << queryDurationMs
              << " constructionTimeMilliseconds=" << constructionDurationMs
              << " bitsPerElement=" << (double) rs.getBits() / numObjects
              << std::endl;
}

template <template<size_t> class RecSplit, typename hash128_t, size_t I>
void dispatchLeafSize(size_t param) {
    if constexpr (I <= 2) {
        std::cerr<<"The parameter "<<param<<" for the leaf size was not compiled into this binary."<<std::endl;
    } else if (I == param) {
        construct<RecSplit<I>, hash128_t>();
    } else {
        dispatchLeafSize<RecSplit, hash128_t, I - 1>(param);
    }
}

int constructAll(int argc, const char* const* argv) {
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numObjects", numObjects, "Number of objects to construct with");
    cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to measure");
    cmd.add_string('L', "leafMethod", leafMethod, "Method to use for leaf bijections");
    cmd.add_bytes('l', "leafSize", leafSize, "Leaf size to construct");
    cmd.add_bytes('b', "bucketSize", bucketSize, "Bucket size to construct");
    cmd.add_bytes('t', "numThreads", numThreads, "Threads to use for construction");

    if (!cmd.process(argc, argv)) {
        return 1;
    }

    if (leafMethod == "bruteforce") {
        dispatchLeafSize<RecSplit, bez::function::hash128_t, 10>(leafSize);
    } else if (leafMethod == "rotations") {
        //dispatchLeafSize<RecSplitRotate, bez::function::hash128_t, bez::function::MAX_LEAF_SIZE>(leafSize);
    }
#ifdef SHOCKHASH_ENABLED
    else if (leafMethod == "cuckoo") {
        dispatchLeafSize<shockhash::ShockHash, sux::function::hash128_t, shockhash::MAX_LEAF_SIZE>(leafSize);
    }
#endif
    else {
        std::cerr<<"Invalid leaf mode argument: "<<leafMethod<<std::endl;
    }
    return 0;
}
