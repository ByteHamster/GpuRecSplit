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
template<size_t LEAF_SIZE, bool USE_BIJECTIONS_ROTATE>
using RecSplit = bez::function::SIMDRecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, USE_BIJECTIONS_ROTATE>;
std::string name = "SimdRecSplit";
#elif defined(GPU)
#include <function/GPURecSplit.cuh>
template<size_t LEAF_SIZE, bool USE_BIJECTIONS_ROTATE>
using RecSplit = bez::function::GPURecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, USE_BIJECTIONS_ROTATE>;
std::string name = "GpuRecSplit";
#else
#include <function/RecSplit.hpp>
template<size_t LEAF_SIZE, bool USE_BIJECTIONS_ROTATE>
using RecSplit = bez::function::RecSplit<LEAF_SIZE, bez::util::AllocType::MALLOC, USE_BIJECTIONS_ROTATE>;
std::string name = "RecSplit";
#endif

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory");

size_t numObjects = 1e6;
size_t numQueries = 1e6;
bool rotations = false;
size_t leafSize = 8;
size_t bucketSize = 1000;

template<typename RecSplit>
void construct() {
    auto time = std::chrono::system_clock::now();
    long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    std::cout<<"Generating input data (Seed: "<<seed<<")"<<std::endl;
    util::XorShift64 prng(seed);
	std::vector<bez::function::hash128_t> keys;
    for (size_t i = 0; i < numObjects; i++) {
        keys.push_back(bez::function::hash128_t(prng(), prng()));
    }

    std::cout<<"Constructing"<<std::endl;
    auto beginConstruction = std::chrono::high_resolution_clock::now();
    RecSplit rs(keys, bucketSize
                        #if defined(SIMD)
                            , num_threads
                        #endif
                        );
    unsigned long constructionDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginConstruction).count();

    std::cout<<"Querying"<<std::endl;
    uint64_t h = 0;
    auto beginQueries = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numQueries; i++) {
        h ^= rs(bez::function::hash128_t(prng(), prng() ^ h));
    }
    auto queryDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginQueries).count();
    DO_NOT_OPTIMIZE(h);

    std::cout << "RESULT"
              << " name=" << name
              << " l=" << leafSize
              << " b=" << bucketSize
              << " rotations=" << rotations
              << " numObjects=" << numObjects
              << " numQueries=" << numQueries
              << " queryDurationMs=" << queryDurationMs
              << " constructionDurationMs=" << constructionDurationMs
              << std::endl;
}

template <size_t I>
void dispatchLeafSize(size_t param) {
    if constexpr (I <= 2) {
        std::cerr<<"The parameter "<<param<<" for the leaf size was not compiled into this binary."<<std::endl;
    } else if (I == param) {
        if (rotations) {
            construct<RecSplit<I, true>>();
        } else {
            construct<RecSplit<I, false>>();
        }
    } else {
        dispatchLeafSize<I - 1>(param);
    }
}

int constructAll(int argc, const char* const* argv) {
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numObjects", numObjects, "Number of objects to construct with");
    cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to measure");
    cmd.add_bool('r', "rotations", rotations, "Enable rotations");
    cmd.add_bytes('l', "leafSize", leafSize, "Leaf size to construct");
    cmd.add_bytes('b', "bucketSize", bucketSize, "Bucket size to construct");

    if (!cmd.process(argc, argv)) {
        return 1;
    }
    dispatchLeafSize<bez::function::MAX_LEAF_SIZE>(leafSize);
    return 0;
}
