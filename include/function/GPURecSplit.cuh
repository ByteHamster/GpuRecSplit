/*
 * Sux: Succinct data structures
 *
 * Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional permissions
 * described in the GCC Runtime Library Exception, version 3.1, as published by
 * the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and a copy of
 * the GCC Runtime Library Exception along with this program; see the files
 * COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*
 * Differences to RecSplit:
 * - The leaves, lower_aggr and upper_aggr use a constant start_seed independent of the length
 *   from the root. This simplifies the GPU implementation and reduces the query time by a little
 *   bit. It doesn't have a real disadvantage, but the resulting MPHF is different than the one
 *   produced by RecSplit.
 */

#pragma once

#include <sux/util/Vector.hpp>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>

#include "DoubleEF.hpp"
#include "RiceBitVector.hpp"
#include "AbstractParallelRecSplit.hpp"

// Define constexpr namespace ce
#include <gcem.hpp>
namespace ce = gcem;

namespace cg = cooperative_groups;

namespace bez::function {

using namespace std;
using namespace std::chrono;

static_assert(sizeof(uint64_t) == sizeof(unsigned long long int),
    "Since atomicMin does not accept unsigned long int, we need to cast uint64_t to unsigned long long int!");
using ulli = unsigned long long int;

#define checkCudaError(cmd) {_checkCudaError(cmd, #cmd, __FILE__, __LINE__);}
void _checkCudaError(cudaError_t error, const char* cmd, const char* file, int line) {
    if (error != cudaSuccess) {
        printf("CUDA Error %s by command %s in file %s line %d\n", cudaGetErrorString(error), cmd, file, line);
        cudaDeviceReset();
        exit(error);
    }
}

__forceinline__ __device__
static uint32_t remixAndRemap32(uint64_t x, uint32_t n) {
    constexpr int masklen = 16;
    constexpr uint32_t mask = (uint32_t(1) << masklen) - 1;
    return ((remix32(uint32_t(x >> 32) ^ uint32_t(x)) & mask) * n) >> masklen;
}

// Computes the point at which one should stop to test whether
// bijection extraction failed on the GPU (around the square root of the leaf size).
static constexpr array<uint8_t, MAX_LEAF_SIZE + 1> fill_gpu_bij_midstop() {
    array<uint8_t, MAX_LEAF_SIZE + 1> memo{ 0 };
    for (int s = 0; s < MAX_LEAF_SIZE + 1; ++s) memo[s] = s < (int)ce::ceil(3 * ce::sqrt(s)) ? s : (int)ce::ceil(3 * ce::sqrt(s));
    return memo;
}

static constexpr __constant__ auto gpu_bij_midstop = fill_gpu_bij_midstop(); // use flag --expt-relaxed-constexpr for this

static constexpr uint64_t INVALID_RESULT = numeric_limits<uint64_t>::max();

// Only use the sync_frequency in leafKernel for LEAF_SIZE >= SYNC_THRESHOLD.
static constexpr int SYNC_THRESHOLD = 14;

/**
  * Computes one splitting of the higher levels with fanout 2 using a single thread block.
  */
__global__ void higherLevelKernel(uint64_t *__restrict__ device_keys, uint64_t *__restrict__ bucket_size_acc,
                                  uint32_t *__restrict__ results,
                                  size_t *__restrict__ bucketIdxes, size_t offsetInBucket, size_t offsetInResults,
                                  const uint32_t size, const uint32_t split, const uint64_t seed, size_t maxResultsSize) {
    __builtin_assume(size < MAX_BUCKET_SIZE);
    __builtin_assume(split < size);
    __builtin_assume(split >= 8); // since upper_aggr >= 8

    extern __shared__ uint64_t s_bucket[];
    __shared__ uint64_t s_result;
    __shared__ uint32_t writeback_pos[2];

    auto block = cg::this_thread_block();
    size_t bucketIdx = bucketIdxes[blockIdx.y];
    uint64_t *__restrict__ bucket = device_keys + bucket_size_acc[bucketIdx] + offsetInBucket;
    cooperative_groups::memcpy_async(block, s_bucket, bucket, size * sizeof(uint64_t));
    if (threadIdx.x == 0) {
        s_result = INVALID_RESULT;
        writeback_pos[0] = 0;
        writeback_pos[1] = split;
    }
    uint64_t x = seed + threadIdx.x;
    uint32_t count = 0;
    cg::wait(block);

    for (;;) {
        for (uint32_t i = 0; i < size; ++i) {
            count += remixAndRemap32(x + s_bucket[i], size) < split;
        }
        if (count == split) {
            atomicMin(reinterpret_cast<ulli*>(&s_result), static_cast<ulli>(x - seed));
        }
        block.sync();
        if (s_result != INVALID_RESULT) {
            break;
        }
        x += blockDim.x;
        count = 0;
    }

    if (threadIdx.x == 0)
        results[bucketIdx * maxResultsSize + offsetInResults] = s_result;

    for (uint32_t i = threadIdx.x; i < size; i += blockDim.x) {
        const uint64_t value = s_bucket[i];
        const uint32_t child_pos = remixAndRemap32(seed + s_result + value, size) >= split;
        const uint32_t pos = atomicAdd(&writeback_pos[child_pos], 1);
        bucket[pos] = value;
    }
}

/**
  * Computes splittings of the first and second aggregation layer with one thread block
  * per splitting.
  * 
  * \tparam _MAX_FANOUT An upper limit for the given fanout
  */
template<int _MAX_FANOUT, int BLOCK_SIZE, uint32_t SPLIT, uint64_t SEED>
__global__ void aggrKernel(uint64_t *__restrict__ g_keys, uint64_t *__restrict__ bucket_size_acc,
                           uint32_t *__restrict__ results,
                           size_t *bucketIdxes, uint64_t offsetInBucket, uint64_t offsetInResults,
                           const uint32_t size, const uint32_t fanout, size_t maxResultsSize) {
    static_assert(_MAX_FANOUT <= 9, "aggrKernel only works for _MAX_FANOUT <= 9!");
    __builtin_assume(fanout > 1);
    __builtin_assume(fanout <= _MAX_FANOUT);
    __builtin_assume(size <= _MAX_FANOUT * SPLIT);
    __builtin_assume(SPLIT < size);

    extern __shared__ uint64_t s_bucket[];
    // If SPLIT is less than 255 (which is always true given LEAF_SIZE <= 24) 8 bits are sufficient
    // to hold the counts. Note that it may overflow, but in that case the count cannot reach SPLIT
    // again and therefore is correctly identified as an unsuitable split.
    static_assert(SPLIT < 255, "SPLIT must be less than 255 for aggrKernel to work correctly!"
        "Note that less than 256 is not enough since an overflow may carry to the next count.");
    constexpr uint32_t FULL_FOUND_32 = (SPLIT << 24) | (SPLIT << 16) | (SPLIT << 8) | SPLIT;
    __shared__ uint64_t s_result;
    __shared__ uint32_t writeback_pos[_MAX_FANOUT];

    auto block = cg::this_thread_block();
    size_t bucketIdx = bucketIdxes[blockIdx.y];
    uint64_t *__restrict__ bucket = g_keys + bucket_size_acc[bucketIdx] + offsetInBucket + blockIdx.x * size;
    cooperative_groups::memcpy_async(block, s_bucket, bucket, size * sizeof(uint64_t));
    if (threadIdx.x < fanout) {
        if (threadIdx.x == 0)
            s_result = INVALID_RESULT;
        writeback_pos[threadIdx.x] = threadIdx.x * SPLIT;
    }
    uint64_t x = SEED + threadIdx.x;

    if (fanout <= 5) {
        uint32_t count = 0;
        uint32_t found = SPLIT;
        for (uint32_t i = 1; i < fanout - 1; ++i)
            found |= SPLIT << (8 * i);
        if (fanout <= 4)
            found |= (size - (fanout - 1) * SPLIT) << (8 * (fanout - 1));
        cg::wait(block);

        for (;;) {
            for (uint32_t i = 0; i < size; ++i) {
                uint32_t shift = ((remixAndRemap32(x + s_bucket[i], size) / SPLIT) << 3);
                if constexpr (_MAX_FANOUT >= 5)
                    count += __funnelshift_lc(0U, 1U, shift); // is zero for shift >= 32
                else
                    count += 1U << shift;
            }
            if (__builtin_expect(count == found, 0)) {
                atomicMin(reinterpret_cast<ulli*>(&s_result), static_cast<ulli>(x - SEED));
            }

            block.sync();
            if (s_result != INVALID_RESULT) {
                break;
            }
            x += blockDim.x;
            count = 0;
        }
    } else {
        uint64_t count = 0;
        uint32_t found_high = SPLIT;
        for (uint32_t i = 1; i < fanout - 5; ++i)
            found_high |= SPLIT << (8 * i);
        if (fanout <= 8)
            found_high |= (size - (fanout - 1) * SPLIT) << (8 * (fanout - 5));
        uint64_t found = ((uint64_t)found_high << 32) | FULL_FOUND_32;
        cg::wait(block);

        for (;;) {
            for (uint32_t i = 0; i < size; ++i) {
                uint32_t shift = ((remixAndRemap32(x + s_bucket[i], size) / SPLIT) << 3);
                if constexpr (_MAX_FANOUT == 9)
                    count += shift < 64 ? 1ULL << shift : 0;
                else
                    count += 1ULL << shift;
            }
            if (__builtin_expect(count == found, 0)) {
                atomicMin(reinterpret_cast<ulli *>(&s_result), static_cast<ulli>(x - SEED));
            }

            block.sync();
            if (s_result != INVALID_RESULT) {
                break;
            }
            x += blockDim.x;
            count = 0;
        }
    }

    if (threadIdx.x == 0)
        results[bucketIdx * maxResultsSize + offsetInResults + blockIdx.x] = s_result;

    for (uint32_t i = threadIdx.x; i < size; i += blockDim.x) {
        const uint64_t value = s_bucket[i];
        const uint32_t child_pos = remixAndRemap32(SEED + s_result + value, size) / SPLIT;
        const uint32_t pos = atomicAdd(&writeback_pos[child_pos], 1);
        bucket[pos] = value;
    }
}

template<uint32_t SIZE>
__forceinline__ __host__ __device__
static constexpr uint32_t rotate(uint32_t val, uint32_t x) {
    return ((val << x) | (val >> (SIZE - x))) & ((1 << SIZE) - 1);
}

template<uint32_t MAX_SIZE, bool USE_ROTATE>
__global__ void leafKernel(const uint64_t *__restrict__ g_keys, const uint64_t *__restrict__ bucket_size_acc,
                           uint32_t *__restrict__ results,
                           size_t *bucketIdxes, size_t offsetInBucket, size_t offsetInResults, const uint32_t size,
                           const int sync_frequency, size_t maxResultsSize) {
    assert(!USE_ROTATE || size == MAX_SIZE);
    __builtin_assume(!USE_ROTATE || size == MAX_SIZE);
    __builtin_assume(size <= MAX_SIZE);
    __builtin_assume(size > 1);
    __builtin_assume(MAX_SIZE >= SYNC_THRESHOLD || sync_frequency == 0); // not useful for small leaves

    constexpr uint64_t SEED = start_seed[NUM_START_SEEDS - 1];
    __shared__ uint64_t s_bucket[MAX_SIZE];
    __shared__ uint64_t s_result;

    [[maybe_unused]] __shared__ int size_left;
    [[maybe_unused]] __shared__ int size_right;
    [[maybe_unused]] __shared__ uint64_t s_bucket_left[USE_ROTATE ? MAX_SIZE : 1]; // Unfortunately, size 0 is not allowed
    [[maybe_unused]] __shared__ uint64_t s_bucket_right[USE_ROTATE ? MAX_SIZE : 1];

    auto block = cg::this_thread_block();
    size_t bucketIdx = bucketIdxes[blockIdx.y];
    const uint64_t *__restrict__ bucket = g_keys + bucket_size_acc[bucketIdx] + offsetInBucket + blockIdx.x * size;
    cooperative_groups::memcpy_async(block, s_bucket, bucket, size * sizeof(uint64_t));
    if (threadIdx.x == 0) {
        s_result = INVALID_RESULT;
        if constexpr (USE_ROTATE) {
            size_left = size_right = 0;
        }
    }
    uint64_t x = SEED + (USE_ROTATE ? MAX_SIZE : 1U) * threadIdx.x;
    uint32_t mask_left = 0;
    [[maybe_unused]] uint32_t mask_right = 0;
    const uint32_t found = (1U << size) - 1;
    [[maybe_unused]] int sync_counter = 0;
    cg::wait(block);

    if constexpr (USE_ROTATE) {
        assert(blockDim.x >= MAX_SIZE);
        uint64_t element;
        if (threadIdx.x < MAX_SIZE) {
            element = s_bucket[threadIdx.x];
            if (element % 2 == 0) {
                int pos = atomicAdd(&size_left, 1);
                s_bucket_left[pos] = element;
            } else {
                int pos = atomicAdd(&size_right, 1);
                s_bucket_right[pos] = element;
            }
        }
        block.sync();

        for (;;) {
            for (int i = 0; i < size_left; ++i) {
                mask_left |= 1U << remixAndRemap32(x + s_bucket_left[i], MAX_SIZE);
            }
            if (__popc(mask_left) == size_left) {
                for (int i = 0; i < size_right; ++i) {
                    mask_right |= 1U << remixAndRemap32(x + s_bucket_right[i], MAX_SIZE);
                }
                if (__popc(mask_right) == size_right) {
                    // Try to rotate right part to see if both together form a bijection
                    for (uint32_t rotations = 0; rotations < MAX_SIZE; ++rotations) {
                        if ((mask_left | mask_right) == found) {
                            atomicMin(reinterpret_cast<ulli*>(&s_result), static_cast<ulli>(x + rotations - SEED));
                            break;
                        }
                        mask_right = rotate<MAX_SIZE>(mask_right, 1);
                    }
                }
            }

            block.sync();
            if (s_result != INVALID_RESULT) {
                break;
            }
            x += blockDim.x * MAX_SIZE;
            mask_left = mask_right = 0;
        }
    } else if constexpr (MAX_SIZE < 14) { // no midstop, no sync_counter
        for (;;) {
            for (uint32_t i = 0; i < size; ++i) {
                mask_left |= 1U << remixAndRemap32(x + s_bucket[i], size);
            }
            if (__builtin_expect(mask_left == found, 0)) {
                atomicMin(reinterpret_cast<ulli*>(&s_result), static_cast<ulli>(x - SEED));
            }

            block.sync();
            if (s_result != INVALID_RESULT) {
                break;
            }
            x += blockDim.x;
            mask_left = 0;
        }
    } else {
        const uint32_t midstop = gpu_bij_midstop[size];
        for (;;) {
            uint32_t i;
            for (i = 0; i < midstop; ++i) {
                mask_left |= 1U << remixAndRemap32(x + s_bucket[i], size);
            }
            if (__builtin_expect(__popc(mask_left) == midstop, 0)) {
                for (; i < size; ++i) {
                    mask_left |= 1U << remixAndRemap32(x + s_bucket[i], size);
                }
                if (mask_left == found) {
                    atomicMin(reinterpret_cast<ulli*>(&s_result), static_cast<ulli>(x - SEED));
                }
            }
            if (++sync_counter >= sync_frequency) {
                sync_counter = 0;
                block.sync();
                if (s_result != INVALID_RESULT) {
                    break;
                }
            }
            x += blockDim.x;
            mask_left = 0;
        }
    }

    if (threadIdx.x == 0)
        results[bucketIdx * maxResultsSize + offsetInResults + blockIdx.x] = s_result;
}

/**
 *
 * A class for storing minimal perfect hash functions. The template
 * parameter decides how large a leaf will be. Larger leaves imply
 * slower construction, but less space and faster evaluation.
 *
 * @tparam LEAF_SIZE the size of a leaf; typicals value range from 6 to 8
 * for fast, small maps, or up to 16 for very compact functions.
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC, bool USE_BIJECTIONS_ROTATE = true>
class GPURecSplit
    : public AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, false> {
    using Superclass = AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, false>;
    using Superclass::SplitStrat;
    using Superclass::_leaf;
    using Superclass::lower_aggr;
    using Superclass::upper_aggr;
    using Superclass::memo;
    using Superclass::use_bijections_rotate;

    static constexpr int NUM_STREAMS = 128;
    static constexpr int HIGHER_LEVEL_BLOCK_SIZE = 64;
    static constexpr int UPPER_AGGR_BLOCK_SIZE = 256;
    static constexpr int LOWER_AGGR_BLOCK_SIZE = 256;
    static constexpr int LEAF_BLOCK_SIZE = 512;
    static constexpr size_t CUDA_MAX_GRID = 65535;

    const size_t numThreads = 8;
    uint64_t *device_keys;
    uint64_t *device_bucket_size_acc;
    uint32_t *device_results;
  public:
    GPURecSplit() {}

    /** Builds a GPURecSplit instance using a given list of keys and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param keys a vector of strings.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    GPURecSplit(const vector<string> &keys, const size_t bucket_size, const int numThreads)
            : numThreads(numThreads) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash128_t *h = (hash128_t *)malloc(this->keys_count * sizeof(hash128_t));

        if (numThreads == 1) {
            for (size_t i = 0; i < this->keys_count; ++i) {
                h[i] = first_hash(keys[i].c_str(), keys[i].size());
            }
        } else {
            size_t keysPerThread = this->keys_count / numThreads + 1;
            std::vector<std::thread> threads;
            for (size_t thread = 0; thread < numThreads; thread++) {
                threads.emplace_back([&, thread] {
                    size_t from = thread * keysPerThread;
                    size_t to = std::min(this->keys_count, (thread + 1) * keysPerThread);
                    for (size_t i = from; i < to; ++i) {
                        h[i] = first_hash(keys[i].c_str(), keys[i].size());
                    }
                });
            }
            for (std::thread &t : threads) {
                t.join();
            }
        }

        hash_gen(h);
        free(h);
    }

    /** Builds a GPURecSplit instance using a given list of 128-bit hashes and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * Note that this constructor is mainly useful for benchmarking.
     * @param keys a vector of 128-bit hashes.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    GPURecSplit(vector<hash128_t> &keys, const size_t bucket_size, const int numThreads)
            : numThreads(numThreads) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash_gen(&keys[0]);
    }

    /** Builds a GPURecSplit instance using a list of keys returned by a stream and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param input an open input stream returning a list of keys, one per line.
     * @param bucket_size the desired bucket size.
     */
    GPURecSplit(ifstream& input, const size_t bucket_size, const int numThreads)
            : numThreads(numThreads) {
        this->bucket_size = bucket_size;
        vector<hash128_t> h;
        for(string key; getline(input, key);) h.push_back(first_hash(key.c_str(), key.size()));
        this->keys_count = h.size();
        hash_gen(&h[0]);
    }

  private:
    array<uint32_t, 3> executeBucketKernels(uint32_t size, size_t *bucketIdxes, size_t numBuckets, cudaStream_t stream, size_t maxResultsSize) {

        uint32_t num_results = executeHigherLevels(size, bucketIdxes, numBuckets, 0, 0, stream, maxResultsSize);
        array<uint32_t, 3> result_counts;
        result_counts[0] = num_results;

        constexpr uint32_t upper_fanout = upper_aggr / lower_aggr;
        if (size >= upper_aggr) {
            const int num_blocks = size / upper_aggr;
            aggrKernel<upper_fanout, UPPER_AGGR_BLOCK_SIZE, lower_aggr, start_seed[NUM_START_SEEDS - 3]>
                <<<dim3(num_blocks, numBuckets), UPPER_AGGR_BLOCK_SIZE, upper_aggr * sizeof(uint64_t), stream>>>(
                device_keys, device_bucket_size_acc, device_results, bucketIdxes, 0, num_results, upper_aggr, upper_fanout, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += num_blocks;
        }
        uint32_t remaining = size % upper_aggr;
        if (remaining > lower_aggr) {
            aggrKernel<upper_fanout, UPPER_AGGR_BLOCK_SIZE, lower_aggr, start_seed[NUM_START_SEEDS - 3]>
                <<<dim3(1, numBuckets), UPPER_AGGR_BLOCK_SIZE, remaining * sizeof(uint64_t), stream>>>(device_keys, device_bucket_size_acc,
                 device_results, bucketIdxes, size - remaining,
                num_results, remaining, uint16_t(remaining + lower_aggr - 1) / lower_aggr, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += 1;
        }
        result_counts[1] = num_results;

        constexpr uint32_t lower_fanout = lower_aggr / LEAF_SIZE;
        if (size >= lower_aggr) {
            const int num_blocks = size / lower_aggr;
            aggrKernel<lower_fanout, LOWER_AGGR_BLOCK_SIZE, LEAF_SIZE, start_seed[NUM_START_SEEDS - 2]>
                <<<dim3(num_blocks, numBuckets), LOWER_AGGR_BLOCK_SIZE, lower_aggr * sizeof(uint64_t), stream>>>(
                device_keys, device_bucket_size_acc, device_results, bucketIdxes, 0, num_results, lower_aggr, lower_fanout, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += num_blocks;
        }
        remaining = size % lower_aggr;
        if (remaining > LEAF_SIZE) {
            aggrKernel<lower_fanout, LOWER_AGGR_BLOCK_SIZE, LEAF_SIZE, start_seed[NUM_START_SEEDS - 2]>
                <<<dim3(1, numBuckets), LOWER_AGGR_BLOCK_SIZE, remaining * sizeof(uint64_t), stream>>>(device_keys,
                     device_bucket_size_acc, device_results, bucketIdxes, size - remaining,
                     num_results, remaining, uint16_t(remaining + LEAF_SIZE - 1) / LEAF_SIZE, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += 1;
        }
        result_counts[2] = num_results;

        constexpr uint32_t SYNC_CONST = LEAF_SIZE < SYNC_THRESHOLD ? 10'000'000 : 32;
        if (size >= LEAF_SIZE) {
            const int num_blocks = size / LEAF_SIZE;
            constexpr uint32_t sync_frequency = (1U << golomb_param(LEAF_SIZE)) / LEAF_BLOCK_SIZE / SYNC_CONST;
            // for __builtin_assume in leafKernel
            static_assert(LEAF_SIZE >= SYNC_THRESHOLD || sync_frequency == 0, "sync_frequency must be 0 for LEAF_SIZE < SYNC_THRESHOLD!");
            leafKernel<LEAF_SIZE, use_bijections_rotate><<<dim3(num_blocks, numBuckets), LEAF_BLOCK_SIZE, 0, stream>>>(device_keys, device_bucket_size_acc,
                 device_results, bucketIdxes, 0, num_results, LEAF_SIZE, sync_frequency, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += num_blocks;
        }
        remaining = size % LEAF_SIZE;
        if (remaining > 1) {
            const uint32_t sync_frequency = (1U << golomb_param(remaining)) / LEAF_BLOCK_SIZE / SYNC_CONST;
            assert(LEAF_SIZE >= SYNC_THRESHOLD || sync_frequency == 0); // for __builtin_assume in leafKernel
            leafKernel<LEAF_SIZE, false><<<dim3(1, numBuckets), LEAF_BLOCK_SIZE, 0, stream>>>(device_keys, device_bucket_size_acc,
                device_results, bucketIdxes, size - remaining,
                num_results, remaining, sync_frequency, maxResultsSize);
            checkCudaError(cudaGetLastError());
            num_results += 1;
        }
        assert(num_results < maxResultsSize);
        return result_counts;
    }

    uint32_t executeHigherLevels(size_t size, size_t *bucketIdxes, size_t numBuckets,
                     size_t offsetInBucket, size_t offsetInResults, cudaStream_t stream, size_t maxResultsSize, const int level = 0) {
        if (size > upper_aggr) {
            const uint32_t split = ((uint16_t(size / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
            higherLevelKernel<<<dim3(1, numBuckets), HIGHER_LEVEL_BLOCK_SIZE, size * sizeof(uint64_t), stream>>>(device_keys, device_bucket_size_acc, device_results,
                bucketIdxes, offsetInBucket, offsetInResults, size, split, start_seed[level], maxResultsSize);
            checkCudaError(cudaGetLastError());
            uint32_t skip_results = 1 + executeHigherLevels(split, bucketIdxes, numBuckets, offsetInBucket, offsetInResults + 1, stream, maxResultsSize, level + 1);
            skip_results += executeHigherLevels(size - split, bucketIdxes, numBuckets, offsetInBucket + split, offsetInResults + skip_results, stream, maxResultsSize, level + 1);
            return skip_results;
        }
        return 0;
    }

    void unpackResults(const size_t size, uint32_t *result, typename RiceBitVector<AT>::Builder &builder,
        vector<uint32_t> &unary, array<uint32_t, 3> &result_counts) {
        int higher_level_pos = 0;
        unpackResults(size, result, builder, unary, result_counts, higher_level_pos);
    }

    void unpackResults(const size_t size, uint32_t *result, typename RiceBitVector<AT>::Builder &builder,
                       vector<uint32_t> &unary, array<uint32_t, 3>& result_counts, int &higher_level_pos) {
        const auto log2golomb = golomb_param(size);
        if (size > upper_aggr) {
            const uint64_t x = result[higher_level_pos++];
            builder.appendFixed(x, log2golomb);
            unary.push_back(x >> log2golomb);
            const size_t split = ((uint16_t(size / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
            unpackResults(split, result, builder, unary, result_counts, higher_level_pos);
            if (size - split > 1)
                unpackResults(size - split, result, builder, unary, result_counts, higher_level_pos);
        } else if (size > lower_aggr) {
            const uint64_t x = result[result_counts[0]++];
            builder.appendFixed(x, log2golomb);
            unary.push_back(x >> log2golomb);
            size_t i;
            for (i = 0; i < size - lower_aggr; i += lower_aggr) {
                unpackResults(lower_aggr, result, builder, unary, result_counts, higher_level_pos);
            }
            if (size - i > 1)
                unpackResults(size - i, result, builder, unary, result_counts, higher_level_pos);
        } else if (size > LEAF_SIZE) {
            const uint64_t x = result[result_counts[1]++];
            builder.appendFixed(x, log2golomb);
            unary.push_back(x >> log2golomb);
            size_t i;
            for (i = 0; i < size - LEAF_SIZE; i += LEAF_SIZE) {
                unpackResults(LEAF_SIZE, result, builder, unary, result_counts, higher_level_pos);
            }
            if (size - i > 1)
                unpackResults(size - i, result, builder, unary, result_counts, higher_level_pos);
        } else {
            const uint64_t x = result[result_counts[2]++];
            builder.appendFixed(x, log2golomb);
            unary.push_back(x >> log2golomb);
        }
    }

    void hash_gen(hash128_t *hashes) {
#ifdef MORESTATS
        time_bij = 0;
        memset(time_split, 0, sizeof time_split);
        split_unary = split_fixed = 0;
        bij_unary = bij_fixed = 0;
        min_split_code = 1ULL << 63;
        max_split_code = sum_split_codes = 0;
        min_bij_code = 1ULL << 63;
        max_bij_code = sum_bij_codes = 0;
        sum_depths = 0;
        size_t minsize = this->keys_count, maxsize = 0;
        double ub_split_bits = 0, ub_bij_bits = 0;
        double ub_split_evals = 0;

        auto total_start_time = high_resolution_clock::now();
#endif

#ifndef __SIZEOF_INT128__
        if (this->keys_count > (1ULL << 32)) {
            fprintf(stderr, "For more than 2^32 keys, you need 128-bit integer support.\n");
            abort();
        }
#endif

        this->nbuckets = max(1, (this->keys_count + this->bucket_size - 1) / this->bucket_size);
        auto bucket_size_acc = vector<uint64_t>(this->nbuckets + 1);
        auto bucket_pos_acc = vector<uint64_t>(this->nbuckets + 1);
        auto sorted_keys = vector<uint64_t>(this->keys_count);

#ifdef MORESTATS
        auto sorting_start_time = high_resolution_clock::now();
#endif
        this->parallelPartition(hashes, sorted_keys, bucket_size_acc, numThreads);
#ifdef MORESTATS
        auto sorting_end_time = high_resolution_clock::now();
        auto sorting_time = duration_cast<nanoseconds>(sorting_end_time - sorting_start_time).count();
#endif
        std::vector<std::vector<size_t>> bucketsBySize;
        bucketsBySize.resize(MAX_BUCKET_SIZE);

        size_t maxBucketSize = 0;
        size_t maxNumBucketsWithSameSize = 0;
        for (size_t i = 0; i < this->nbuckets; i++) {
            size_t s = bucket_size_acc.at(i + 1) - bucket_size_acc.at(i);
            bucketsBySize.at(s).push_back(i);
            maxBucketSize = std::max(maxBucketSize, s);
            maxNumBucketsWithSameSize = std::max(maxNumBucketsWithSameSize, bucketsBySize.at(s).size());
#ifdef MORESTATS
            /*auto upper_leaves = (s + _leaf - 1) / _leaf;
            auto upper_height = ceil(log(upper_leaves) / log(2)); // TODO: check
            auto upper_s = _leaf * pow(2, upper_height);
            ub_split_bits += (double)upper_s / (_leaf * 2) * log2(2 * M_PI * _leaf) - .5 * log2(2 * M_PI * upper_s);
            ub_bij_bits += upper_leaves * _leaf * (log2e - .5 / _leaf * log2(2 * M_PI * _leaf));
            ub_split_evals += 4 * upper_s * sqrt(pow(2 * M_PI * upper_s, 2 - 1) / pow(2, 2));
            minsize = min(minsize, s);
            maxsize = max(maxsize, s);*/
#endif
        }

#ifdef MORESTATS
        auto bucket_by_size_end_time = high_resolution_clock::now();
        auto bucket_by_size_time = duration_cast<nanoseconds>(bucket_by_size_end_time - sorting_end_time).count();
#endif

        typename RiceBitVector<AT>::Builder builder;

        std::vector<cudaStream_t> streams;
        streams.resize(NUM_STREAMS);
        for (cudaStream_t &stream : streams) {
            checkCudaError(cudaStreamCreate(&stream));
        }

        size_t maxResultsSize = 2 * (maxBucketSize / LEAF_SIZE + 1);
        const size_t results_size = this->nbuckets * maxResultsSize * sizeof(uint32_t);
        uint32_t *host_results;
        checkCudaError(cudaMallocHost(&host_results, results_size));
        checkCudaError(cudaMalloc(&device_keys, this->keys_count * sizeof(uint64_t)));
        checkCudaError(cudaMalloc(&device_bucket_size_acc, this->nbuckets * sizeof(uint64_t)));
        checkCudaError(cudaMalloc(&device_results, results_size));
        size_t *bucketIdxes;
        size_t bucketIdxesStreamOffset = std::min(maxNumBucketsWithSameSize, CUDA_MAX_GRID);
        checkCudaError(cudaMalloc(&bucketIdxes, bucketIdxesStreamOffset * NUM_STREAMS * sizeof(size_t)));

        checkCudaError(cudaMemcpyAsync(device_keys, sorted_keys.data(), this->keys_count * sizeof(uint64_t), cudaMemcpyHostToDevice, streams[0]));
        checkCudaError(cudaMemcpyAsync(device_bucket_size_acc, bucket_size_acc.data(), this->nbuckets * sizeof(uint64_t), cudaMemcpyHostToDevice, streams[1]));
        cudaDeviceSynchronize();

#ifdef MORESTATS
        auto alloc_end_time = high_resolution_clock::now();
        auto alloc_time = duration_cast<nanoseconds>(alloc_end_time - bucket_by_size_end_time).count();
#endif

        vector<array<uint32_t, 3>> result_counts(MAX_BUCKET_SIZE);
        for (size_t bucketSize = 2; bucketSize < MAX_BUCKET_SIZE; bucketSize++) {
            if (bucketsBySize.at(bucketSize).empty()) {
                continue;
            }
            auto &buckets = bucketsBySize.at(bucketSize);
            size_t numBuckets = buckets.size();
            for (size_t from = 0; from < numBuckets; from += CUDA_MAX_GRID) {
                size_t batchSize = std::min(numBuckets - from, CUDA_MAX_GRID);
                size_t streamIdx = bucketSize % NUM_STREAMS; // Round-robin
                checkCudaError(cudaMemcpyAsync(bucketIdxes + streamIdx * bucketIdxesStreamOffset, buckets.data() + from, batchSize * sizeof(size_t), cudaMemcpyHostToDevice, streams[streamIdx]));
                result_counts[bucketSize] = executeBucketKernels(bucketSize, bucketIdxes + streamIdx * bucketIdxesStreamOffset,
                                     batchSize, streams[streamIdx], maxResultsSize);
            }
        }

#ifdef MORESTATS
        auto cuda_submit_end_time = high_resolution_clock::now();
        auto cuda_submit_time = duration_cast<nanoseconds>(cuda_submit_end_time - alloc_end_time).count();
#endif

        checkCudaError(cudaDeviceSynchronize());
        checkCudaError(cudaMemcpy(host_results, device_results, results_size, cudaMemcpyDeviceToHost));

#ifdef MORESTATS
        auto cuda_end_time = high_resolution_clock::now();
        auto cuda_time = duration_cast<nanoseconds>(cuda_end_time - cuda_submit_end_time).count();
#endif

        bucket_pos_acc[0] = 0;
        vector<uint32_t> unary;
        for (size_t i = 0; i < this->nbuckets; i++) {
            size_t bucketSize = bucket_size_acc.at(i + 1) - bucket_size_acc.at(i);
            auto result_counts_copy = result_counts[bucketSize];
            unpackResults(bucketSize, host_results + i * maxResultsSize, builder, unary, result_counts_copy);
            builder.appendUnaryAll(unary);
            unary.clear();
            bucket_pos_acc[i + 1] = builder.getBits();
        }

#ifdef MORESTATS
        auto unpack_end_time = high_resolution_clock::now();
        auto unpack_time = duration_cast<nanoseconds>(unpack_end_time - cuda_end_time).count();
#endif

        checkCudaError(cudaFreeHost(host_results));
        checkCudaError(cudaDeviceReset());
        builder.appendFixed(1, 1); // Sentinel (avoids checking for parts of size 1)
        this->descriptors = builder.build();
        this->ef = DoubleEF<AT, true>(bucket_size_acc, bucket_pos_acc);

#ifdef STATS
        // Evaluation purposes only
        double ef_sizes = (double)this->ef.bitCountCumKeys() / this->keys_count;
        double ef_bits = (double)this->ef.bitCountPosition() / this->keys_count;
        double rice_desc = (double)builder.getBits() / this->keys_count;
        printf("Elias-Fano cumul sizes:  %f bits/bucket\n", (double)this->ef.bitCountCumKeys() / this->nbuckets);
        printf("Elias-Fano cumul bits:   %f bits/bucket\n", (double)this->ef.bitCountPosition() / this->nbuckets);
        printf("Elias-Fano cumul sizes:  %f bits/key\n", ef_sizes);
        printf("Elias-Fano cumul bits:   %f bits/key\n", ef_bits);
        printf("Rice-Golomb descriptors: %f bits/key\n", rice_desc);
        printf("Total bits:              %f bits/key\n", ef_sizes + ef_bits + rice_desc);
#endif
#ifdef MORESTATS

        printf("\n");
        printf("Min bucket size: %lu\n", minsize);
        printf("Max bucket size: %lu\n", maxsize);

        printf("\n");
        printf("Total time:          %11.3f ms\n", duration_cast<nanoseconds>(high_resolution_clock::now() - total_start_time).count() * 1E-6);
        printf("Sorting time:        %11.3f ms\n", sorting_time * 1E-6);
        printf("Bucket sorting time: %11.3f ms\n", bucket_by_size_time * 1E-6);
        printf("Alloc time:          %11.3f ms\n", alloc_time * 1E-6);
        printf("Cuda submit time:    %11.3f ms\n", cuda_submit_time * 1E-6);
        printf("Cuda wait time:      %11.3f ms\n", cuda_time * 1E-6);
        printf("Unpack time:         %11.3f ms\n", unpack_time * 1E-6);
        printf("Bijections:          %11.3f ms\n", time_bij * 1E-6);
        for (int i = 0; i < MAX_LEVEL_TIME; i++) {
            if (time_split[i] > 0) {
                printf("Split level %d:      %11.3f ms\n", i, time_split[i] * 1E-6);
            }
        }

#endif
    }

    friend ostream &operator<<(ostream &os, const GPURecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE> &rs) {
        const size_t leaf_size = LEAF_SIZE;
        os.write((char *)&leaf_size, sizeof(leaf_size));
        os.write((char *)&rs.bucket_size, sizeof(rs.bucket_size));
        os.write((char *)&rs.keys_count, sizeof(rs.keys_count));
        os << rs.descriptors;
        os << rs.ef;
        return os;
    }

    friend istream &operator>>(istream &is, GPURecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE> &rs) {
        size_t leaf_size;
        is.read((char *)&leaf_size, sizeof(leaf_size));
        if (leaf_size != LEAF_SIZE) {
            fprintf(stderr, "Serialized leaf size %d, code leaf size %d\n", int(leaf_size), int(LEAF_SIZE));
            abort();
        }
        is.read((char *)&rs.bucket_size, sizeof(rs.bucket_size));
        is.read((char *)&rs.keys_count, sizeof(rs.keys_count));
        rs.nbuckets = max(1, (rs.keys_count + rs.bucket_size - 1) / rs.bucket_size);

        is >> rs.descriptors;
        is >> rs.ef;
        return is;
    }
};

} // namespace sux::function
