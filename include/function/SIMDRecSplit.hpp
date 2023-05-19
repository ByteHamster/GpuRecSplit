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

#include "../util/Vector.hpp"
#include "../util/SimdUtils.hpp"
#include "DoubleEF.hpp"
#include "RiceBitVector.hpp"
#include "AbstractParallelRecSplit.hpp"
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <mutex>
#include <condition_variable>
#include <thread>

// Define constexpr namespace ce
#include <gcem.hpp>
namespace ce = gcem;

#ifndef SIMD
#error Need to compile SIMDRecSplit with -DSIMD
#endif

#ifdef SIMDRS_512_BIT
constexpr double MIDSTOP_FACTOR = 2.6;
#else
constexpr double MIDSTOP_FACTOR = 2.3;
#endif

#ifdef NDEBUG
#if defined(__clang__)
#define SIMDRS_ASSUME(condition) __builtin_assume(condition)
#elif defined(_MSC_VER)
#define SIMDRS_ASSUME(condition) __assume(condition)
#else
// Use with caution! Especially, condition must not have side effects!
// https://stackoverflow.com/a/26195434
#define SIMDRS_ASSUME(condition) do { if (!(condition)) __builtin_unreachable(); } while (0)
#endif
#else
#define SIMDRS_ASSUME(condition) assert(condition)
#endif

#if defined(__GNUC__) || defined(__GNUG__) || defined(__clang__)
#define SIMDRS_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define SIMDRS_RESTRICT __restrict
#else
#define SIMDRS_RESTRICT /* no op */
#endif

namespace bez::function {

using namespace std;
using namespace std::chrono;

/**
  * 32-bit finalizer function in Austin Appleby's MurmurHash3 (https://github.com/aappleby/smhasher).
  */
FullVecUi inline remix32(FullVecUi z) {
    z ^= z >> 16;
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae35;
    z ^= z >> 16;
    return z;
}

FullVecUi remap(FullVecUq x, FullVecUq y, uint32_t n) {
    constexpr int masklen = 16;
    constexpr uint32_t mask = (uint32_t(1) << masklen) - 1;
    //FullVecUi combined(compress(x >> 32, y >> 32) ^ compress(x, y)); // This is a bit slower than below
    const FullVecUi xx(x);
    const FullVecUi yy(y);
#ifdef SIMDRS_512_BIT
    FullVecUi combined = blend16<0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30>(xx, yy)
        ^ blend16<1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>(xx, yy);
#else
    FullVecUi combined = blend8<0, 2, 4, 6, 8, 10, 12, 14>(xx, yy) ^ blend8<1, 3, 5, 7, 9, 11, 13, 15>(xx, yy);
#endif
    return ((remix32(combined) & mask) * n) >> masklen;
}

// Computes the point at which one should stop to test whether
// bijection extraction failed (around the square root of the leaf size).
template<int N>
static constexpr array<uint8_t, N> fill_bij_midstop() {
    array<uint8_t, N> memo{ 0 };
    for (int s = 0; s < N; ++s) memo[s] = s < (int)ce::ceil(MIDSTOP_FACTOR * ce::sqrt(s)) ?
        s : (int)ce::ceil(MIDSTOP_FACTOR * ce::sqrt(s));
    return memo;
}

// Is zero beyond 2^31 (if N > 31)
template<size_t N>
static constexpr array<uint32_t, N> fill_powers_of_two() {
    array<uint32_t, N> memo{ 0 };
    for (size_t s = 0; s < min(N, 32); ++s) memo[s] = 1 << s;
    return memo;
}

// The first 4 entries are 0, then there are the first 4 powers of 256
// and then a zero padding.
static constexpr array<uint32_t, 4 + MAX_FANOUT + FULL_VEC_32_COUNT> fill_aggr_level_count_lookup() {
    array<uint32_t, 4 + MAX_FANOUT + FULL_VEC_32_COUNT> memo{ 0 };
    for (size_t s = 0; s < 4; ++s) memo[s + 4] = 1 << (8 * s);
    return memo;
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

template <size_t LEAF_SIZE, util::AllocType AT = util::AllocType::MALLOC, bool USE_BIJECTIONS_ROTATE = true>
class SIMDRecSplit
    : public AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, false> {
    using Superclass = AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, false>;
    using Superclass::SplitStrat;
    using Superclass::_leaf;
    using Superclass::lower_aggr;
    using Superclass::upper_aggr;
    using Superclass::memo;
    using Superclass::use_bijections_rotate;

    static constexpr auto bij_midstop = fill_bij_midstop<LEAF_SIZE + 1>();
    static constexpr auto powers_of_two = fill_powers_of_two<LEAF_SIZE + FULL_VEC_32_COUNT>(); // padding for lookup
    static constexpr auto aggr_level_count_lookup = fill_aggr_level_count_lookup();

  public:
    SIMDRecSplit() {}

    /** Builds a SIMDRecSplit instance using a given list of keys and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param keys a vector of strings.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    SIMDRecSplit(const vector<string> &keys, const size_t bucket_size, int num_threads) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash128_t *h = (hash128_t *)malloc(this->keys_count * sizeof(hash128_t));

        if (num_threads == 1) {
            for (size_t i = 0; i < this->keys_count; ++i) {
                h[i] = first_hash(keys[i].c_str(), keys[i].size());
            }
        } else {
            size_t keysPerThread = this->keys_count / num_threads + 1;
            std::vector<std::thread> threads;
            for (size_t thread = 0; thread < num_threads; thread++) {
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

        hash_gen(h, num_threads);
        free(h);
    }

    /** Builds a SIMDRecSplit instance using a given list of 128-bit hashes and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * Note that this constructor is mainly useful for benchmarking.
     * @param keys a vector of 128-bit hashes.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    SIMDRecSplit(vector<hash128_t> &keys, const size_t bucket_size, int num_threads) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash_gen(&keys[0], num_threads);
    }

    /** Builds a SIMDRecSplit instance using a list of keys returned by a stream and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param input an open input stream returning a list of keys, one per line.
     * @param bucket_size the desired bucket size.
     */
    SIMDRecSplit(ifstream& input, const size_t bucket_size, int num_threads) {
        this->bucket_size = bucket_size;
        vector<hash128_t> h;
        for(string key; getline(input, key);) h.push_back(first_hash(key.c_str(), key.size()));
        this->keys_count = h.size();
        hash_gen(&h[0], num_threads);
    }

  private:
    uint64_t higherLevel(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m, size_t split, typename RiceBitVector<AT>::Builder &builder,
        vector<uint32_t> &unary, const int level) {
#ifdef MORESTATS
        auto start_time = high_resolution_clock::now();
#endif
        const auto end = start + m;
        SIMDRS_ASSUME(m > upper_aggr);
        SIMDRS_ASSUME(m <= MAX_BUCKET_SIZE);
        uint64_t x = start_seed[level];
#ifdef SIMDRS_512_BIT
        FullVecUq xVec(0, 1, 2, 3, 4, 5, 6, 7);
#else
        FullVecUq xVec(0, 1, 2, 3);
#endif
        xVec += x;

        FullVecUi counter;
        FullVecIb found_result;
        for (;;) {
            counter = 0;
            for (size_t i = start; i < end; i++) {
                const FullVecUq first = bucket[i] + xVec;
                // seems equal to "counter += remap(first, first + FULL_VEC_64_COUNT, m) < split" and faster than if_sub
                counter = if_add(remap(first, first + FULL_VEC_64_COUNT, m) < split, counter, uint32_t(-1));
#ifdef MORESTATS
                ++num_split_evals;
#endif
            }
            found_result = counter == -split; // -split since true is represented as -1 in vectors
            if (horizontal_or(found_result)) break;
            x += FULL_VEC_32_COUNT;
            xVec += FULL_VEC_32_COUNT;
        }

        const auto found_idx = horizontal_find_first(found_result);
        SIMDRS_ASSUME(found_idx != -1);
        x += found_idx;

        size_t count[2];
        count[0] = 0;
        count[1] = split;
        size_t i;
        for (i = start; i + FULL_VEC_32_COUNT <= end; i += FULL_VEC_32_COUNT) {
            FullVecUq first, second;
            first.load(&bucket[i]);
            second.load(&bucket[i + FULL_VEC_64_COUNT]);
            auto bits = to_bits(remap(first + x, second + x, m) >= split); // Fast for AVX-512
            for (size_t j = 0; j < FULL_VEC_32_COUNT; ++j)
                temp[count[(bits >> j) & 1]++] = bucket[i + j]; // TODO: Vectorize? Probably hard!
        }
        FullVecUq first, second;
        first.load(&bucket[i]);
        second.load(&bucket[i + FULL_VEC_64_COUNT]);
        auto bits = to_bits(remap(first + x, second + x, m) >= split);
        for (size_t j = 0; j < end - i; ++j)
            temp[count[(bits >> j) & 1]++] = bucket[i + j];
        copy(&temp[0], &(temp.data()[m]), &bucket[start]);

        x -= start_seed[level];
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);

#ifdef MORESTATS
        time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        return x;
    }

    template<uint32_t _MAX_FANOUT, uint32_t SPLIT, uint64_t SEED>
    uint64_t aggrLevel(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m, typename RiceBitVector<AT>::Builder &builder,
        vector<uint32_t> &unary, [[maybe_unused]] const int level) {
#ifdef MORESTATS
        auto start_time = high_resolution_clock::now();
#endif
        const auto end = start + m;
        uint64_t x = SEED;
        const uint32_t fanout = uint16_t(m + SPLIT - 1) / SPLIT;
        SIMDRS_ASSUME(m > LEAF_SIZE);
        SIMDRS_ASSUME(m <= _MAX_FANOUT * SPLIT);
        SIMDRS_ASSUME(fanout >= 2);
        SIMDRS_ASSUME(fanout <= _MAX_FANOUT);
        size_t i;
        static_assert(_MAX_FANOUT <= MAX_FANOUT, "_MAX_FANOUT must be at most MAX_FANOUT!");
        static_assert(SPLIT < 255, "SPLIT must be less than 255 for aggrLevel to work correctly!"
            "Note that less than 256 is not enough since an overflow may carry to the next count.");
#ifdef SIMDRS_512_BIT
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3, SEED + 4, SEED + 5, SEED + 6, SEED + 7);
#else
        FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3);
#endif
        FullVecIb found_result;
        if (fanout <= 5) {
            FullVecUi count;
            uint32_t found = SPLIT;
            for (i = 1; i < fanout - 1; ++i)
                found |= SPLIT << (8 * i);
            if (fanout <= 4)
                found |= (m - (fanout - 1) * SPLIT) << (8 * i);
            for (;;) {
                count = 0;
                for (size_t i = start; i < end; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m) / const_uint(SPLIT);
                    count += lookup<min(_MAX_FANOUT, 5)>(remapped, &aggr_level_count_lookup[4]);
#ifdef MORESTATS
                    ++num_split_evals;
#endif
                }
                found_result = count == found;
                if (horizontal_or(found_result)) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
        } else {
            FullVecUi count_low;
            FullVecUi count_high;
            constexpr uint32_t found_low = (SPLIT << 24) | (SPLIT << 16) | (SPLIT << 8) | SPLIT;
            uint32_t found_high = SPLIT;
            for (i = 1; i < fanout - 5; ++i)
                found_high |= SPLIT << (8 * i);
            if (fanout <= 8)
                found_high |= (m - (fanout - 1) * SPLIT) << (8 * i);
            for (;;) {
                count_low = count_high = 0;
                for (size_t i = start; i < end; i++) {
                    const FullVecUq first = bucket[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m) / const_uint(SPLIT);
                    count_low += lookup<_MAX_FANOUT>(remapped, &aggr_level_count_lookup[4]);
                    count_high += lookup<_MAX_FANOUT>(remapped, &aggr_level_count_lookup[0]);
#ifdef MORESTATS
                    ++num_split_evals;
#endif
                }
                found_result = count_low == found_low & count_high == found_high;
                if (horizontal_or(found_result)) break;
                x += FULL_VEC_32_COUNT;
                xVec += FULL_VEC_32_COUNT;
            }
        }

        const auto found_idx = horizontal_find_first(found_result);
        SIMDRS_ASSUME(found_idx != -1);
        x += found_idx;

        uint64_t *SIMDRS_RESTRICT temp_c = &temp[m];
        for (size_t i = 0, c = 0; i < fanout; i++, c += SPLIT) temp_c[i] = c;
        uint32_t remapped[FULL_VEC_32_COUNT];
        for (i = start; i + FULL_VEC_32_COUNT <= end; i += FULL_VEC_32_COUNT) {
            FullVecUq first, second;
            first.load(&bucket[i]);
            second.load(&bucket[i + FULL_VEC_64_COUNT]);
            (remap(first + x, second + x, m) / const_uint(SPLIT)).store(remapped);
            for (size_t j = 0; j < FULL_VEC_32_COUNT; ++j)
                temp[temp_c[remapped[j]]++] = bucket[i + j]; // TODO: Vectorize? Probably hard!
        }
        FullVecUq first, second;
        first.load(&bucket[i]);
        second.load(&bucket[i + FULL_VEC_64_COUNT]);
        (remap(first + x, second + x, m) / const_uint(SPLIT)).store(remapped);
        for (size_t j = 0; j < end - i; ++j)
            temp[temp_c[remapped[j]]++] = bucket[i + j];
        copy(&temp[0], &(temp.data()[m]), &bucket[start]);

        x -= SEED;
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);

#ifdef MORESTATS
        time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        return x;
    }

    static FullVecUi powerOfTwo(FullVecUi x) {
        if constexpr (NO_AVX) {
            return FullVecUi(lookup<LEAF_SIZE>(x, &powers_of_two[0]));
        } else {
#ifdef SIMDRS_512_BIT
            return _mm512_sllv_epi32(FullVecUi(1), x);
#else
            return _mm256_sllv_epi32(FullVecUi(1), x);
#endif
        }
    }

    static FullVecUi rotateRight(FullVecUi val, uint32_t x) {
        return ((val << uint32_t(LEAF_SIZE - x)) | (val >> x)) & ((1 << LEAF_SIZE) - 1);
    }

    void leafLevel(vector<uint64_t> &bucket, size_t start, size_t m, typename RiceBitVector<AT>::Builder &builder,
            vector<uint32_t> &unary, [[maybe_unused]] const int level) {
        SIMDRS_ASSUME(m >= 2);
        SIMDRS_ASSUME(m <= LEAF_SIZE);
        const auto end = start + m;
        constexpr uint64_t SEED = start_seed[NUM_START_SEEDS - 1];
        uint64_t x = SEED;
#ifdef MORESTATS
        sum_depths += m * level;
        auto start_time = high_resolution_clock::now();
#endif
        if (use_bijections_rotate && m == LEAF_SIZE) {
            int items_left_count = 0;
            int items_right_count = 0;
            uint64_t items_left[LEAF_SIZE];
            uint64_t items_right[LEAF_SIZE];
            for (size_t i = start; i < end; i++) {
                if (bucket[i] % 2 == 0) {
                    items_left[items_left_count] = bucket[i];
                    items_left_count++;
                } else {
                    items_right[items_right_count] = bucket[i];
                    items_right_count++;
                }
            }

            FullVecUi mask_left, mask_right;
            constexpr uint32_t found = uint32_t(1 << LEAF_SIZE) - 1;
#ifdef SIMDRS_512_BIT
            FullVecUq xVec(SEED, SEED + 1 * LEAF_SIZE, SEED + 2 * LEAF_SIZE, SEED + 3 * LEAF_SIZE,
                SEED + 4 * LEAF_SIZE, SEED + 5 * LEAF_SIZE, SEED + 6 * LEAF_SIZE, SEED + 7 * LEAF_SIZE);
#else
            FullVecUq xVec(SEED, SEED + 1 * LEAF_SIZE, SEED + 2 * LEAF_SIZE, SEED + 3 * LEAF_SIZE);
#endif
            for (;; x += FULL_VEC_32_COUNT * LEAF_SIZE, xVec += FULL_VEC_32_COUNT * LEAF_SIZE) {
                mask_left = mask_right = 0;
                for (int i = 0; i < items_left_count; i++) {
                    const FullVecUq first = items_left[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT * LEAF_SIZE, LEAF_SIZE);
                    mask_left |= powerOfTwo(remapped);
                }
                for (int i = 0; i < items_right_count; i++) {
                    const FullVecUq first = items_right[i] + xVec;
                    const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT * LEAF_SIZE, LEAF_SIZE);
                    mask_right |= powerOfTwo(remapped);
                }
#ifdef SIMDRS_512_BIT_POPCNT
                if constexpr (LEAF_SIZE >= 11) {
                    if (horizontal_and((FullVecUi(_mm512_popcnt_epi32(mask_right)) != items_right_count)
                                     | (FullVecUi(_mm512_popcnt_epi32(mask_left)) != items_left_count))) {
                        continue; // Collisions in either part
                    }
                }
#endif
                // Try to rotate right part to see if both together form a bijection
                size_t rotations;
                FullVecUi result = LEAF_SIZE;
                // We start with the highest rotation since we want to find the lowest which satisfies the condition.
                // The condition is false after rotations == 0 because of the underflow.
                for (rotations = LEAF_SIZE - 1; rotations < LEAF_SIZE; rotations--) {
                    mask_right = rotateRight(mask_right, 1);
                    result = select((mask_left | mask_right) == found, rotations, result);
                }
                const auto found_idx = horizontal_find_first(result < LEAF_SIZE);
                if (found_idx != -1) {
                    x += found_idx * LEAF_SIZE + result[found_idx];
                    break;
                }
            }
        } else {
            FullVecUi mask;
            const uint32_t found = (1 << m) - 1;
#ifdef SIMDRS_512_BIT
            FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3, SEED + 4, SEED + 5, SEED + 6, SEED + 7);
            constexpr size_t midstop_threshold = 12;
#else
            FullVecUq xVec(SEED, SEED + 1, SEED + 2, SEED + 3);
            constexpr size_t midstop_threshold = 11;
#endif
            if (m < midstop_threshold) {
                for (;;) {
                    mask = 0;
                    for (size_t i = start; i < end; i++) {
                        const FullVecUq first = bucket[i] + xVec;
                        const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m);
                        mask |= powerOfTwo(remapped);
                    }
#ifdef MORESTATS
                    num_bij_evals[m] += m;
#endif
                    if (horizontal_or(mask == found)) break;
                    x += FULL_VEC_32_COUNT;
                    xVec += FULL_VEC_32_COUNT;
                }
                const auto found_idx = horizontal_find_first(mask == found);
                SIMDRS_ASSUME(found_idx != -1);
                x += found_idx;
            } else {
                const uint32_t midstop = bij_midstop[m];
                uint32_t backlog_masks[FULL_VEC_32_COUNT];
                uint64_t backlog_x[FULL_VEC_32_COUNT];
                int backlog_size = 0;
                for (;;) {
                    mask = 0;
                    size_t i;
                    for (i = start; i < start + midstop; i++) {
                        const FullVecUq first = bucket[i] + xVec;
                        const FullVecUi remapped = remap(first, first + FULL_VEC_64_COUNT, m);
                        mask |= powerOfTwo(remapped);
                    }
#ifdef MORESTATS
                    num_bij_evals[m] += midstop;
#endif
                    bool brk = false;
#ifdef SIMDRS_512_BIT_POPCNT
                    FullVecUi popcounts = _mm512_popcnt_epi32(mask);
                    uint32_t bits = to_bits(popcounts == midstop);
                    while (bits != 0) {
                        uint32_t j = __builtin_ctz(bits); // We don't care for MSVC here since it doesn't define __AVX512VPOPCNTDQ__
                        bits = clear_rho(bits);
                        uint32_t interm_mask = mask[j];
#else
                    uint32_t intermediate[FULL_VEC_32_COUNT];
                    mask.store(intermediate);
                    for (uint32_t j = 0; j < FULL_VEC_32_COUNT; ++j) {
                        uint32_t interm_mask = intermediate[j];
                        if (nu(interm_mask) == (int)midstop) {
#endif
                            backlog_masks[backlog_size] = interm_mask;
                            backlog_x[backlog_size] = x + j;
                            if (++backlog_size == FULL_VEC_32_COUNT) {
                                FullVecUq first, second;
                                first.load(backlog_x);
                                second.load(backlog_x + FULL_VEC_64_COUNT);
                                FullVecUi backlog_mask;
                                backlog_mask.load(backlog_masks);
                                for (; i < end; i++) {
                                    const FullVecUi remapped = remap(first + bucket[i], second + bucket[i], m);
                                    backlog_mask |= powerOfTwo(remapped);
                                }
                                if (horizontal_or(backlog_mask == found)) {
                                    brk = true;
                                    const auto found_idx = horizontal_find_first(backlog_mask == found);
                                    SIMDRS_ASSUME(found_idx != -1);
                                    x = backlog_x[found_idx];
                                    break;
                                }
                                backlog_size = 0;
                            }
#ifndef SIMDRS_512_BIT_POPCNT
                        }
#endif
                    }
                    if (brk) break;
                    x += FULL_VEC_32_COUNT;
                    xVec += FULL_VEC_32_COUNT;
                }
            }
        }
#ifdef MORESTATS
        time_bij += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
        x -= SEED;
        const auto log2golomb = golomb_param(m);
        builder.appendFixed(x, log2golomb);
        unary.push_back(x >> log2golomb);
#ifdef MORESTATS
        bij_count[m]++;
        num_bij_trials[m] += x + 1;
        bij_unary += 1 + (x >> log2golomb);
        bij_fixed += log2golomb;

        min_bij_code = min(min_bij_code, x);
        max_bij_code = max(max_bij_code, x);
        sum_bij_codes += x;

        auto b = bij_memo_golomb[m];
        auto log2b = lambda(b);
        bij_unary_golomb += x / b + 1;
        bij_fixed_golomb += x % b < ((1 << (log2b + 1)) - b) ? log2b : log2b + 1;
#endif
    }

    void recSplit(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t m,
            typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, const int level = 0) {
        assert(m > 1);
        if (m <= _leaf) {
            leafLevel(bucket, start, m, builder, unary, level);
        } else {
            [[maybe_unused]] uint64_t x;
            if (m > upper_aggr) { // fanout = 2
                const size_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
                x = higherLevel(bucket, temp, start, m, split, builder, unary, level);
                recSplit(bucket, temp, start, split, builder, unary, level + 1);
                if (m - split > 1) recSplit(bucket, temp, start + split, m - split, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else if (m > lower_aggr) { // 2nd aggregation level
                x = aggrLevel<upper_aggr / lower_aggr, lower_aggr, start_seed[NUM_START_SEEDS - 3]>(bucket, temp, start, m, builder, unary, level);
                size_t i;
                for (i = 0; i < m - lower_aggr; i += lower_aggr) {
                    recSplit(bucket, temp, start + i, lower_aggr, builder, unary, level + 1);
                }
                if (m - i > 1) recSplit(bucket, temp, start + i, m - i, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else { // First aggregation level, m <= lower_aggr
                x = aggrLevel<lower_aggr / _leaf, _leaf, start_seed[NUM_START_SEEDS - 2]>(bucket, temp, start, m, builder, unary, level);
                size_t i;
                for (i = 0; i < m - _leaf; i += _leaf) {
                    leafLevel(bucket, start + i, _leaf, builder, unary, level + 1);
                }
                if (m - i > 1) leafLevel(bucket, start + i, m - i, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            }

#ifdef MORESTATS
            ++split_count;
            num_split_trials += x + 1;
            double e_trials = 1;
            size_t aux = m;
            SplittingStrategy<LEAF_SIZE> strat(m);
            auto v = strat.begin();
            for (int i = 0; i < strat.fanout(); ++i, ++v) {
                e_trials *= pow((double)m / *v, *v);
                for (size_t j = *v; j > 0; --j, --aux) {
                    e_trials *= (double)j / aux;
                }
            }
            expected_split_trials += (size_t)e_trials;
            expected_split_evals += (size_t)e_trials * m;
            const auto log2golomb = golomb_param(m);
            split_unary += 1 + (x >> log2golomb);
            split_fixed += log2golomb;

            min_split_code = min(min_split_code, x);
            max_split_code = max(max_split_code, x);
            sum_split_codes += x;

            auto b = split_golomb_b<LEAF_SIZE>(m);
            auto log2b = lambda(b);
            split_unary_golomb += x / b + 1;
            split_fixed_golomb += x % b < ((1ULL << (log2b + 1)) - b) ? log2b : log2b + 1;
#endif
        }
    }

    void compute_thread(int tid, int num_threads, mutex &mtx, condition_variable &condition,
                        vector<uint64_t> &bucket_size_acc, vector<uint64_t> &bucket_pos_acc,
                        vector<uint64_t> &sorted_keys, int &next_thread_to_append_builder,
                        typename RiceBitVector<AT>::Builder &builder) {
        typename RiceBitVector<AT>::Builder local_builder;
        vector<uint32_t> unary;
        vector<uint64_t> temp(MAX_BUCKET_SIZE);
        size_t begin = tid * this->nbuckets / num_threads;
        size_t end = std::min(this->nbuckets, (tid + 1) * this->nbuckets / num_threads);
        if (tid == num_threads - 1) {
            end = this->nbuckets;
        }
        for (size_t i = begin; i < end; ++i) {
            const size_t s = bucket_size_acc[i + 1] - bucket_size_acc[i];
            if (s > 1) {
                recSplit(sorted_keys, temp, bucket_size_acc[i], s, local_builder, unary);
                local_builder.appendUnaryAll(unary);
                unary.clear();
            }
            bucket_pos_acc[i + 1] = local_builder.getBits();
#ifdef MORESTATS
            auto upper_leaves = (s + _leaf - 1) / _leaf;
            auto upper_height = ceil(log(upper_leaves) / log(2)); // TODO: check
            auto upper_s = _leaf * pow(2, upper_height);
            ub_split_bits += (double)upper_s / (_leaf * 2) * log2(2 * M_PI * _leaf) - .5 * log2(2 * M_PI * upper_s);
            ub_bij_bits += upper_leaves * _leaf * (log2e - .5 / _leaf * log2(2 * M_PI * _leaf));
            ub_split_evals += 4 * upper_s * sqrt(pow(2 * M_PI * upper_s, 2 - 1) / pow(2, 2));
            minsize = min(minsize, s);
            maxsize = max(maxsize, s);
#endif
        }
        if (tid == 0) {
            builder = std::move(local_builder);
            lock_guard<mutex> lock(mtx);
            next_thread_to_append_builder = 1;
            condition.notify_all();
        } else {
            uint64_t prev_bucket_pos;
            {
                unique_lock<mutex> lock(mtx);
                condition.wait(lock, [&] { return next_thread_to_append_builder == tid; });
                prev_bucket_pos = builder.getBits();
                builder.appendRiceBitVector(local_builder);
                next_thread_to_append_builder = tid + 1;
                condition.notify_all();
            }
            for (size_t i = begin + 1; i < end + 1; ++i) {
                bucket_pos_acc[i] += prev_bucket_pos;
            }
        }
    }

    void hash_gen(hash128_t *hashes, int num_threads) {
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
        minsize = this->keys_count;
        maxsize = 0;
        ub_split_bits = 0;
        ub_bij_bits = 0;
        ub_split_evals = 0;

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
        auto sorted_keys = vector<uint64_t>(this->keys_count + FULL_VEC_32_COUNT); // Padding for vectorized access

#ifdef MORESTATS
        auto sorting_start_time = high_resolution_clock::now();
#endif
        this->parallelPartition(hashes, sorted_keys, bucket_size_acc, num_threads);
#ifdef MORESTATS
        auto sorting_time = duration_cast<nanoseconds>(high_resolution_clock::now() - sorting_start_time).count();
#endif
        typename RiceBitVector<AT>::Builder builder;

        vector<thread> threads;
        threads.reserve(num_threads);
        mutex mtx;
        condition_variable condition;
        int next_thread_to_append_builder = 0;
        bucket_pos_acc[0] = 0;
        if (num_threads == 1) {
             compute_thread(0, num_threads, mtx, condition,
                   bucket_size_acc, bucket_pos_acc, sorted_keys,
                   next_thread_to_append_builder, builder);
        } else {
            for (int tid = 0; tid < num_threads; ++tid) {
                threads.emplace_back([&, tid] {
                    compute_thread(tid, num_threads, mtx, condition,
                                   bucket_size_acc, bucket_pos_acc, sorted_keys,
                                   next_thread_to_append_builder, builder);
                });
            }
            for (auto &thread: threads) {
                thread.join();
            }
        }
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
        printf("Total time: %13.3f ms\n", duration_cast<nanoseconds>(high_resolution_clock::now() - total_start_time).count() * 1E-6);
        printf("Sorting time: %11.3f ms\n", sorting_time * 1E-6);
        printf("Bijections: %13.3f ms\n", time_bij * 1E-6);
        for (int i = 0; i < MAX_LEVEL_TIME; i++) {
            if (time_split[i] > 0) {
                printf("Split level %d: %10.3f ms\n", i, time_split[i] * 1E-6);
            }
        }

        uint64_t fact = 1, tot_bij_count = 0, tot_bij_evals = 0;
        printf("\n");
        printf("Bij               count              trials                 exp               evals                 exp           tot evals\n");
        for (int i = 0; i < MAX_LEAF_SIZE; i++) {
            if (num_bij_trials[i] != 0) {
                tot_bij_count += bij_count[i];
                tot_bij_evals += num_bij_evals[i];
                printf("%-3d%20d%20.2f%20.2f%20.2f%20.2f%20lld\n", i, bij_count[i], (double)num_bij_trials[i] / bij_count[i], pow(i, i) / fact, (double)num_bij_evals[i] / bij_count[i],
                       (_leaf <= 8 ? i : bij_midstop[i]) * pow(i, i) / fact, num_bij_evals[i]);
            }
            fact *= (i + 1);
        }

        printf("\n");
        printf("Split count:       %16zu\n", split_count);

        printf("Total split evals: %16lld\n", num_split_evals);
        printf("Total bij evals:   %16lld\n", tot_bij_evals);
        printf("Total evals:       %16lld\n", num_split_evals + tot_bij_evals);

        printf("\n");
        printf("Average depth:        %f\n", (double)sum_depths / this->keys_count);
        printf("\n");
        printf("Trials per split:     %16.3f\n", (double)num_split_trials / split_count);
        printf("Exp trials per split: %16.3f\n", (double)expected_split_trials / split_count);
        printf("Evals per split:      %16.3f\n", (double)num_split_evals / split_count);
        printf("Exp evals per split:  %16.3f\n", (double)expected_split_evals / split_count);

        printf("\n");
        printf("Unary bits per bij: %10.5f\n", (double)bij_unary / tot_bij_count);
        printf("Fixed bits per bij: %10.5f\n", (double)bij_fixed / tot_bij_count);
        printf("Total bits per bij: %10.5f\n", (double)(bij_unary + bij_fixed) / tot_bij_count);

        printf("\n");
        printf("Unary bits per split: %10.5f\n", (double)split_unary / split_count);
        printf("Fixed bits per split: %10.5f\n", (double)split_fixed / split_count);
        printf("Total bits per split: %10.5f\n", (double)(split_unary + split_fixed) / split_count);
        printf("Total bits per key:   %10.5f\n", (double)(bij_unary + bij_fixed + split_unary + split_fixed) / this->keys_count);

        printf("\n");
        printf("Unary bits per bij (Golomb): %10.5f\n", (double)bij_unary_golomb / tot_bij_count);
        printf("Fixed bits per bij (Golomb): %10.5f\n", (double)bij_fixed_golomb / tot_bij_count);
        printf("Total bits per bij (Golomb): %10.5f\n", (double)(bij_unary_golomb + bij_fixed_golomb) / tot_bij_count);

        printf("\n");
        printf("Unary bits per split (Golomb): %10.5f\n", (double)split_unary_golomb / split_count);
        printf("Fixed bits per split (Golomb): %10.5f\n", (double)split_fixed_golomb / split_count);
        printf("Total bits per split (Golomb): %10.5f\n", (double)(split_unary_golomb + split_fixed_golomb) / split_count);
        printf("Total bits per key (Golomb):   %10.5f\n", (double)(bij_unary_golomb + bij_fixed_golomb + split_unary_golomb + split_fixed_golomb) / this->keys_count);

        printf("\n");

        printf("Total split bits        %16.3f\n", (double)split_fixed + split_unary);
        printf("Upper bound split bits: %16.3f\n", ub_split_bits);
        printf("Total bij bits:         %16.3f\n", (double)bij_fixed + bij_unary);
        printf("Upper bound bij bits:   %16.3f\n\n", ub_bij_bits);
#endif
    }
};

} // namespace sux::function
