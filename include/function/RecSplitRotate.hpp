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
 * Should work with GCC and MSVC (GCC is recommended). They may produce different
 * MPHFs due to the implementation of remap128.
 */

#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <util/Hash128.h>
#include <thread>
#include <condition_variable>
#include <sux/util/Vector.hpp>
#include "DoubleEF.hpp"
#include "RiceBitVector.hpp"
#include "golombRiceMemo.hpp"
#include "AbstractParallelRecSplit.hpp"

namespace bez::function::recsplit_rotate {

using namespace std;
using namespace std::chrono;

static const int MAX_LEAF_SIZE = 24;
static const int MAX_FANOUT = 32;

#if defined(MORESTATS) && !defined(STATS)
#define STATS
#endif

#ifdef MORESTATS

#define MAX_LEVEL_TIME (20)

static constexpr double log2e = 1.44269504089;
static uint64_t num_bij_trials[MAX_LEAF_SIZE], num_split_trials;
static uint64_t num_bij_evals[MAX_LEAF_SIZE], num_split_evals;
static uint64_t bij_count[MAX_LEAF_SIZE], split_count;
static uint64_t expected_split_trials, expected_split_evals;
static uint64_t bij_unary, bij_fixed, bij_unary_golomb, bij_fixed_golomb;
static uint64_t split_unary, split_fixed, split_unary_golomb, split_fixed_golomb;
static uint64_t max_split_code, min_split_code, sum_split_codes;
static uint64_t max_bij_code, min_bij_code, sum_bij_codes;
static uint64_t sum_depths;
static uint64_t time_bij;
static uint64_t time_split[MAX_LEVEL_TIME];
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846 // This is somehow necessary with MSVC
#endif

#ifdef MORESTATS
// Optimal Golomb code moduli for leaves (for stats).
static constexpr uint64_t bij_memo_golomb[] = {0,        0,        1,         3,         7,          18,         45,          113,         288,        740,
                                               1910,     4954,     12902,     33714,     88350,      232110,     611118,      1612087,     4259803,    11273253,
                                               29874507, 79265963, 210551258, 559849470, 1490011429, 3968988882, 10580669970, 28226919646, 75354118356};
#endif

// Computes the point at which one should stop to test whether
// bijection extraction failed (around the square root of the leaf size).

static constexpr array<uint8_t, MAX_LEAF_SIZE> fill_bij_midstop() {
    array<uint8_t, MAX_LEAF_SIZE> memo{0};
    for (int s = 0; s < MAX_LEAF_SIZE; ++s) memo[s] = s < (int)ce::ceil(2 * ce::sqrt(s)) ? s : (int)ce::ceil(2 * ce::sqrt(s));
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

template <size_t LEAF_SIZE, sux::util::AllocType AT = sux::util::AllocType::MALLOC, bool USE_BIJECTIONS_ROTATE = true>
class RecSplit : public AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, true> {
    using Superclass = AbstractParallelRecSplit<LEAF_SIZE, AT, USE_BIJECTIONS_ROTATE, true>;
    using Superclass::SplitStrat;
    using Superclass::_leaf;
    using Superclass::lower_aggr;
    using Superclass::upper_aggr;
    using Superclass::memo;
    using Superclass::use_bijections_rotate;
    using Superclass::remixAndRemap;

    static constexpr array<uint8_t, MAX_LEAF_SIZE> bij_midstop = fill_bij_midstop();

    using Superclass::ef;
    using Superclass::descriptors;
    using Superclass::nbuckets;
    using Superclass::keys_count;
    using Superclass::bucket_size;
    using Superclass::hash128_to_bucket;
  public:
    RecSplit() {}

    /** Builds a RecSplit instance using a given list of keys and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param keys a vector of strings.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    RecSplit(const vector<string> &keys, const size_t bucket_size, const size_t num_threads) {
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

    /** Builds a RecSplit instance using a given list of 128-bit hashes and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * Note that this constructor is mainly useful for benchmarking.
     * @param keys a vector of 128-bit hashes.
     * @param bucket_size the desired bucket size; typical sizes go from
     * 100 to 2000, with smaller buckets giving slightly larger but faster
     * functions.
     */
    RecSplit(vector<hash128_t> &keys, const size_t bucket_size, const size_t num_threads) {
        this->bucket_size = bucket_size;
        this->keys_count = keys.size();
        hash_gen(&keys[0], num_threads);
    }

    /** Builds a RecSplit instance using a list of keys returned by a stream and bucket size.
     *
     * **Warning**: duplicate keys will cause this method to never return.
     *
     * @param input an open input stream returning a list of keys, one per line.
     * @param bucket_size the desired bucket size.
     */
    RecSplit(ifstream& input, const size_t bucket_size, const size_t num_threads) {
        this->bucket_size = bucket_size;
        vector<hash128_t> h;
        for(string key; getline(input, key);) h.push_back(first_hash(key.c_str(), key.size()));
        this->keys_count = h.size();
        hash_gen(&h[0], num_threads);
    }

    /** Returns the number of keys used to build this RecSplit instance. */
    inline size_t size() { return this->keys_count; }

  private:
    static constexpr uint32_t rotate(uint32_t val, uint8_t x) {
        return ((val << x) | (val >> (LEAF_SIZE - x))) & ((1 << LEAF_SIZE) - 1);
    }

    void recSplit(vector<uint64_t> &bucket, vector<uint64_t> &temp, size_t start, size_t end,
                  typename RiceBitVector<AT>::Builder &builder, vector<uint32_t> &unary, const int level = 0) {
        const auto m = end - start;
        assert(m > 1);
        uint64_t x;

        if (m <= _leaf) {
            x = start_seed[NUM_START_SEEDS - 1];
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

                constexpr uint32_t found = uint32_t(1 << LEAF_SIZE) - 1;
                for (;; x += LEAF_SIZE) {
                    uint32_t mask_left = 0;
                    uint32_t mask_right = 0;
                    for (int i = 0; i < items_left_count; i++) {
                        mask_left |= uint32_t(1) << remixAndRemap(items_left[i] + x, LEAF_SIZE);
                    }
                    if (nu(mask_left) != items_left_count) {
                        continue; // Collisions in left part
                    }
                    for (int i = 0; i < items_right_count; i++) {
                        mask_right |= uint32_t(1) << remixAndRemap(items_right[i] + x, LEAF_SIZE);
                    }
                    if (nu(mask_right) != items_right_count) {
                        continue; // Collisions in right part
                    }
                    // Try to rotate right part to see if both together form a bijection
                    size_t rotations;
                    for (rotations = 0; rotations < LEAF_SIZE; rotations++) {
                        if ((mask_left | mask_right) == found) {
                            x += rotations;
                            break;
                        }
                        mask_right = rotate(mask_right, 1);
                    }
                    if (rotations < LEAF_SIZE)
                        break;
                }
            } else {
                uint32_t mask;
                const uint32_t found = (1 << m) - 1;
                if constexpr (_leaf <= 8) {
                    for (;;) {
                        mask = 0;
                        for (size_t i = start; i < end; i++) mask |= uint32_t(1) << remixAndRemap(bucket[i] + x, m);
#ifdef MORESTATS
                        num_bij_evals[m] += m;
#endif
                        if (mask == found) break;
                        x++;
                    }
                } else {
                    const int midstop = bij_midstop[m];
                    for (;;) {
                        mask = 0;
                        size_t i;
                        for (i = start; i < start + midstop; i++) mask |= uint32_t(1) << remixAndRemap(bucket[i] + x, m);
#ifdef MORESTATS
                        num_bij_evals[m] += midstop;
#endif
                        if (nu(mask) == midstop) {
                            for (; i < end; i++) mask |= uint32_t(1) << remixAndRemap(bucket[i] + x, m);
#ifdef MORESTATS
                            num_bij_evals[m] += m - midstop;
#endif
                            if (mask == found) break;
                        }
                        x++;
                    }
                }
            }
#ifdef MORESTATS
            time_bij += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
            x -= start_seed[NUM_START_SEEDS - 1];
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
        } else {
#ifdef MORESTATS
            auto start_time = high_resolution_clock::now();
#endif
            if (m > upper_aggr) { // fanout = 2
                x = start_seed[level];
                const size_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;

                size_t count[2];
                for (;;) {
                    count[0] = 0;
                    for (size_t i = start; i < end; i++) {
                        count[remixAndRemap(bucket[i] + x, m) >= split]++;
#ifdef MORESTATS
                        ++num_split_evals;
#endif
                    }
                    if (count[0] == split) break;
                    x++;
                }

                count[0] = 0;
                count[1] = split;
                for (size_t i = start; i < end; i++) {
                    temp[count[remixAndRemap(bucket[i] + x, m) >= split]++] = bucket[i];
                }
                copy(&temp[0], &(temp.data()[m]), &bucket[start]);
                x -= start_seed[level];
                const auto log2golomb = golomb_param(m);
                builder.appendFixed(x, log2golomb);
                unary.push_back(x >> log2golomb);

#ifdef MORESTATS
                time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
                recSplit(bucket, temp, start, start + split, builder, unary, level + 1);
                if (m - split > 1) recSplit(bucket, temp, start + split, end, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else if (m > lower_aggr) { // 2nd aggregation level
                x = start_seed[NUM_START_SEEDS - 3];
                const size_t fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
                size_t count[MAX_FANOUT];
                for (;;) {
                    memset(count, 0, sizeof count - sizeof *count);
                    for (size_t i = start; i < end; i++) {
                        count[uint16_t(remixAndRemap(bucket[i] + x, m)) / lower_aggr]++;
#ifdef MORESTATS
                        ++num_split_evals;
#endif
                    }
                    size_t broken = 0;
                    for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - lower_aggr;
                    if (!broken) break;
                    x++;
                }

                for (size_t i = 0, c = 0; i < fanout; i++, c += lower_aggr) count[i] = c;
                for (size_t i = start; i < end; i++) {
                    temp[count[uint16_t(remixAndRemap(bucket[i] + x, m)) / lower_aggr]++] = bucket[i];
                }
                copy(&temp[0], &(temp.data()[m]), &bucket[start]);

                x -= start_seed[NUM_START_SEEDS - 3];
                const auto log2golomb = golomb_param(m);
                builder.appendFixed(x, log2golomb);
                unary.push_back(x >> log2golomb);

#ifdef MORESTATS
                time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
                size_t i;
                for (i = 0; i < m - lower_aggr; i += lower_aggr) {
                    recSplit(bucket, temp, start + i, start + i + lower_aggr, builder, unary, level + 1);
                }
                if (m - i > 1) recSplit(bucket, temp, start + i, end, builder, unary, level + 1);
#ifdef MORESTATS
                else
                    sum_depths += level;
#endif
            } else { // First aggregation level, m <= lower_aggr
                x = start_seed[NUM_START_SEEDS - 2];
                const size_t fanout = uint16_t(m + _leaf - 1) / _leaf;
                size_t count[MAX_FANOUT];
                for (;;) {
                    memset(count, 0, sizeof count - sizeof *count);
                    for (size_t i = start; i < end; i++) {
                        count[uint16_t(remixAndRemap(bucket[i] + x, m)) / _leaf]++;
#ifdef MORESTATS
                        ++num_split_evals;
#endif
                    }
                    size_t broken = 0;
                    for (size_t i = 0; i < fanout - 1; i++) broken |= count[i] - _leaf;
                    if (!broken) break;
                    x++;
                }
                for (size_t i = 0, c = 0; i < fanout; i++, c += _leaf) count[i] = c;
                for (size_t i = start; i < end; i++) {
                    temp[count[uint16_t(remixAndRemap(bucket[i] + x, m)) / _leaf]++] = bucket[i];
                }
                copy(&temp[0], &(temp.data()[m]), &bucket[start]);

                x -= start_seed[NUM_START_SEEDS - 2];
                const auto log2golomb = golomb_param(m);
                builder.appendFixed(x, log2golomb);
                unary.push_back(x >> log2golomb);

#ifdef MORESTATS
                time_split[min(MAX_LEVEL_TIME, level)] += duration_cast<nanoseconds>(high_resolution_clock::now() - start_time).count();
#endif
                size_t i;
                for (i = 0; i < m - _leaf; i += _leaf) {
                    recSplit(bucket, temp, start + i, start + i + _leaf, builder, unary, level + 1);
                }
                if (m - i > 1) recSplit(bucket, temp, start + i, end, builder, unary, level + 1);
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
            SplittingStrategy<LEAF_SIZE> strat{m};
            auto v = strat.begin();
            for (size_t i = 0; i < strat.fanout(); ++i, ++v) {
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

    void compute_thread(int tid, int num_threads, mutex &mtx, std::condition_variable &condition,
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
                recSplit(sorted_keys, temp, bucket_size_acc[i], bucket_size_acc[i + 1], local_builder, unary);
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
        if (keys_count > (1ULL << 32)) {
            fprintf(stderr, "For more than 2^32 keys, you need 128-bit integer support.\n");
            abort();
        }
#endif
        nbuckets = max(1, (keys_count + bucket_size - 1) / bucket_size);
        auto bucket_size_acc = vector<uint64_t>(nbuckets + 1);
        auto bucket_pos_acc = vector<uint64_t>(nbuckets + 1);
        auto sorted_keys = vector<uint64_t>(keys_count);

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
        std::condition_variable condition;
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
        descriptors = builder.build();
        ef = DoubleEF<AT, true>(bucket_size_acc, bucket_pos_acc);

#ifdef STATS
        // Evaluation purposes only
        double ef_sizes = (double)ef.bitCountCumKeys() / keys_count;
        double ef_bits = (double)ef.bitCountPosition() / keys_count;
        double rice_desc = (double)builder.getBits() / keys_count;
        printf("Elias-Fano cumul sizes:  %f bits/bucket\n", (double)ef.bitCountCumKeys() / nbuckets);
        printf("Elias-Fano cumul bits:   %f bits/bucket\n", (double)ef.bitCountPosition() / nbuckets);
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
                printf("%-3d%20lu%20.2f%20.2f%20.2f%20.2f%20lu\n", i, bij_count[i], (double)num_bij_trials[i] / bij_count[i], pow(i, i) / fact, (double)num_bij_evals[i] / bij_count[i],
                       (_leaf <= 8 ? i : bij_midstop[i]) * pow(i, i) / fact, num_bij_evals[i]);
            }
            fact *= (i + 1);
        }

        printf("\n");
        printf("Split count:       %16zu\n", split_count);

        printf("Total split evals: %16lu\n", num_split_evals);
        printf("Total bij evals:   %16lu\n", tot_bij_evals);
        printf("Total evals:       %16lu\n", num_split_evals + tot_bij_evals);

        printf("\n");
        printf("Average depth:        %f\n", (double)sum_depths / keys_count);
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
        printf("Total bits per key:   %10.5f\n", (double)(bij_unary + bij_fixed + split_unary + split_fixed) / keys_count);

        printf("\n");
        printf("Unary bits per bij (Golomb): %10.5f\n", (double)bij_unary_golomb / tot_bij_count);
        printf("Fixed bits per bij (Golomb): %10.5f\n", (double)bij_fixed_golomb / tot_bij_count);
        printf("Total bits per bij (Golomb): %10.5f\n", (double)(bij_unary_golomb + bij_fixed_golomb) / tot_bij_count);

        printf("\n");
        printf("Unary bits per split (Golomb): %10.5f\n", (double)split_unary_golomb / split_count);
        printf("Fixed bits per split (Golomb): %10.5f\n", (double)split_fixed_golomb / split_count);
        printf("Total bits per split (Golomb): %10.5f\n", (double)(split_unary_golomb + split_fixed_golomb) / split_count);
        printf("Total bits per key (Golomb):   %10.5f\n", (double)(bij_unary_golomb + bij_fixed_golomb + split_unary_golomb + split_fixed_golomb) / keys_count);

        printf("\n");

        printf("Total split bits        %16.3f\n", (double)split_fixed + split_unary);
        printf("Upper bound split bits: %16.3f\n", ub_split_bits);
        printf("Total bij bits:         %16.3f\n", (double)bij_fixed + bij_unary);
        printf("Upper bound bij bits:   %16.3f\n\n", ub_bij_bits);
#endif
    }
};

} // namespace sux::function
