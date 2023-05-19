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
#include "DoubleEF.hpp"
#include "RiceBitVector.hpp"
#include "golombRiceMemo.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <util/Sorter.hpp>
#include <util/Hash128.h>

// Define constexpr namespace ce
#include <gcem.hpp>
namespace ce = gcem;

namespace bez::function {

using namespace std;

// Assumed *maximum* size of a bucket. Works with high probability up to average bucket size ~2000.
static constexpr int MAX_BUCKET_SIZE = 3000;

static constexpr int MAX_LEAF_SIZE = 24;
static constexpr int MAX_FANOUT = 9;

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
size_t minsize = 0, maxsize = 0;
double ub_split_bits = 0, ub_bij_bits = 0;
double ub_split_evals = 0;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846 // This is somehow necessary with MSVC
#endif

// Starting seed at given distance from the root (extracted at random).
static constexpr uint64_t start_seed[] = {0x106393c187cae21a, 0x6453cec3f7376937, 0x643e521ddbd2be98, 0x3740c6412f6572cb, 0x717d47562f1ce470, 0x4cd6eb4c63befb7c, 0x9bfd8c5e18c8da73,
                                          0x082f20e10092a9a3, 0x2ada2ce68d21defc, 0xe33cb4f3e7c6466b, 0x3980be458c509c59, 0xc466fd9584828e8c, 0x45f0aabe1a61ede6, 0xf6e7b8b33ad9b98d,
                                          0x4ef95e25f4b4983d, 0x81175195173b92d3, 0x4e50927d8dd15978, 0x1ea2099d1fafae7f, 0x425c8a06fbaaa815, 0xcd4216006c74052a};
static constexpr int NUM_START_SEEDS = sizeof(start_seed) / sizeof(uint64_t);

// Quick replacements for min/max on not-so-large integers.

static constexpr inline uint64_t min(int64_t x, int64_t y) { return y + ((x - y) & ((x - y) >> 63)); }
static constexpr inline uint64_t max(int64_t x, int64_t y) { return x - ((x - y) & ((x - y) >> 63)); }

// Optimal Golomb-Rice parameters for leaves.
static constexpr uint8_t bij_memo[] = {0, 0, 0, 1, 3, 4, 5, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 21, 22, 23, 25, 26, 28, 29, 30};

#ifdef MORESTATS
// Optimal Golomb code moduli for leaves (for stats).
static constexpr uint64_t bij_memo_golomb[] = {0,        0,        1,         3,         7,          18,         45,          113,         288,        740,
                                               1910,     4954,     12902,     33714,     88350,      232110,     611118,      1612087,     4259803,    11273253,
                                               29874507, 79265963, 210551258, 559849470, 1490011429, 3968988882, 10580669970, 28226919646, 75354118356};
#endif

/** A class emboding the splitting strategy of RecSplit.
 *
 *  Note that this class is used _for statistics only_. The splitting strategy is embedded
 *  into the generation code, which uses only the public fields SplittingStrategy::lower_aggr and SplittingStrategy::upper_aggr.
 */

template <size_t LEAF_SIZE> class SplittingStrategy {
    static constexpr size_t _leaf = LEAF_SIZE;
    static_assert(_leaf >= 2);
    static_assert(_leaf <= MAX_LEAF_SIZE);
    size_t m, curr_unit, curr_index, last_unit;
    size_t _fanout;
    size_t unit;

    inline size_t part_size() const { return (curr_index < _fanout - 1) ? unit : last_unit; }

  public:
    /** The lower bound for primary (lower) key aggregation. */
    static constexpr size_t lower_aggr = _leaf * max(2, ce::ceil(0.35 * _leaf + 1. / 2));
    /** The lower bound for secondary (upper) key aggregation. */
    static constexpr size_t upper_aggr = lower_aggr * (_leaf < 7 ? 2 : ce::ceil(0.21 * _leaf + 9. / 10));

    static inline constexpr void split_params(const size_t m, size_t &fanout, size_t &unit) {
        if (m > upper_aggr) { // High-level aggregation (fanout 2)
            unit = upper_aggr * (uint16_t(m / 2 + upper_aggr - 1) / upper_aggr);
            fanout = 2;
        } else if (m > lower_aggr) { // Second-level aggregation
            unit = lower_aggr;
            fanout = uint16_t(m + lower_aggr - 1) / lower_aggr;
        } else { // First-level aggregation
            unit = _leaf;
            fanout = uint16_t(m + _leaf - 1) / _leaf;
        }
    }

    // Note that you can call this iterator only *once*.
    class split_iterator {
        SplittingStrategy *strat;

      public:
        using value_type = size_t;
        using difference_type = ptrdiff_t;
        using pointer = size_t *;
        using reference = size_t &;
        using iterator_category = input_iterator_tag;

        split_iterator(SplittingStrategy *strat) : strat(strat) {}
        size_t operator*() const { return strat->curr_unit; }
        size_t *operator->() const { return &strat->curr_unit; }
        split_iterator &operator++() {
            ++strat->curr_index;
            strat->curr_unit = strat->part_size();
            strat->last_unit -= strat->curr_unit;
            return *this;
        }
        bool operator==(const split_iterator &other) const { return strat == other.strat; }
        bool operator!=(const split_iterator &other) const { return !(*this == other); }
    };

    explicit SplittingStrategy(size_t m) : m(m), last_unit(m), curr_index(0), curr_unit(0) {
        split_params(m, _fanout, unit);
        this->curr_unit = part_size();
        this->last_unit -= this->curr_unit;
    }

    split_iterator begin() { return split_iterator(this); }
    split_iterator end() { return split_iterator(nullptr); }

    inline size_t fanout() { return this->_fanout; }
};

// Generates the precomputed table of 32-bit values holding the Golomb-Rice code
// of a splitting (upper 5 bits), the number of nodes in the associated subtree
// (following 11 bits) and the sum of the Golomb-Rice codelengths in the same
// subtree (lower 16 bits).

template <size_t LEAF_SIZE> static constexpr void _fill_golomb_rice(const int m, array<uint32_t, MAX_BUCKET_SIZE> *memo) {
    array<int, MAX_FANOUT> k{0};

    size_t fanout = 0, unit = 0;
    SplittingStrategy<LEAF_SIZE>::split_params(m, fanout, unit);

    k[fanout - 1] = m;
    for (size_t i = 0; i < fanout - 1; ++i) {
        k[i] = unit;
        k[fanout - 1] -= k[i];
    }

    double sqrt_prod = 1;
    for (size_t i = 0; i < fanout; ++i) sqrt_prod *= ce::sqrt(k[i]);

    const double p = ce::sqrt(m) / (ce::pow(2 * M_PI, (fanout - 1.) / 2) * sqrt_prod);
    auto golomb_rice_length = (uint32_t)ce::ceil(ce::log2(-ce::log((ce::sqrt(5) + 1) / 2) / ce::log1p(-p))); // log2 Golomb modulus

    assert(golomb_rice_length <= 0x1F); // Golomb-Rice code, stored in the 5 upper bits
    (*memo)[m] = golomb_rice_length << 27;
    for (size_t i = 0; i < fanout; ++i) golomb_rice_length += (*memo)[k[i]] & 0xFFFF;
    assert(golomb_rice_length <= 0xFFFF); // Sum of Golomb-Rice codeslengths in the subtree, stored in the lower 16 bits
    (*memo)[m] |= golomb_rice_length;

    uint32_t nodes = 1;
    for (size_t i = 0; i < fanout; ++i) nodes += ((*memo)[k[i]] >> 16) & 0x7FF;
    assert(LEAF_SIZE < 3 || nodes <= 0x7FF); // Number of nodes in the subtree, stored in the middle 11 bits
    (*memo)[m] |= nodes << 16;
}

template <size_t LEAF_SIZE> static constexpr array<uint32_t, MAX_BUCKET_SIZE> fill_golomb_rice() {
    array<uint32_t, MAX_BUCKET_SIZE> memo{0};
    size_t s = 0;
    for (; s <= LEAF_SIZE; ++s) memo[s] = bij_memo[s] << 27 | (s > 1) << 16 | bij_memo[s];
    for (; s < MAX_BUCKET_SIZE; ++s) _fill_golomb_rice<LEAF_SIZE>(s, &memo);
    return memo;
}

// Computes the Golomb modulu of a splitting (for statistics purposes only)

template <size_t LEAF_SIZE> static constexpr uint64_t split_golomb_b(const int m) {
    array<int, MAX_FANOUT> k{0};

    size_t fanout = 0, unit = 0;
    SplittingStrategy<LEAF_SIZE>::split_params(m, fanout, unit);

    k[fanout - 1] = m;
    for (size_t i = 0; i < fanout - 1; ++i) {
        k[i] = unit;
        k[fanout - 1] -= k[i];
    }

    double sqrt_prod = 1;
    for (size_t i = 0; i < fanout; ++i) sqrt_prod *= ce::sqrt(k[i]);

    const double p = ce::sqrt(m) / (ce::pow(2 * M_PI, (fanout - 1.) / 2) * sqrt_prod);
    return ce::ceil(-log(2 - p) / log1p(-p)); // Golomb modulus
}

#define first_hash(k, len) spooky(k, len, 0)
#define golomb_param(m) (memo[m] >> 27)
#define skip_bits(m) (memo[m] & 0xFFFF)
#define skip_nodes(m) ((memo[m] >> 16) & 0x7FF)

/** David Stafford's (http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)
 * 13th variant of the 64-bit finalizer function in Austin Appleby's
 * MurmurHash3 (https://github.com/aappleby/smhasher).
 *
 * @param z a 64-bit integer.
 * @return a 64-bit integer obtained by mixing the bits of `z`.
 */
#ifdef __CUDACC__
__forceinline__ __host__ __device__
#endif
uint64_t remix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

/**
  * 32-bit finalizer function in Austin Appleby's MurmurHash3 (https://github.com/aappleby/smhasher).
  */
#ifdef __CUDACC__
__forceinline__ __host__ __device__
#endif
uint32_t remix32(uint32_t z) {
    z ^= z >> 16;
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae35;
    z ^= z >> 16;
    return z;
}

#ifdef __CUDACC__
__forceinline__ __host__ __device__
#endif
uint32_t remap(uint64_t x, uint32_t n) {
    constexpr bool USE_64_BITS = false;
    if constexpr (USE_64_BITS) {
        constexpr int masklen = 48;
        constexpr uint64_t mask = (uint64_t(1) << masklen) - 1;
        return ((remix(x) & mask) * n) >> masklen;
    } else {
        constexpr int masklen = 16;
        constexpr uint32_t mask = (uint32_t(1) << masklen) - 1;
        return ((remix32(uint32_t(x >> 32) ^ uint32_t(x)) & mask) * n) >> masklen;
    }
}

/**
 *
 * A class for storing minimal perfect hash functions. The template
 * parameter decides how large a leaf will be. Larger leaves imply
 * slower construction, but less space and faster evaluation.
 * This is the super class for GPURecSplit and SIMDRecSplit.
 * It is not the super class of RecSplit since it does some things
 * different, e.g., the queries.
 *
 * @tparam LEAF_SIZE the size of a leaf; typicals value range from 6 to 8
 * for fast, small maps, or up to 16 for very compact functions.
 * @tparam AT a type of memory allocation out of sux::util::AllocType.
 */

template <size_t LEAF_SIZE, util::AllocType AT = util::AllocType::MALLOC, bool USE_BIJECTIONS_ROTATE = false>
class AbstractParallelRecSplit {
  protected:
    using SplitStrat = SplittingStrategy<LEAF_SIZE>;

    static constexpr size_t _leaf = LEAF_SIZE;
    static constexpr size_t lower_aggr = SplitStrat::lower_aggr;
    static constexpr size_t upper_aggr = SplitStrat::upper_aggr;

    // For each bucket size, the Golomb-Rice parameter (upper 8 bits) and the number of bits to
    // skip in the fixed part of the tree (lower 24 bits).
    static constexpr array<uint32_t, MAX_BUCKET_SIZE> memo = get_golomb_rice_memo(LEAF_SIZE);

    static constexpr bool use_bijections_rotate = USE_BIJECTIONS_ROTATE;

    size_t bucket_size;
    size_t nbuckets;
    size_t keys_count;
    RiceBitVector<AT> descriptors;
    DoubleEF<AT, true> ef;

    // Maps a 128-bit to a bucket using the first 64-bit half.
    inline uint64_t hash128_to_bucket(const hash128_t &hash) const { return remap128(hash.first, nbuckets); }

    void parallelPartition(hash128_t *input, vector<uint64_t> &sorted,
                           vector<uint64_t> &bucket_size_acc, size_t num_threads) {
        assert(sorted.size() >= keys_count);
        assert(bucket_size_acc.size() == nbuckets + 1);

        sorter::sortParallel_hash128_t(input, keys_count, num_threads);

        // For most reasonable input sizes, doing this sequentially is faster
#ifndef PARALLEL_PARTITION_PREFIX
        size_t i = 0;
        const hash128_t *it = input;
        const hash128_t *end = input + keys_count;
        for (size_t bucket = 0; bucket < nbuckets; bucket++) {
            bucket_size_acc.at(bucket) = i;
            while (hash128_to_bucket(*it) == bucket && it != end) {
                sorted[i] = it->second;
                i++;
                it++;
            }
        }
        bucket_size_acc[nbuckets] = keys_count;
#else
        const size_t keysPerThread = keys_count / num_threads + 1;
        std::vector<std::thread> threads;
        for (size_t t = 0; t < num_threads; t++) {
            threads.emplace_back([&, t] {
                size_t begin = std::min(keysPerThread * t, keys_count);
                size_t end = std::min(keysPerThread * (t + 1), keys_count);

                // Align begin and end to next bucket borders
                while (hash128_to_bucket(input[end]) == hash128_to_bucket(input[end + 1]) && end < keys_count) {
                    end++;
                }
                if (t != 0) {
                    while (hash128_to_bucket(input[begin - 1]) == hash128_to_bucket(input[begin]) && begin < end) {
                        begin++;
                    }
                }
                size_t firstBucket = (t == 0) ? 0 : hash128_to_bucket(input[begin]);
                size_t lastBucket = (t == num_threads - 1) ? nbuckets : hash128_to_bucket(input[end]);

                size_t previousBucket = firstBucket - 1;
                for (size_t i = begin; i < end; i++) {
                    sorted[i] = input[i].second;
                    size_t bucket = hash128_to_bucket(input[i]);
                    while (previousBucket != bucket) {
                        previousBucket++;
                        bucket_size_acc[previousBucket] = i;
                    }
                }
                while (previousBucket < lastBucket) {
                    previousBucket++;
                    bucket_size_acc[previousBucket] = keys_count;
                }
            });
        }
        for (size_t t = 0; t < num_threads; t++) {
            threads.at(t).join();
        }
#endif
    }

  public:
    // TODO: why isn't this const?
    /** Returns the value associated with the given 128-bit hash.
     *
     * Note that this method is mainly useful for benchmarking.
     * @param hash a 128-bit hash.
     * @return the associated value.
     */
    size_t operator()(const hash128_t &hash) {
        const size_t bucket = hash128_to_bucket(hash);
        uint64_t cum_keys, cum_keys_next, bit_pos;
        ef.get(bucket, cum_keys, cum_keys_next, bit_pos);

        // Number of keys in this bucket
        size_t m = cum_keys_next - cum_keys;
        auto reader = descriptors.reader();
        reader.readReset(bit_pos, skip_bits(m));
        int level = 0;

        while (m > upper_aggr) { // fanout = 2
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remap(hash.second + d + start_seed[level], m);

            const uint32_t split = ((uint16_t(m / 2 + upper_aggr - 1) / upper_aggr)) * upper_aggr;
            if (hmod < split) {
                m = split;
            } else {
                reader.skipSubtree(skip_nodes(split), skip_bits(split));
                m -= split;
                cum_keys += split;
            }
            level++;
        }
        if (m > lower_aggr) {
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remap(hash.second + d + start_seed[NUM_START_SEEDS - 3], m);

            const int part = uint16_t(hmod) / lower_aggr;
            m = min(lower_aggr, m - part * lower_aggr);
            cum_keys += lower_aggr * part;
            if (part) reader.skipSubtree(skip_nodes(lower_aggr) * part, skip_bits(lower_aggr) * part);
        }

        if (m > _leaf) {
            const auto d = reader.readNext(golomb_param(m));
            const size_t hmod = remap(hash.second + d + start_seed[NUM_START_SEEDS - 2], m);

            const int part = uint16_t(hmod) / _leaf;
            m = min(_leaf, m - part * _leaf);
            cum_keys += _leaf * part;
            if (part) reader.skipSubtree(part, skip_bits(_leaf) * part);
        }

        const auto b = reader.readNext(golomb_param(m));
        const uint64_t x = hash.second + b + start_seed[NUM_START_SEEDS - 1];
        uint64_t leaf_val;
        if (use_bijections_rotate && m == _leaf) {
            const uint64_t rotation = b % _leaf;
            leaf_val = remap(x - rotation, _leaf);
            if (hash.second % 2 != 0) // is not left => need to rotate
                leaf_val = (leaf_val + rotation) % _leaf;
        } else {
            leaf_val = remap(x, m);
        }
        return cum_keys + leaf_val;
    }

    /** Returns the value associated with the given key.
     *
     * @param key a key.
     * @return the associated value.
     */
    size_t operator()(const string &key) { return operator()(first_hash(key.c_str(), key.size())); }

    /** Returns an estimate of the size in bits of this structure. */
    size_t getBits() {
        return ef.bitCountCumKeys() + ef.bitCountPosition() + descriptors.getBits() + 8 * sizeof(*this);
    }
};

} // namespace sux::function
