#pragma once

#include <sux/support/SpookyV2.hpp>

namespace bez::function {

/** 128-bit hashes.
 *
 * In the construction of GPURecSplit, keys are replaced with instances
 * of this class using SpookyHash, first thing.
 * Moreover, it is possible to build and query GPURecSplit instances using 128-bit
 * random hashes only (mainly for benchmarking purposes).
 */
typedef struct __hash128_t {
    uint64_t first, second;

    bool operator<(const __hash128_t &o) const { return first < o.first || second < o.second; }

    __hash128_t(const uint64_t first, const uint64_t second) {
        this->first = first;
        this->second = second;
    }
} hash128_t;

/** Convenience function hashing a key a returning a __hash128_t
 *
 * @param data a pointer to the key.
 * @param length the length in bytes of the key.
 * @param seed an additional seed.
 */

hash128_t inline spooky(const void *data, const size_t length, const uint64_t seed) {
    uint64_t h0 = seed, h1 = seed;
    SpookyHash::Hash128(data, length, &h0, &h1);
    return {h1, h0};
}

}
