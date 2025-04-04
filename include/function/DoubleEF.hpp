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

#pragma once

#include <sux/util/Vector.hpp>
#include "../support/common.hpp"
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>
#include <cassert>

#ifdef SIMD
#include "../util/SimdUtils.hpp"
#endif

#ifndef LOG2Q
#define LOG2Q 8
#endif

namespace bez::function {

using namespace bez;

/** A double Elias-Fano list.
 *
 * This class exists solely to implement RecSplit.
 * @tparam AT a type of memory allocation out of util::AllocType.
 * @tparam FAST_CONSTRUCTION true if the construction should be accelerated by
 *           using SIMD and other techniques.
 */

template <sux::util::AllocType AT = sux::util::AllocType::MALLOC, bool FAST_CONSTRUCTION = false>
class DoubleEF {
  private:
    static constexpr uint64_t log2q = LOG2Q;
    static constexpr uint64_t q = 1 << log2q;
    static constexpr uint64_t q_mask = q - 1;
    static constexpr uint64_t super_q = 1 << 14;
    static constexpr uint64_t super_q_mask = super_q - 1;
    static constexpr uint64_t q_per_super_q = super_q / q;
    static constexpr uint64_t super_q_size = 1 + q_per_super_q / 4;
    sux::util::Vector<uint64_t, AT> lower_bits, upper_bits_position, upper_bits_cum_keys, jump;
    uint64_t lower_bits_mask_cum_keys, lower_bits_mask_position;

    uint64_t num_buckets, u_cum_keys, u_position;
    uint64_t l_position, l_cum_keys;
    int64_t cum_keys_min_delta, min_diff;
    uint64_t bits_per_key_fixed_point;

    __inline static void set(sux::util::Vector<uint64_t, AT> &bits, const uint64_t pos) { bits[pos / 64] |= 1ULL << pos % 64; }

    __inline static void set_bits(sux::util::Vector<uint64_t, AT> &bits, const uint64_t start, const int width, const uint64_t value) {
        const uint64_t mask = ((UINT64_C(1) << width) - 1) << start % 8;
        uint64_t t;
        memcpy(&t, (uint8_t *)&bits + start / 8, 8);
        t = (t & ~mask) | value << start % 8;
        memcpy((uint8_t *)&bits + start / 8, &t, 8);
    }

    __inline size_t lower_bits_size_words() const { return ((num_buckets + 1) * (l_cum_keys + l_position) + 63) / 64 + 1; }

    __inline size_t cum_keys_size_words() const { return (num_buckets + 1 + (u_cum_keys >> l_cum_keys) + 63) / 64; }

    __inline size_t position_size_words() const { return (num_buckets + 1 + (u_position >> l_position) + 63) / 64; }

    __inline size_t jump_size_words() const {
        size_t size = (num_buckets / super_q) * super_q_size * 2;                                        // Whole blocks
        if (num_buckets % super_q != 0) size += (1 + ((num_buckets % super_q + q - 1) / q + 3) / 4) * 2; // Partial block
        return size;
    }

    friend std::ostream &operator<<(std::ostream &os, const DoubleEF<AT, FAST_CONSTRUCTION> &ef) {
        os.write((char *)&ef.num_buckets, sizeof(ef.num_buckets));
        os.write((char *)&ef.u_cum_keys, sizeof(ef.u_cum_keys));
        os.write((char *)&ef.u_position, sizeof(ef.u_position));
        os.write((char *)&ef.cum_keys_min_delta, sizeof(ef.cum_keys_min_delta));
        os.write((char *)&ef.min_diff, sizeof(ef.min_diff));
        os.write((char *)&ef.bits_per_key_fixed_point, sizeof(ef.bits_per_key_fixed_point));

        os << ef.lower_bits;
        os << ef.upper_bits_cum_keys;
        os << ef.upper_bits_position;
        os << ef.jump;
        return os;
    }

    friend std::istream &operator>>(std::istream &is, DoubleEF<AT, FAST_CONSTRUCTION> &ef) {
        is.read((char *)&ef.num_buckets, sizeof(ef.num_buckets));
        is.read((char *)&ef.u_cum_keys, sizeof(ef.u_cum_keys));
        is.read((char *)&ef.u_position, sizeof(ef.u_position));
        is.read((char *)&ef.cum_keys_min_delta, sizeof(ef.cum_keys_min_delta));
        is.read((char *)&ef.min_diff, sizeof(ef.min_diff));
        is.read((char *)&ef.bits_per_key_fixed_point, sizeof(ef.bits_per_key_fixed_point));

        ef.l_position = ef.u_position / (ef.num_buckets + 1) == 0 ? 0 : lambda(ef.u_position / (ef.num_buckets + 1));
        ef.l_cum_keys = ef.u_cum_keys / (ef.num_buckets + 1) == 0 ? 0 : lambda(ef.u_cum_keys / (ef.num_buckets + 1));
        assert(ef.l_cum_keys * 2 + ef.l_position <= 56);

        ef.lower_bits_mask_cum_keys = (UINT64_C(1) << ef.l_cum_keys) - 1;
        ef.lower_bits_mask_position = (UINT64_C(1) << ef.l_position) - 1;

        is >> ef.lower_bits;
        is >> ef.upper_bits_cum_keys;
        is >> ef.upper_bits_position;
        is >> ef.jump;
        return is;
    }

  public:
    DoubleEF() {}

    DoubleEF(const std::vector<uint64_t> &cum_keys, const std::vector<uint64_t> &position) {
        assert(cum_keys.size() == position.size());
        num_buckets = cum_keys.size() - 1;

        bits_per_key_fixed_point = (uint64_t(1) << 20) * (position[num_buckets] / (double)cum_keys[num_buckets]);

        min_diff = std::numeric_limits<int64_t>::max() / 2;
        cum_keys_min_delta = std::numeric_limits<int64_t>::max() / 2;
        int64_t prev_bucket_bits = 0;

        uint64_t i = 1;
#ifdef SIMD
        if constexpr (FAST_CONSTRUCTION) {
            FullVecQ cum_keys_min_delta_vec = cum_keys_min_delta;
            FullVecQ min_diff_vec = min_diff;
            FullVecQ prev_bucket_bits_vec = 0;
            for (; i + FULL_VEC_64_COUNT - 1 <= num_buckets; i += FULL_VEC_64_COUNT) {
                FullVecQ cc, nkeys_delta, positions;
                cc.load(&cum_keys[i]);
                nkeys_delta.load(&cum_keys[i - 1]);
                nkeys_delta = cc - nkeys_delta;
                cum_keys_min_delta_vec = min(cum_keys_min_delta_vec, nkeys_delta);
                positions.load(&position[i]);
                FullVecQ bucket_bits = positions - (bits_per_key_fixed_point * cc >> 20);
                prev_bucket_bits_vec |= shift_bytes_up<8>(FullVecC(bucket_bits));
                min_diff_vec = min(min_diff_vec, bucket_bits - prev_bucket_bits_vec);
                prev_bucket_bits_vec = shift_bytes_down<FULL_VEC_64_COUNT * 8 - 8>(FullVecC(bucket_bits));
            }
            cum_keys_min_delta = horizontal_min(cum_keys_min_delta_vec);
            min_diff = horizontal_min(min_diff_vec);
            prev_bucket_bits = prev_bucket_bits_vec[0];
        }
#endif // SIMD
        for (; i <= num_buckets; ++i) {
            const int64_t nkeys_delta = cum_keys[i] - cum_keys[i - 1];
            cum_keys_min_delta = min(cum_keys_min_delta, nkeys_delta);
            const int64_t bucket_bits = int64_t(position[i]) - int64_t(bits_per_key_fixed_point * cum_keys[i] >> 20);
            min_diff = min(min_diff, bucket_bits - prev_bucket_bits);
            prev_bucket_bits = bucket_bits;
        }

        u_position = int64_t(position[num_buckets]) - int64_t(bits_per_key_fixed_point * cum_keys[num_buckets] >> 20) - int64_t(num_buckets * min_diff) + 1;
        l_position = u_position / (num_buckets + 1) == 0 ? 0 : lambda(u_position / (num_buckets + 1));
        u_cum_keys = cum_keys[num_buckets] - num_buckets * cum_keys_min_delta + 1;
        l_cum_keys = u_cum_keys / (num_buckets + 1) == 0 ? 0 : lambda(u_cum_keys / (num_buckets + 1));
        assert(l_cum_keys * 2 + l_position <= 56); // To be able to perform a single unaligned read

#ifdef MORESTATS
        printf("Elias-Fano l (cumulative): %lu\n", l_cum_keys);
        printf("Elias-Fano l (positions): %lu\n", l_position);
        printf("Elias-Fano u (cumulative): %lu\n", u_cum_keys);
        printf("Elias-Fano u (positions): %lu\n", u_position);
#endif

        lower_bits_mask_cum_keys = (UINT64_C(1) << l_cum_keys) - 1;
        lower_bits_mask_position = (UINT64_C(1) << l_position) - 1;

        const uint64_t words_lower_bits = lower_bits_size_words();
        lower_bits.size(words_lower_bits);
        const uint64_t words_cum_keys = cum_keys_size_words();
        upper_bits_cum_keys.size(words_cum_keys);
        const uint64_t words_position = position_size_words();
        upper_bits_position.size(words_position);

        i = 0;
        uint64_t cum_delta = 0, bit_delta = 0;
#ifdef SIMD
        if constexpr (FAST_CONSTRUCTION) {
#ifdef SIMDRS_512_BIT
            FullVecUq iVec(0, 1, 2, 3, 4, 5, 6, 7);
#else
            FullVecUq iVec(0, 1, 2, 3);
#endif
            FullVecUq cum_delta_vec = cum_keys_min_delta * iVec;
            FullVecUq bit_delta_vec = min_diff * iVec;
            for (; i + FULL_VEC_64_COUNT - 1 <= num_buckets; i += FULL_VEC_64_COUNT, iVec += FULL_VEC_64_COUNT,
                cum_delta_vec += FULL_VEC_64_COUNT * cum_keys_min_delta, bit_delta_vec += FULL_VEC_64_COUNT * min_diff) {
                FullVecUq cc;
                cc.load(&cum_keys[i]);
                FullVecQ positions;
                positions.load(&position[i]);
                const FullVecUq pval_delta = FullVecUq(positions - FullVecQ(bits_per_key_fixed_point * cc >> 20)) - bit_delta_vec;
                const FullVecUq cc_delta = cc - cum_delta_vec;
                if (l_cum_keys + l_position != 0) {
                    const FullVecUq start = iVec * (l_cum_keys + l_position);
                    const FullVecUq startMod8 = start & 7;
                    FullVecUq lower_bits_vec = 0;
                    if (l_cum_keys != 0) {
                        const FullVecUq value = cc_delta & lower_bits_mask_cum_keys;
                        lower_bits_vec = shift(value, startMod8);
                    }
                    if (l_position != 0) {
                        const FullVecUq value = pval_delta & lower_bits_mask_position;
                        lower_bits_vec |= shift(value, startMod8 + l_cum_keys);
                    }
                    const FullVecUq pointer = (uint64_t)&lower_bits + (start >> 3);
                    uint64_t values[FULL_VEC_64_COUNT];
                    lower_bits_vec.store(values);
                    uint8_t *pointers[FULL_VEC_64_COUNT];
                    pointer.store(pointers);
                    for (uint32_t j = 0; j < FULL_VEC_64_COUNT; ++j) {
                        uint64_t t;
                        memcpy(&t, pointers[j], 8);
                        t |= values[j]; // no mask necessary
                        memcpy(pointers[j], &t, 8);
                    }
                }
                {
                    FullVecUq pos = (cc_delta >> (uint32_t)l_cum_keys) + iVec;
                    FullVecUq pointer = (uint64_t)&upper_bits_cum_keys + ((pos >> 6) << 3);
                    uint64_t *pointers[FULL_VEC_64_COUNT];
                    pointer.store(pointers);
                    FullVecUq value = shift(FullVecUq(1), pos & 63);
                    uint64_t values[FULL_VEC_64_COUNT];
                    value.store(values);
                    for (uint32_t j = 0; j < FULL_VEC_64_COUNT; ++j)
                        *(pointers[j]) |= values[j];
                }
                {
                    FullVecUq pos = (pval_delta >> (uint32_t)l_position) + iVec;
                    FullVecUq pointer = (uint64_t)&upper_bits_position + ((pos >> 6) << 3);
                    uint64_t *pointers[FULL_VEC_64_COUNT];
                    pointer.store(pointers);
                    FullVecUq value = shift(FullVecUq(1), pos & 63);
                    uint64_t values[FULL_VEC_64_COUNT];
                    value.store(values);
                    for (uint32_t j = 0; j < FULL_VEC_64_COUNT; ++j)
                        *(pointers[j]) |= values[j];
                }
            }
            cum_delta = cum_keys_min_delta * i;
            bit_delta = min_diff * i;
        }
#endif // SIMD
        for (; i <= num_buckets; i++, cum_delta += cum_keys_min_delta, bit_delta += min_diff) {
            if (l_cum_keys != 0) set_bits(lower_bits, i * (l_cum_keys + l_position), l_cum_keys, (cum_keys[i] - cum_delta) & lower_bits_mask_cum_keys);
            set(upper_bits_cum_keys, ((cum_keys[i] - cum_delta) >> l_cum_keys) + i);

            const auto pval = int64_t(position[i]) - int64_t(bits_per_key_fixed_point * cum_keys[i] >> 20);
            if (l_position != 0) set_bits(lower_bits, i * (l_cum_keys + l_position) + l_cum_keys, l_position, (pval - bit_delta) & lower_bits_mask_position);
            set(upper_bits_position, ((pval - bit_delta) >> l_position) + i);
        }

        const uint64_t jump_words = jump_size_words();
        jump.size(jump_words);

        if constexpr (FAST_CONSTRUCTION) {
            static_assert(q <= super_q && super_q % q == 0 && q > 64);
            for (uint64_t i = 0, c = 0, q_counter = q_mask, last_super_q = 0; i < words_cum_keys; i++) {
                const uint64_t bits = upper_bits_cum_keys[i];
                const uint64_t popcount = nu(bits);
                q_counter += popcount;
                c += popcount;
                if (q_counter >= q) {
                    q_counter -= q;
                    const uint64_t bit_num = i * 64 + select64(bits, popcount - q_counter - 1);
                    const uint64_t idx = (c / super_q) * (super_q_size * 2);
                    uint64_t exact_c = c - q_counter - 1;
                    if ((exact_c & super_q_mask) == 0) jump[idx] = last_super_q = bit_num;
                    const uint64_t offset = bit_num - last_super_q;
                    assert(offset < (1 << 16));
                    ((uint16_t *)(&jump + idx + 2))[2 * ((exact_c & super_q_mask) / q)] = offset;
                }
            }

            for (uint64_t i = 0, c = 0, q_counter = q_mask, last_super_q = 0; i < words_position; i++) {
                const uint64_t bits = upper_bits_position[i];
                const uint64_t popcount = nu(bits);
                q_counter += popcount;
                c += popcount;
                if (q_counter >= q) {
                    q_counter -= q;
                    const uint64_t bit_num = i * 64 + select64(bits, popcount - q_counter - 1);
                    const uint64_t idx = (c / super_q) * (super_q_size * 2);
                    uint64_t exact_c = c - q_counter - 1;
                    if ((exact_c & super_q_mask) == 0) jump[idx + 1] = last_super_q = bit_num;
                    const uint64_t offset = bit_num - last_super_q;
                    assert(offset < (1 << 16));
                    ((uint16_t *)(&jump + idx + 2))[2 * ((exact_c & super_q_mask) / q) + 1] = offset;
                }
            }
        } else {
            for (uint64_t i = 0, c = 0, last_super_q = 0; i < words_cum_keys; i++) {
                for (int b = 0; b < 64; b++) {
                    if (upper_bits_cum_keys[i] & UINT64_C(1) << b) {
                        if ((c & super_q_mask) == 0) jump[(c / super_q) * (super_q_size * 2)] = last_super_q = i * 64 + b;
                        if ((c & q_mask) == 0) {
                            const uint64_t offset = i * 64 + b - last_super_q;
                            if (offset >= (1 << 16)) abort();
                            ((uint16_t *)(&jump + (c / super_q) * (super_q_size * 2) + 2))[2 * ((c % super_q) / q)] = offset;
                        }
                        c++;
                    }
                }
            }

            for (uint64_t i = 0, c = 0, last_super_q = 0; i < words_position; i++) {
                for (int b = 0; b < 64; b++) {
                    if (upper_bits_position[i] & UINT64_C(1) << b) {
                        if ((c & super_q_mask) == 0) jump[(c / super_q) * (super_q_size * 2) + 1] = last_super_q = i * 64 + b;
                        if ((c & q_mask) == 0) {
                            const uint64_t offset = i * 64 + b - last_super_q;
                            if (offset >= (1 << 16)) abort();
                            ((uint16_t *)(&jump + (c / super_q) * (super_q_size * 2) + 2))[2 * ((c % super_q) / q) + 1] = offset;
                        }
                        c++;
                    }
                }
            }
        }

#ifndef NDEBUG
        for (uint64_t i = 0; i < num_buckets; i++) {
            uint64_t x, x2, y;

            get(i, x, x2, y);
            assert(x == cum_keys[i]);
            assert(x2 == cum_keys[i + 1]);
            assert(y == position[i]);

            get(i, x, y);
            assert(x == cum_keys[i]);
            assert(y == position[i]);
        }
#endif
    }

    void get(const uint64_t i, uint64_t &cum_keys, uint64_t &cum_keys_next, uint64_t &position) {
        const uint64_t pos_lower = i * (l_cum_keys + l_position);
        uint64_t lower;
        memcpy(&lower, (uint8_t *)&lower_bits + pos_lower / 8, 8);
        lower >>= pos_lower % 8;

        const uint64_t jump_super_q = (i / super_q) * super_q_size * 2;
        const uint64_t jump_inside_super_q = (i % super_q) / q;
        const uint64_t jump_cum_keys = jump[jump_super_q] + ((uint16_t *)(&jump + jump_super_q + 2))[2 * jump_inside_super_q];
        const uint64_t jump_position = jump[jump_super_q + 1] + ((uint16_t *)(&jump + jump_super_q + 2))[2 * jump_inside_super_q + 1];

        uint64_t curr_word_cum_keys = jump_cum_keys / 64;
        uint64_t curr_word_position = jump_position / 64;
        uint64_t window_cum_keys = upper_bits_cum_keys[curr_word_cum_keys] & UINT64_C(-1) << jump_cum_keys % 64;
        uint64_t window_position = upper_bits_position[curr_word_position] & UINT64_C(-1) << jump_position % 64;
        uint64_t delta_cum_keys = i & q_mask;
        uint64_t delta_position = i & q_mask;

        for (uint64_t bit_count; (bit_count = nu(window_cum_keys)) <= delta_cum_keys; delta_cum_keys -= bit_count) window_cum_keys = upper_bits_cum_keys[++curr_word_cum_keys];
        for (uint64_t bit_count; (bit_count = nu(window_position)) <= delta_position; delta_position -= bit_count) window_position = upper_bits_position[++curr_word_position];

        const uint64_t select_cum_keys = select64(window_cum_keys, delta_cum_keys);
        const int64_t cum_delta = i * cum_keys_min_delta;
        cum_keys = ((curr_word_cum_keys * 64 + select_cum_keys - i) << l_cum_keys | (lower & lower_bits_mask_cum_keys)) + cum_delta;

        lower >>= l_cum_keys;
        const int64_t bit_delta = i * min_diff;
        position = ((curr_word_position * 64 + select64(window_position, delta_position) - i) << l_position | (lower & lower_bits_mask_position)) + bit_delta +
                   int64_t(bits_per_key_fixed_point * cum_keys >> 20);

        window_cum_keys &= (-1ULL << select_cum_keys) << 1;
        while (window_cum_keys == 0) window_cum_keys = upper_bits_cum_keys[++curr_word_cum_keys];

        lower >>= l_position;
        cum_keys_next = ((curr_word_cum_keys * 64 + rho(window_cum_keys) - i - 1) << l_cum_keys | (lower & lower_bits_mask_cum_keys)) + cum_delta + cum_keys_min_delta;
    }

    void get(const uint64_t i, uint64_t &cum_keys, uint64_t &position) {
        const uint64_t pos_lower = i * (l_cum_keys + l_position);
        uint64_t lower;
        memcpy(&lower, (uint8_t *)&lower_bits + pos_lower / 8, 8);
        lower >>= pos_lower % 8;

        const uint64_t jump_super_q = (i / super_q) * super_q_size * 2;
        const uint64_t jump_inside_super_q = (i % super_q) / q;
        const uint64_t jump_cum_keys = jump[jump_super_q] + ((uint16_t *)(&jump + jump_super_q + 2))[2 * jump_inside_super_q];
        const uint64_t jump_position = jump[jump_super_q + 1] + ((uint16_t *)(&jump + jump_super_q + 2))[2 * jump_inside_super_q + 1];

        uint64_t curr_word_cum_keys = jump_cum_keys / 64;
        uint64_t curr_word_position = jump_position / 64;
        uint64_t window_cum_keys = upper_bits_cum_keys[curr_word_cum_keys] & UINT64_C(-1) << jump_cum_keys % 64;
        uint64_t window_position = upper_bits_position[curr_word_position] & UINT64_C(-1) << jump_position % 64;
        uint64_t delta_cum_keys = i & q_mask;
        uint64_t delta_position = i & q_mask;

        for (uint64_t bit_count; (bit_count = nu(window_cum_keys)) <= delta_cum_keys; delta_cum_keys -= bit_count) window_cum_keys = upper_bits_cum_keys[++curr_word_cum_keys];
        for (uint64_t bit_count; (bit_count = nu(window_position)) <= delta_position; delta_position -= bit_count) window_position = upper_bits_position[++curr_word_position];

        const uint64_t select_cum_keys = select64(window_cum_keys, delta_cum_keys);
        const size_t cum_delta = i * cum_keys_min_delta;
        cum_keys = ((curr_word_cum_keys * 64 + select_cum_keys - i) << l_cum_keys | (lower & lower_bits_mask_cum_keys)) + cum_delta;

        lower >>= l_cum_keys;
        const int64_t bit_delta = i * min_diff;
        position = ((curr_word_position * 64 + select64(window_position, delta_position) - i) << l_position | (lower & lower_bits_mask_position)) + bit_delta +
                   int64_t(bits_per_key_fixed_point * cum_keys >> 20);
    }

    uint64_t bitCountCumKeys() { return (num_buckets + 1) * l_cum_keys + num_buckets + 1 + (u_cum_keys >> l_cum_keys) + jump_size_words() / 2; }

    uint64_t bitCountPosition() { return (num_buckets + 1) * l_position + num_buckets + 1 + (u_position >> l_position) + jump_size_words() / 2; }
};

} // namespace sux::function
