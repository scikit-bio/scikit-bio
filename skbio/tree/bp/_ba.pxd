# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Derived from improved-octo-waddle (https://github.com/biocore/improved-octo-waddle)
# originally authored by Daniel McDonald, distributed under the Modified BSD License.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

cdef extern from "<inttypes.h>":
    ctypedef unsigned int uint64_t

cdef extern from "bit_array.h":
    struct BIT_ARRAY:
        pass
    ctypedef uint64_t bit_index_t

    # allocations
    BIT_ARRAY* bit_array_create(bit_index_t nbits)
    void bit_array_free(BIT_ARRAY* bitarray)
    BIT_ARRAY* bit_array_clone(const BIT_ARRAY* bitarr)

    # utility
    char* bit_array_to_str(const BIT_ARRAY* bitarr, char* str)
    bit_index_t bit_array_length(const BIT_ARRAY* bit_arr)
    bit_index_t bit_array_num_bits_set(const BIT_ARRAY* bitarr)

    # bit juggling
    void bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b) nogil
    void bit_array_toggle_bit(BIT_ARRAY* bitarr, bit_index_t b) nogil
    char bit_array_get_bit(const BIT_ARRAY* bitarr, bit_index_t b) nogil
    void bit_array_clear_bit(BIT_ARRAY* bitarr, bit_index_t b)

    # logical operations
    void bit_array_and(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
    void bit_array_or(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
    void bit_array_xor(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
    void bit_array_not(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2)

    # cyclic shifting
    void bit_array_cycle_right(BIT_ARRAY* bitarr, bit_index_t dist)
    void bit_array_cycle_left (BIT_ARRAY* bitarr, bit_index_t dist)
