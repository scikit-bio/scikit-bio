# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np
cimport numpy as cnp
from skbio.alignment import TabularMSA, AlignPath
from skbio.sequence import DNA
from numpy cimport int32_t, uint8_t, int8_t, int64_t, float32_t
cnp.import_array()

# Scoring constants
cdef uint8_t GAP = 1
cdef uint8_t GAP_SCORE = 255
cdef int32_t NEG_INF = -2147483647

# Alignment score wrapper function for export
def align_score(seq1, seq2, subMatrix, gap_open, gap_extend):
    # Has sequences through substitution matrix
    seq1_idx = subMatrix._char_hash[DNA(seq1)._bytes]
    seq2_idx = subMatrix._char_hash[DNA(seq2)._bytes]
    # Run though helper function to get ouput
    return _align_score(seq1_idx, seq2_idx, subMatrix._data.astype(np.float32), gap_open, gap_extend)

# Alignment score main function
cdef float32_t _align_score(const uint8_t[::1] seq1, const uint8_t[::1] seq2, const float32_t[:, :] subMatrix, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
    # Initialize  variables
    cdef uint8_t i
    cdef uint8_t state = 0
    cdef float32_t score = 0

    # Iterate through sequences
    for i in range(seq1.shape[0]):
        # Gap in seq1
        if seq1[i] == GAP_SCORE:
            # Either no current gap or gap in seq2
            score += gap_extend
            if state == 0 or state == 2:
                score += gap_open
                state = 1
        # Gap in seq2
        elif seq2[i] == GAP_SCORE:
            # Either no current gap or gap in seq1
            score += gap_extend
            if state == 0 or state == 1:
                score += gap_open
                state = 2
        # No gaps
        else:
            score += subMatrix[seq1[i], seq2[i]] 
            state = 0
    # Return final result
    return score