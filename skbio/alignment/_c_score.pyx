# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION


def _trim_terminal_gaps(
    unsigned char[:, :] bits,
    Py_ssize_t[::1] starts,
    Py_ssize_t[::1] ends,
):
    r"""Identify the terminal gap-free region of a multiple alignment.

    Parameters
    ----------
    bits : ndarray of uint8 of shape (n_sequences, n_positions)
        Bit array representing gaps in the alignment.
    starts : ndarray of int of shape (n_sequences,)
        Start position of terminal gap-free region of each sequence.
        (i.e., index of the first position within non-gap)
    ends : ndarray of int of shape (n_sequences,)
        End position of terminal gap-free region of each sequence.
        (i.e., index of the first position after non-gap)

    Notes
    -----
    If a sequence only contains gaps, double zeros will be assigned.

    This function works for both the original alignment (columns are positions) and
    the alignment path (columns are segments).

    """
    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t n = bits.shape[1]

    # TODO: Parallelization is possible.
    for i in range(bits.shape[0]):
        k = n
        for j in range(n):
            if bits[i, j] == 0:
                starts[i] = k = j
                break
        if k == n:  # gap-only sequence
            starts[i] = ends[i] = 0
            continue
        for j in range(n - 1, k - 1, -1):
            if bits[i, j] == 0:
                ends[i] = j + 1
                break


def _multi_align_score(
    unsigned char[:, :] seqs,
    unsigned char[:, :] bits,
    Py_ssize_t[:] lens,
    Py_ssize_t[:] starts,
    Py_ssize_t[:] ends,
    double[:, :] submat,
    double match,
    double mismatch,
    double gap_open,
    double gap_extend,
    bint terminal_gaps,
):
    """Calculate sum-of-pairs (SP) alignment score of aligned sequences.
    """
    cdef Py_ssize_t i1, i2, j, k
    cdef Py_ssize_t start, end, pos
    cdef int L, cumL, prev, curr

    cdef bint is_submat = submat.shape[0] > 0

    cdef double score = 0  # cumulative alignment score
    cdef Py_ssize_t n = seqs.shape[0]  # number of sequences

    # calculate alignment score of each pair of sequences and sum up
    # TODO: This can be parallelized.
    for i1 in range(n):
        for i2 in range(i1 + 1, n):

            # determine start and end segment indices to iterate over
            if terminal_gaps:
                start = min(starts[i1], starts[i2])
                end = max(ends[i1], ends[i2])
            else:
                start = max(starts[i1], starts[i2])
                end = min(ends[i1], ends[i2])

            # determine start position in the alignment
            pos = 0
            for j in range(start):
                pos += lens[j]

            prev = 0  # previous state
            cumL = 0  # cumulative gap length

            # iterate over segments
            for j in range(start, end):
                L = lens[j]
                curr = bits[i1, j] + bits[i2, j] * 2

                # gap in both sequences: ignore
                if curr == 3:
                    pos += L
                    continue

                # end of previous gap
                if prev and curr != prev:
                    score -= gap_open + (cumL - 1) * gap_extend

                # non-gap in both sequences
                if curr == 0:

                    # iterate by position within segment
                    if is_submat:
                        for k in range(pos, pos + L):
                            score += submat[seqs[i1, k], seqs[i2, k]]
                    else:
                        for k in range(pos, pos + L):
                            if seqs[i1, k] == seqs[i2, k]:
                                score += match
                            else:
                                score += mismatch

                # gap in either sequence
                else:

                    # gap in the same sequence continues
                    if curr == prev:
                        cumL += L

                    # gap switches to a different sequence
                    else:
                        cumL = L

                if curr != prev:
                    prev = curr

                pos += L

            # handle last gap
            if prev:
                score -= gap_open + (cumL - 1) * gap_extend

    return score
