# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True

from cython cimport floating


def _trim_end_gaps(
    const unsigned char[:, :] bits,
    Py_ssize_t[::1] starts,
    Py_ssize_t[::1] stops,
):
    r"""Identify the terminal gap-free region of a multiple alignment.

    Parameters
    ----------
    bits : ndarray of uint8 of shape (n_sequences, n_positions)
        Bit array representing gaps in the alignment.
    starts : ndarray of int of shape (n_sequences,)
        Start position of terminal gap-free region of each sequence.
        (i.e., index of the first position within non-gap)
    stops : ndarray of int of shape (n_sequences,)
        Stop position of terminal gap-free region of each sequence.
        (i.e., index of the first position after non-gap)

    Notes
    -----
    If a sequence only contains gaps, double zeros will be assigned.

    This function works for both the original alignment (columns are positions) and
    the alignment path (columns are segments).

    Parallelization is possible, although perhaps unnecessary.

    """
    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t n = bits.shape[1]

    for i in range(bits.shape[0]):
        k = n
        for j in range(n):
            if bits[i, j] == 0:
                starts[i] = k = j
                break
        if k == n:  # gap-only sequence
            starts[i] = stops[i] = 0
            continue
        for j in range(n - 1, k - 1, -1):
            if bits[i, j] == 0:
                stops[i] = j + 1
                break


def _multi_align_score(
    const Py_ssize_t[:, :] seqs,
    const unsigned char[:, :] bits,
    const Py_ssize_t[:] lens,
    const Py_ssize_t[::1] starts,
    const Py_ssize_t[::1] stops,
    const floating[:, :] submat,
    floating match,
    floating mismatch,
    floating gap_open,
    floating gap_extend,
    bint free_ends,
):
    """Calculate sum-of-pairs (SP) alignment score of aligned sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Character array represented by ASCII code or alphabet index.
    bits : ndarray of uint8 of shape (n_sequences, n_segments)
        Bit array representing gap status in the alignment path.
    lens : ndarray of int of shape (n_sequences, n_segments)
        Lengths of segments in the alignment path.
    starts : ndarray of int of shape (n_sequences,)
        Start position of terminal gap-free region of each sequence.
    stops : ndarray of int of shape (n_sequences,)
        Stop position of terminal gap-free region of each sequence.
    submat : ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.
    match : float
        Match score.
    mismatch : float
        Mismatch score.
    gap_open : float
        Gap opening penalty.
    gap_extend : float
        Gap extension penalty.
    free_ends : bool
        Whether terminal gaps are free from penalization.

    Returns
    -------
    float
        Alignment score.

    """
    # TODO: Some array parameters can be [::1].

    # This function employs an algorithm that is more complex than intuition. Instead
    # of iterating over all positions and accumulatively adding score or cost of each
    # position, it operates on the alignment path, which divides the alignment into
    # segments representing altering status. This design permits the calculation of
    # gap costs on the entire contiguous gap rather than by each gap position. It not
    # only saves compute, but also enables complex gap penalty schemes, such as convex
    # and dual affine penalties, although they are not currently implemented.

    # This algorithm calculates the alignment score between each pair of sequences and
    # sums the results. Therefore, it has a time complexity of O(n^2), which isn't
    # ideal especially when there are many sequences. Alternatively, one can iterate
    # over positions, calculate character frequencies, then calculate the overall
    # score accordingly. However, this design only works for linear gap penalty.

    # TODO: Implement a separate algorithm for linear gap penalty on many sequences.

    cdef Py_ssize_t i1, i2, j, k
    cdef Py_ssize_t start, end, pos
    cdef int L, cumL, prev, curr

    cdef bint is_submat = submat.shape[0] > 0

    cdef floating score = 0  # cumulative alignment score
    cdef Py_ssize_t n = seqs.shape[0]  # number of sequences

    # calculate alignment score of each pair of sequences and sum up
    # TODO: The current algorithm can be parallelized.
    for i1 in range(n):
        for i2 in range(i1 + 1, n):

            # determine start and end segment indices to iterate over
            if free_ends:
                start = max(starts[i1], starts[i2])
                end = min(stops[i1], stops[i2])
            else:
                start = min(starts[i1], starts[i2])
                end = max(stops[i1], stops[i2])

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
                    score -= gap_open + cumL * gap_extend

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
                score -= gap_open + cumL * gap_extend

    return score
