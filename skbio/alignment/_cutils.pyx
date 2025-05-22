# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True

from cython cimport floating
from libc.math cimport INFINITY, fabs, fabsf
from libc.float cimport FLT_MAX, DBL_MAX


cdef inline floating xabs(floating x) noexcept nogil:
    """Calculate absolute value of a floating-point number.

    It takes a fused floating type variable (float or double) and dispatches the
    corresponding C function (fabsf or fabs). Since it is inlined, there is little
    overhead compared with directly calling the C function.

    """
    if floating is float:
        return fabsf(x)
    else:
        return fabs(x)


cdef inline floating lbound(bint local, floating gap_extend) noexcept nogil:
    """Get lower bound of alignment score.

    `gap_extend` is supplied to let the compiler know the floating type (float or
    double). Otherwise it will complain: "Return type is not specified as argument
    type".

    Currently it is not used, because `-INFINITY` is good enough.

    """
    if local:
        return 0
    elif floating is float:
        return -FLT_MAX
    else:
        return -DBL_MAX


def _fill_linear_matrix(
    floating[:, ::1] scomat,
    const floating[:, ::1] query,
    const Py_ssize_t[::1] target,
    floating gap_extend,
    bint local,
):
    """Calculate optimal scores over the alignment matrix with linear gap penalty.

    Parameters:
    ----------
    scomat : memoryview of ndarray of shape (m, n)
        Main matrix.
    query : memoryview of ndarray of float of shape (m, n_symbols)
        Query profile.
    target : memoryview of ndarray of int of shape (n,)
        Target sequence.
    gap_extend : floating
        Gap extension penalty.
    local : bint
        Local (True) or global (False) alignment.

    Notes
    -----
    Matrix filling is the most computationally intensive step in pairwise sequence
    alignment using dynamic programming. It is O(mn). The current algorithm adopts
    the classic and simple double loops method.

    Efficiency is achieved through memory layout and compiler's auto-optimization.
    All arrays are C-contiguous, facilitating row-wise iteration, which is what the
    inner loop (j-indexed) does. The loop content is simple as possible, involving
    a single `max` and no branches. The linear and affine versions are separated,
    whereas the global and local versions are unified with a lower bound constant
    defined prior to the loops. Directions are not saved during matrix filling, but
    instead re-calculated during traceback, which is O(m + n) anyway.

    `query` is a pre-computed profile of the query sequence. Each row represents a
    character in the query sequence and each column represents a symbol in the
    subsitution matrix. `target` is a vector of indices in the subsitution matrix,
    and each elements represents a character in the target sequence. Therefore,
    `query[i - 1, target[j - 1]]` directly looks up the substitution score between
    i-th character in query and j-th character in target.

    Although this algorithm does not utilize various manual optimization techniques
    that have been developed and shown effective in pairwise sequence alignment, such
    as wavefront sweep, loop tiling, SIMD striping and prefix scan, they are noted
    here for reference.

    """
    cdef floating bound = 0 if local else -INFINITY
    cdef Py_ssize_t m1 = scomat.shape[0], n1 = scomat.shape[1]
    cdef Py_ssize_t i, j
    cdef floating* row

    for i in range(1, m1):
        row = &query[i - 1, 0]
        for j in range(1, n1):
            scomat[i, j] = max(
                scomat[i - 1, j - 1] + row[target[j - 1]],
                scomat[i, j - 1] - gap_extend,
                scomat[i - 1, j] - gap_extend,
                bound,
            )


def _fill_affine_matrices(
    floating[:, ::1] scomat,
    floating[:, ::1] insmat,
    floating[:, ::1] delmat,
    const floating[:, ::1] query,
    const Py_ssize_t[::1] target,
    floating gap_open,
    floating gap_extend,
    bint local,
):
    """Calculate optimal scores over the alignment matrix with affine gap penalty.

    Parameters:
    ----------
    scomat : memoryview of ndarray of shape (m, n)
        Main matrix.
    insmat : memoryview of ndarray of shape (m, n)
        Insertion matrix.
    delmat : memoryview of ndarray of shape (m, n)
        Deletion matrix.
    query : memoryview of ndarray of float of shape (m, n_symbols)
        Query profile.
    target : memoryview of ndarray of int of shape (n,)
        Target sequence.
    gap_open : floating
        Gap opening penalty.
    gap_extend : floating
        Gap extension penalty.
    local : bint
        Local (True) or global (False) alignment.

    See Also
    --------
    _fill_linear_matrix

    """
    cdef floating bound = 0 if local else -INFINITY
    cdef floating gap_open_extend = gap_open + gap_extend
    cdef Py_ssize_t m1 = scomat.shape[0], n1 = scomat.shape[1]
    cdef Py_ssize_t i, j
    cdef floating sub_, ins_, del_
    cdef floating* row

    for i in range(1, m1):
        row = &query[i - 1, 0]
        for j in range(1, n1):

            # substitution (diagonal)
            sub_ = scomat[i - 1, j - 1] + row[target[j - 1]]

            # open a new insertion or extend a previous insertion (horizontal)
            ins_ = insmat[i, j] = max(
                scomat[i, j - 1] - gap_open_extend,
                insmat[i, j - 1] - gap_extend,
            )

            # open a new deletion or extend a previous deletion (vertical)
            del_ = delmat[i, j] = max(
                scomat[i - 1, j] - gap_open_extend,
                delmat[i - 1, j] - gap_extend,
            )

            scomat[i, j] = max(sub_, ins_, del_, bound)


def _trace_one_linear(
    unsigned char[::1] path,
    Py_ssize_t pos,
    Py_ssize_t i,
    Py_ssize_t j,
    floating[:, ::1] scomat,
    floating gap_extend,
    bint local,
    floating eps,
):
    """Traceback across matrix body with linear gap penalty.

    Parameters:
    ----------
    path : memoryview of ndarray of shape (n_positions,)
        Dense alignment path.
    pos : int
        Current start position of the path.
    i : int
        Current row index in the matrix.
    j : int
        Current column index in the matrix.
    scomat : memoryview of ndarray of shape (m, n)
        Main matrix.
    gap_extend : floating
        Gap extension penalty.
    local : bint
        Local (True) or global (False) alignment.
    eps : floating
        Absolute tolerance.

    Returns
    -------
    int
        Updated start position of the path.
    int
        Updated row index in the matrix.
    int
        Updated column index in the matrix.

    Notes
    -----
    Only one optimal alignment path is traced. The priority is:

        deletion > insertion > substitution

    This function does not involve the original sequences or the substitution matrix,
    thereby improving memory efficiency.

    Altschul & Erickson (1986) [1]_ pointed out that the original Gotoh algorithm
    could fail in some situations. Specifically, if the traceback process does not
    jump between matrices, it won't be able to distinguish ties between creating a
    new gap vs. extending an existing gap. The current function follows the modified
    Gotoh algorithm.

    References
    ----------
    .. [1] Altschul, S. F., & Erickson, B. W. (1986). Optimal sequence alignment using
       affine gap costs. Bull Math Biol, 48, 603-616.

    """
    cdef floating score, gap_score

    # will stop when reaching either edge of the matrix
    while i and j:
        score = scomat[i, j]
        if local and xabs(score) <= eps:
            break
        gap_score = gap_extend + score
        pos -= 1

        # deletion (vertical; gap in seq2)
        if xabs(scomat[i - 1, j] - gap_score) <= eps:
            path[pos] = 2
            i -= 1
        # insertion (horizontal; gap in seq1)
        elif xabs(scomat[i, j - 1] - gap_score) <= eps:
            path[pos] = 1
            j -= 1
        # substitution (diagonal; no gap)
        else:
            path[pos] = 0
            i -= 1
            j -= 1

    return pos, i, j


def _trace_one_affine(
    unsigned char[::1] path,
    Py_ssize_t pos,
    Py_ssize_t i,
    Py_ssize_t j,
    floating[:, ::1] scomat,
    floating[:, ::1] insmat,
    floating[:, ::1] delmat,
    floating gap_extend,
    bint local,
    floating eps,
):
    """Traceback across matrix body with affine gap penalty.

    Parameters:
    ----------
    path : memoryview of ndarray of shape (n_positions,)
        Dense alignment path.
    pos : int
        Current start position of the path.
    i : int
        Current row index in the matrix.
    j : int
        Current column index in the matrix.
    scomat : memoryview of ndarray of shape (m, n)
        Main matrix.
    insmat : memoryview of ndarray of shape (m, n)
        Insertion matrix.
    delmat : memoryview of ndarray of shape (m, n)
        Deletion matrix.
    gap_extend : floating
        Gap extension penalty.
    local : bint
        Local (True) or global (False) alignment.
    eps : floating
        Absolute tolerance.

    Returns
    -------
    int
        Updated start position of the path.
    int
        Updated row index in the matrix.
    int
        Updated column index in the matrix.

    Notes
    -----
    Only one optimal alignment path is traced. The priority is:
    
        Main matrix: jumping to deletion matrix > jumping to insertion matrix >
          substitution.
        Deletion matrix: staying in deletion matrix > jumping to main matrix.
        Insertion matrix: staying in insertion matrix > jumping to main matrix.

    """
    cdef floating score
    cdef int mat = 0

    while i and j:
        score = scomat[i, j]
        if local and xabs(score) <= eps:
            break

        # deletion matrix (vertical; gap in seq2)
        if mat == 2:
            # extend an existing gap (stay in the current matrix),
            # or open a new gap (jump back to main matrix)
            if xabs(delmat[i - 1, j] - gap_extend - delmat[i, j]) > eps:
                mat = 0
            i -= 1
            pos -= 1
            path[pos] = 2

        # insertion matrix (horizontal; gap in seq1)
        elif mat == 1:
            # same as above
            if xabs(insmat[i, j - 1] - gap_extend - insmat[i, j]) > eps:
                mat = 0
            j -= 1
            pos -= 1
            path[pos] = 1

        # main matrix
        else:
            if xabs(delmat[i, j] - score) <= eps:  # jump to deletion matrix
                mat = 2
            elif xabs(insmat[i, j] - score) <= eps:  # jump to insertion matrix
                mat = 1
            else:  # substitution (diagonal; no gap)
                i -= 1
                j -= 1
                pos -= 1
                path[pos] = 0

    return pos, i, j


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
    const Py_ssize_t[:, ::1] seqs,
    const unsigned char[:, :] bits,
    const Py_ssize_t[:] lens,
    const Py_ssize_t[::1] starts,
    const Py_ssize_t[::1] stops,
    const floating[:, ::1] submat,
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
    lens : ndarray of int of shape (n_segments,)
        Lengths of segments in the alignment path.
    starts : ndarray of int of shape (n_sequences,)
        Start position of terminal gap-free region of each sequence.
    stops : ndarray of int of shape (n_sequences,)
        Stop position of terminal gap-free region of each sequence.
    submat : ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.
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

                # non-gap in both sequences: iterate by position within segment
                if curr == 0:
                    for k in range(pos, pos + L):
                        score += submat[seqs[i1, k], seqs[i2, k]]

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


# def local_align_fill(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, const cnp.float64_t[:, :] subMatrix, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
#     """
#     Perform traceback for global sequence alignment using precomputed score matrices.

#     Parameters:
#     ----------
#     D_view : float32_t[:, :]
#         Empty match penalty matrix.
#     P_view : float32_t[:, :]
#         Empty insert penalty matrix.
#     Q_view : float32_t[:, :]
#         Empty deletion penalty matrix.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     Nothing
#     """
    
#     # Initialize variables
#     cdef Py_ssize_t i, j, idx
#     cdef Py_ssize_t m = seq1.shape[0]
#     cdef Py_ssize_t n = seq2.shape[0]

#     # Initialize 0, 0 for each matrix
#     D_view[0, 0] = 0
#     P_view[0, 0] = NEG_INF
#     Q_view[0, 0] = NEG_INF

#     # Initialize first row for each matrix
#     for i in range(1, m + 1):
#         D_view[i, 0] = 0
#         Q_view[i, 0] = NEG_INF

#     # Initialize first column for each matrix
#     for j in range(1, n + 1):
#         D_view[0, j] = 0
#         P_view[0, j] = NEG_INF

#     # Main fill loop
#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             # Fill insertion scores
#             P_view[i, j] = max(
#                 D_view[i - 1, j] + gap_open + gap_extend,
#                 P_view[i - 1, j] + gap_extend
#             )
#             # Fill deletion scores
#             Q_view[i, j] = max(
#                 D_view[i, j - 1] + gap_open + gap_extend,
#                 Q_view[i, j - 1] + gap_extend
#             )
#             # Fill match/mismatch scores, can't go below 0
#             D_view[i, j] = max(
#                 D_view[i - 1, j - 1] + subMatrix[seq1[i-1], seq2[j-1]],
#                 P_view[i, j],
#                 Q_view[i, j],
#                 0
#             )


# import numpy as np
# cimport numpy as cnp
# from skbio.alignment import TabularMSA, AlignPath
# from skbio.sequence import DNA
# from numpy cimport int32_t, uint8_t, int8_t, float32_t, float64_t
# cnp.import_array()

# # Scoring constants
# cdef uint8_t GAP = 1
# cdef uint8_t GAP_SCORE = 255
# cdef int32_t NEG_INF = -2147483647

# # Structure definitions
# cdef struct Index:
#     Py_ssize_t i 
#     Py_ssize_t j

# cdef struct TracebackRes:
#     float32_t score
#     Py_ssize_t idx
#     Index start
#     Index end

# def pairwise_align(seq1, seq2, subMatrix, gap_open, gap_extend, scope):
#     """
#     Public facing align function, runs Cythonized algorithm.

#     Parameters:
#     ----------
#     seq1 : string or any scikit-bio DNA compatible class
#         The first sequence, represented as a string, DNA, or any other DNA-compatible class.
#     seq2 : string or any scikit-bio DNA compatible class
#         The second sequence, represented as a string, DNA, or any other DNA-compatible class.
#     subMatrix : skbio.sequence.SubstitutionMatrix
#         Substitution matrix, taken from the scikit-bio sequence module.
#     gap_open : float
#         The penalty for opening a gap.
#     gap_extend : float
#         The penalty for extending an existing gap.
#     scope : string
#         Either "local" or "global", indicates alignment type

#     Returns:
#     -------
#     TabularMSA
#         TabularMSA that represents final alignment.
#     """

#     # Hash seq1 and seq2 through subMatrix
#     asDNA = [DNA(seq1), DNA(seq2)]
#     seq1_idx = subMatrix._char_hash[asDNA[0]._bytes]
#     seq2_idx = subMatrix._char_hash[asDNA[1]._bytes]
    
#     # Calculate max alignment length
#     max_len = seq1_idx.shape[0] + seq2_idx.shape[0]

#     # Initialize final alignment buffers
#     aligned_seq1 = np.zeros(max_len, dtype=np.uint8)
#     aligned_seq2 = np.zeros(max_len, dtype=np.uint8)

#     # Create tuple of arguments and pass into algorithm
#     args = (seq1_idx, seq2_idx, aligned_seq1, aligned_seq2, subMatrix._data.astype(float), gap_open, gap_extend)
#     if(scope.lower() == "local"):
#         args += (local_align_fill, local_align_trace)
#     else:
#         args += (global_align_fill, global_align_trace)
#     res = arg_pass(*args)

#     # Format result into a usable format and conver tto alignment path
#     result = (aligned_seq1[res['idx']+1:aligned_seq1.size], aligned_seq2[res['idx']+1:aligned_seq2.size])
#     if result[0].size == 0 and result[1].size == 1:
#         raise Exception("No local alignment found!")
#     path = AlignPath.from_bits(np.vstack(result))

#     # Cut out local alignment path from originals
#     if(scope.lower() == "local"):
#         asDNA = [DNA(seq1[res['start']['i']:res['end']['i']]), DNA(seq2[res['start']['j']:res['end']['j']])]
    
#     # Convert to TabularMSA and return
#     return TabularMSA.from_path_seqs(path, asDNA), res['score']

# def arg_pass(const uint8_t[::1] seq1, const uint8_t[::1] seq2, cnp.ndarray aligned_seq1, cnp.ndarray aligned_seq2, const float64_t[:, ::1] subMatrix, const int8_t gap_open, const int8_t gap_extend, fill_func, trace_func):
#     '''
#     Helper function to connect Python code with Cython code
#     '''

#     return _pairwise_align(seq1, seq2, aligned_seq1, aligned_seq2, subMatrix, gap_open, gap_extend, fill_func, trace_func)

# cdef TracebackRes _pairwise_align(const uint8_t[::1] seq1, const uint8_t[::1] seq2, cnp.ndarray aligned_seq1, cnp.ndarray aligned_seq2, const float64_t[:, ::1] subMatrix, const int8_t gap_open, const int8_t gap_extend, fill_func, trace_func) noexcept:
#     """
#     Calculate alignment between seq1 and seq2 using given scope and penalties.

#     Parameters:
#     ----------
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     aligned_seq1 : cnp.ndarray
#         Array for storing the aligned version of `seq1`.
#     aligned_seq2 : cnp.ndarray
#         Array for storing the aligned version of `seq2`.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.
#     fill_func : function
#         Function to fill P, D, and Q matrices.
#     trace_func : function
#         Function to trace result through filled matrices.

#     Returns:
#     -------
#     TracebackRes
#         A structure containing the score as well as the traceback end index.
#     """
    
#     # Initialize variables
#     cdef Py_ssize_t m = seq1.shape[0]
#     cdef Py_ssize_t n = seq2.shape[0]
#     cdef Py_ssize_t max_len = m + n
#     cdef Py_ssize_t idx
#     cdef TracebackRes res

#     # Initialize matrices
#     cdef D = np.empty((m + 1, n + 1), dtype=np.float32)
#     cdef P = np.empty((m + 1, n + 1), dtype=np.float32)
#     cdef Q = np.empty((m + 1, n + 1), dtype=np.float32)

#     # Get memoryview representations for Cython
#     cdef float32_t[:, ::1] D_view = D
#     cdef float32_t[:, ::1] P_view = P
#     cdef float32_t[:, ::1] Q_view = Q
#     cdef uint8_t[::1] aligned_seq1_view = aligned_seq1
#     cdef uint8_t[::1] aligned_seq2_view = aligned_seq2

#     # Run and return fill and trace functions on variables
#     fill_func(D_view, P_view, Q_view, subMatrix, seq1, seq2, gap_open, gap_extend)
#     res = trace_func(D_view, P_view, Q_view, max_len, aligned_seq1_view, aligned_seq2_view, seq1, seq2, subMatrix, gap_open, gap_extend)
#     return res

# cdef void global_align_fill(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, const cnp.float64_t[:, :] subMatrix, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
#     """
#     Perform traceback for global sequence alignment using precomputed score matrices.

#     Parameters:
#     ----------
#     D_view : float32_t[:, :]
#         Empty match penalty matrix.
#     P_view : float32_t[:, :]
#         Empty insert penalty matrix.
#     Q_view : float32_t[:, :]
#         Empty deletion penalty matrix.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     Nothing
#     """

#     # Initialize variables
#     cdef Py_ssize_t i, j
#     cdef Py_ssize_t m = seq1.shape[0]
#     cdef Py_ssize_t n = seq2.shape[0]

#     # Initialize 0, 0 for each matrix
#     D_view[0, 0] = 0
#     P_view[0, 0] = NEG_INF
#     Q_view[0, 0] = NEG_INF

#     # Initialize first row for each matrix
#     for i in range(1, m + 1):
#         D_view[i, 0] = gap_open + gap_extend * i
#         Q_view[i, 0] = NEG_INF

#     # Initialize first column for each matrix
#     for j in range(1, n + 1):
#         D_view[0, j] = gap_open + gap_extend * j
#         P_view[0, j] = NEG_INF
    
#     # Main fill loop
#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             # Fill insertion scores
#             P_view[i, j] = max(
#                 D_view[i - 1, j] + gap_open + gap_extend,
#                 P_view[i - 1, j] + gap_extend
#             )
#             # Fill deletion scores
#             Q_view[i, j] = max(
#                 D_view[i, j - 1] + gap_open + gap_extend,
#                 Q_view[i, j - 1] + gap_extend
#             )
#             # Fill match/mismatch scores
#             D_view[i, j] = max(
#                 D_view[i - 1, j - 1] + subMatrix[seq1[i-1], seq2[j-1]],
#                 P_view[i, j],
#                 Q_view[i, j]
#             )

# cdef void local_align_fill(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, const cnp.float64_t[:, :] subMatrix, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
#     """
#     Perform traceback for global sequence alignment using precomputed score matrices.

#     Parameters:
#     ----------
#     D_view : float32_t[:, :]
#         Empty match penalty matrix.
#     P_view : float32_t[:, :]
#         Empty insert penalty matrix.
#     Q_view : float32_t[:, :]
#         Empty deletion penalty matrix.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     Nothing
#     """
    
#     # Initialize variables
#     cdef Py_ssize_t i, j, idx
#     cdef Py_ssize_t m = seq1.shape[0]
#     cdef Py_ssize_t n = seq2.shape[0]

#     # Initialize 0, 0 for each matrix
#     D_view[0, 0] = 0
#     P_view[0, 0] = NEG_INF
#     Q_view[0, 0] = NEG_INF

#     # Initialize first row for each matrix
#     for i in range(1, m + 1):
#         D_view[i, 0] = 0
#         Q_view[i, 0] = NEG_INF

#     # Initialize first column for each matrix
#     for j in range(1, n + 1):
#         D_view[0, j] = 0
#         P_view[0, j] = NEG_INF

#     # Main fill loop
#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             # Fill insertion scores
#             P_view[i, j] = max(
#                 D_view[i - 1, j] + gap_open + gap_extend,
#                 P_view[i - 1, j] + gap_extend
#             )
#             # Fill deletion scores
#             Q_view[i, j] = max(
#                 D_view[i, j - 1] + gap_open + gap_extend,
#                 Q_view[i, j - 1] + gap_extend
#             )
#             # Fill match/mismatch scores, can't go below 0
#             D_view[i, j] = max(
#                 D_view[i - 1, j - 1] + subMatrix[seq1[i-1], seq2[j-1]],
#                 P_view[i, j],
#                 Q_view[i, j],
#                 0
#             )

# cdef TracebackRes global_align_trace(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, Py_ssize_t max_len, uint8_t[::1] aligned_seq1_view, uint8_t[::1] aligned_seq2_view, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const cnp.float64_t[:, :] subMatrix, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
#     """
#     Perform traceback for global sequence alignment using precomputed score matrices.

#     Parameters:
#     ----------
#     D_view : float32_t[:, :]
#         Filled match penalty matrix.
#     P_view : float32_t[:, :]
#         Filled insertion penalty matrix.
#     Q_view : float32_t[:, :]
#         Filled deletion penalty matrix.
#     max_len : Py_ssize_t
#         The maximum allowed length for the traceback.
#     aligned_seq1_view : uint8_t[::1]
#         Preallocated buffer for storing the aligned version of `seq1`.
#     aligned_seq2_view : uint8_t[::1]
#         Preallocated buffer for storing the aligned version of `seq2`.
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     TracebackRes
#         A structure containing the score as well as the traceback end index.
#     """

#     # Initialize variables
#     cdef float32_t score
#     cdef uint8_t current_matrix
#     cdef TracebackRes res
#     cdef Py_ssize_t idx = max_len - 1
#     cdef Py_ssize_t m = seq1.shape[0]
#     cdef Py_ssize_t n = seq2.shape[0]

#     # Check for matrix with highest score and set as starting point
#     if D_view[m, n] >= P_view[m, n] and D_view[m, n] >= Q_view[m, n]:
#         score = D_view[m, n]
#         current_matrix = 0
#     elif P_view[m, n] >= Q_view[m, n]:
#         score = P_view[m, n]
#         current_matrix = 1
#     else:
#         score = Q_view[m, n]
#         current_matrix = 2

#     # Main traceback while loop
#     while m > 0 and n > 0:
#         if current_matrix == 0:
#             # Go diagonal if match/mismatch
#             if D_view[m, n] == D_view[m-1, n-1] + subMatrix[seq1[m-1], seq2[n-1]]:
#                 m -= 1
#                 n -= 1
#             # Swap to insert or delete matrix otherwise
#             elif D_view[m, n] == P_view[m, n]:
#                 current_matrix = 1
#                 continue
#             else:
#                 current_matrix = 2
#                 continue
#         elif current_matrix == 1:
#             # If gap is opened, go back to match/mismatch
#             aligned_seq2_view[idx] = GAP
#             if P_view[m, n] == D_view[m-1, n] + gap_open + gap_extend:
#                 current_matrix = 0
#             m -= 1
#         else:
#             # If gap is opened, go back to match/mismatch
#             aligned_seq1_view[idx] = GAP
#             if Q_view[m, n] == D_view[m, n-1] + gap_open + gap_extend:
#                 current_matrix = 0 
#             n -= 1
#         # Move back in aligned sequences
#         idx -= 1
    
#     # Set any remaining spaces as gaps
#     while m > 0:
#         aligned_seq2_view[idx] = GAP
#         m -= 1
#         idx -= 1
#     while n > 0:
#         aligned_seq1_view[idx] = GAP
#         n -= 1
#         idx -= 1
    
#     # Store results in return structure
#     res.score = score
#     res.idx = idx
#     return res

# cdef TracebackRes local_align_trace(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, Py_ssize_t max_len, uint8_t[::1] aligned_seq1_view, uint8_t[::1] aligned_seq2_view, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const cnp.float64_t[:, :] subMatrix, const int8_t gap_open, const int8_t gap_extend) noexcept:
#     """
#     Perform traceback for local sequence alignment using precomputed score matrices.

#     Parameters:
#     ----------
#     D_view : float32_t[:, :]
#         Filled match penalty matrix.
#     P_view : float32_t[:, :]
#         Filled insertion penalty matrix.
#     Q_view : float32_t[:, :]
#         Filled deletion penalty matrix.
#     loc : Index
#         Index object, unused in local trace but required for parity.
#     max_len : Py_ssize_t
#         The maximum allowed length for the traceback.
#     aligned_seq1_view : uint8_t[::1]
#         Preallocated buffer for storing the aligned version of `seq1`.
#     aligned_seq2_view : uint8_t[::1]
#         Preallocated buffer for storing the aligned version of `seq2`.
#     seq1 : const cnp.uint8_t[::1]
#         The first sequence hashed through subMatrix.
#     seq2 : const cnp.uint8_t[::1]
#         The second sequence hashed through subMatrix.
#     subMatrix : const cnp.float64_t[:, :]
#         The substitution matrix used for scoring alignments.
#     gap_open : int8_t
#         The penalty for opening a gap.
#     gap_extend : int8_t
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     TracebackRes
#         A structure containing the score as well as the start and end indices of the traceback.
#     """

#     # Initialize variables
#     cdef Index P_loc, D_loc, Q_loc, end
#     cdef float32_t score
#     cdef uint8_t current_matrix
#     cdef TracebackRes res
#     cdef Py_ssize_t idx = max_len - 1

#     # Get highest score from each matrix
#     P_loc.i, P_loc.j = np.unravel_index(np.argmax(P_view), np.shape(P_view))
#     D_loc.i, D_loc.j = np.unravel_index(np.argmax(D_view), np.shape(D_view))
#     Q_loc.i, Q_loc.j = np.unravel_index(np.argmax(Q_view), np.shape(Q_view))

#     # Check for matrix with highest overall score and set as starting point
#     if D_view[D_loc.i, D_loc.j] >= P_view[P_loc.i, P_loc.j] and D_view[D_loc.i, D_loc.j] >= Q_view[Q_loc.i, Q_loc.j]:
#         current_matrix = 0
#         score = D_view[D_loc.i, D_loc.j]
#         loc = D_loc
#     elif P_view[P_loc.i, P_loc.j] >= Q_view[Q_loc.i, Q_loc.j]:
#         current_matrix = 1
#         score = P_view[P_loc.i, P_loc.j]
#         loc = P_loc
#     else:
#         current_matrix = 2
#         score = Q_view[Q_loc.i, Q_loc.j]
#         loc = Q_loc
    
#     # Save endpoint for return
#     end = loc

#     # Main traceback while loop
#     while loc.i > 0 and loc.j > 0:
#         if current_matrix == 0:
#             # Go diagonal if match/mismatch
#             if D_view[loc.i, loc.j] == D_view[loc.i-1, loc.j-1] + subMatrix[seq1[loc.i-1], seq2[loc.j-1]]:
#                 loc.i -= 1
#                 loc.j -= 1
#             elif D_view[loc.i, loc.j] == P_view[loc.i, loc.j]:
#                 current_matrix = 1
#                 continue
#             else:
#                 current_matrix = 2
#                 continue
#             # Break once the end of the traceback is reached
#             if D_view[loc.i, loc.j] == 0:
#                 idx -= 1
#                 break
#         elif current_matrix == 1:
#             # If gap is opened, go back to match/mismatch
#             aligned_seq2_view[idx] = GAP
#             if P_view[loc.i, loc.j] == D_view[loc.i-1, loc.j] + gap_open + gap_extend:
#                 current_matrix = 0
#             loc.i -= 1
#         else:
#             # If gap is opened, go back to match/mismatch
#             aligned_seq1_view[idx] = GAP
#             if Q_view[loc.i, loc.j] == D_view[loc.i, loc.j-1] + gap_open + gap_extend:
#                 current_matrix = 0 
#             loc.j -= 1
#         # Move back in aligned sequences
#         idx -= 1

#     # Store results in return structure
#     res.score = score
#     res.idx = idx
#     res.start = loc
#     res.end = end
#     return res

# # Alignment score wrapper function for export
# def align_score(seq1, seq2, subMatrix, gap_open, gap_extend):
#     """
#     Computes the alignment score between two sequences using the given substitution matrix 
#     and gap penalties.

#     Parameters:
#     ----------
#     seq1 : string or any scikit-bio DNA compatible class
#         The first sequence, represented as a string, DNA, or any other DNA-compatible class.
#     seq2 : string or any scikit-bio DNA compatible class
#         The second sequence, represented as a string, DNA, or any other DNA-compatible class.
#     subMatrix : skbio.sequence.SubstitutionMatrix
#         Substitution matrix, taken from the scikit-bio sequence module.
#     gap_open : float
#         The penalty for opening a gap.
#     gap_extend : float
#         The penalty for extending an existing gap.

#     Returns:
#     -------
#     float
#         The computed alignment score.
#     """
#     # Has sequences through substitution matrix
#     seq1_idx = subMatrix._char_hash[DNA(seq1)._bytes]
#     seq2_idx = subMatrix._char_hash[DNA(seq2)._bytes]
#     # Run though helper function to get ouput
#     return _align_score(seq1_idx, seq2_idx, subMatrix._data.astype(np.float32), gap_open, gap_extend)

# # Alignment score main function
# cdef float32_t _align_score(const uint8_t[::1] seq1, const uint8_t[::1] seq2, const float32_t[:, :] subMatrix, const float32_t gap_open, const float32_t gap_extend) noexcept nogil:
#     # Initialize  variables
#     cdef Py_ssize_t i
#     cdef uint8_t state = 0
#     cdef float32_t score = 0

#     # Iterate through sequences
#     for i in range(seq1.shape[0]):
#         # Gap in seq1
#         if seq1[i] == GAP_SCORE:
#             # Either no current gap or gap in seq2
#             score += gap_extend
#             if state == 0 or state == 2:
#                 score += gap_open
#                 state = 1
#         # Gap in seq2
#         elif seq2[i] == GAP_SCORE:
#             # Either no current gap or gap in seq1
#             score += gap_extend
#             if state == 0 or state == 1:
#                 score += gap_open
#                 state = 2
#         # No gaps
#         else:
#             score += subMatrix[seq1[i], seq2[i]] 
#             state = 0
#     # Return final result
#     return score
