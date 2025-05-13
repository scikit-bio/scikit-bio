# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numbers import Real
from collections import namedtuple

import numpy as np

from skbio._base import SkbioObject
from skbio.util._decorator import classonlymethod
from skbio.sequence import Sequence, SubstitutionMatrix
from skbio.alignment import PairAlignPath
from ._c_pair import _fill_linear_matrix, _fill_affine_matrices


def pair_align(
    seq1,
    seq2,
    local=False,
    sub_score=(1, -1),
    gap_cost=2,
    free_ends=True,
    max_paths=1,
    tolerance=1e-5,
    keep_matrices=False,
):
    r"""Perform pairwise alignment of two sequences.

    .. versionadded:: 0.6.4

    Parameters
    ----------
    seq1 : GrammaredSequence or str
        The first sequence to be aligned.

    seq2 : GrammaredSequence or str
        The second sequence to be aligned.

    local : bool, optional
        Perform global alignment (False, default) or local alignment (True).

    sub_score : tuple of (float, float), SubstitutionMatrix, or str, optional
        Score of a substitution. May be two numbers (match, mismatch), a substitution
        matrix, or its name. See :func:`align_score` for instructions. Default is
        (1, -1).

    gap_cost : float or tuple of (float, float), optional
        Penalty of a gap. May be one (linear) or two numbers (affine). See
        :func:`align_score` for instructions and rationales. Default is -2.

    free_ends : bool, optional
        Whether gaps at the sequence terminals are free from penalization. See
        :func:`align_score` for instructions. Default is True.

    max_paths : int, optional
        Maximum number of alignment paths to return. Default is 1, which is generated
        through a performance-oriented traceback algorithm. A value larger than 1 will
        trigger a less efficient traceback algorithm to enumerate up to this number of
        paths. Setting it as None will return all paths. However, be cautious that the
        total number of paths may be extremely large and could stress the system.
        Setting it as 0 will disable traceback and return no path.

    tolerance : float, optional
        Absolute tolerance in comparing scores of alternative alignment paths. This is
        to ensure floating-point arithmetic safety when ``sub_score`` or ``gap_cost``
        involve decimal numbers. Default is 1e-5. Setting it to 0 or None will slightly
        increase performance, and is usually safe when there are only integers (e.g.,
        2.0), half integers (e.g., 2.5) or numbers with exact binary representation
        involved. Note: relative tolerance is not involved in the calculation.

    keep_matrices : bool, optional
        Whether to include the alignment matrix(ces) in the returned value. They are
        typically for diagnostic or educational purposes. Default is False, which lets
        the memory space free up after the function finishes.

    Returns
    -------
    score : float
        Optimal alignment score.

    paths : list of :class:`~skbio.alignment.PairAlignPath`
        Alignment paths. Up to ``max_paths`` paths will be returned. Note that all
        paths are optimal and share the same alignment score.

    matrices : list of ndarray of float of shape (m + 1, n + 1), optional
        Alignment matrices generated during the computation.

    See Also
    --------
    align_score
    skbio.alignment.PairAlignPath

    Notes
    -----
    This function implements the classical dynamic programming (DP) method for
    pairwise sequence alignment. This method is commonly known as the Needleman-Wunsch
    algorithm [1]_ for global alignment or the Smith-Waterman algorithm [2]_ for local
    alignment. These two algorithms are for linear gap penalty. When affine gap
    penalty is specified, the underlying method is the Gotoh algorithm [3]_, with
    later modifications [4]_.

    References
    ----------
    .. [1] Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the
       search for similarities in the amino acid sequence of two proteins. J Mol Biol,
       48(3), 443-453.

    .. [2] Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular
       subsequences. J Mol Biol, 147(1), 195-197.

    .. [3] Gotoh, O. (1982). An improved algorithm for matching biological sequences.
       J Mol Biol, 162(3), 705-708.

    .. [4] Altschul, S. F., & Erickson, B. W. (1986). Optimal sequence alignment using
       affine gap costs. Bull Math Biol, 48, 603-616.

    """
    # encode sequences
    seq1 = _seq_to_bytes(seq1)
    seq2 = _seq_to_bytes(seq2)
    m = seq1.size
    n = seq2.size

    # substitution matrix
    if isinstance(sub_score, str):
        sub_score = SubstitutionMatrix.by_name(sub_score)
    if isinstance(sub_score, SubstitutionMatrix):
        submat = _submat_from_sm(seq1, seq2, sub_score)
    # match/mismatch scores
    else:
        submat = _submat_from_mm(seq1, seq2, *sub_score)

    # affine or linear gap penalty
    if isinstance(gap_cost, Real):
        gap_open, gap_extend = 0, gap_cost
    else:
        gap_open, gap_extend = gap_cost

    # allocate alignment matrices
    matrices = _alloc_matrices(m, n, gap_open)

    # initialize alignment matrices
    _init_matrices(*matrices, gap_open, gap_extend, local, free_ends)

    # fill alignment matrices (quadratic; compute-intensive)
    if gap_open:
        _fill_affine_matrices(*matrices, submat, gap_open, gap_extend, local)
    else:
        _fill_linear_matrix(matrices[0], submat, gap_extend, local)

    # get optimal alignment score and corresponding stop(s)
    if max_paths == 1:
        score, stops = _one_stop(matrices[0], local, free_ends)
    else:
        score, stops = _all_stops(matrices[0], local, free_ends, tolerance)

    # no path is found
    if local and abs(score) <= tolerance:
        paths, stops = [], np.empty((0, 2), dtype=int)

    # traceback from each stop to reconstruct optimal alignment path(s)
    elif max_paths == 0:
        paths = []
    elif max_paths == 1:
        paths = [_traceback_one(*stops[0], matrices, gap_open, gap_extend, local)]
    else:
        paths = _traceback_all(
            stops, matrices, submat, gap_open, gap_extend, local, max_paths
        )

    # discard or keep matrices
    if not keep_matrices:
        matrices = []
    elif gap_open:
        _fill_nan(matrices)
    else:
        matrices = [matrices[0]]

    return PairAlignResult(float(score), paths, matrices)


# pairwise alignment result
PairAlignResult = namedtuple("PairAlignResult", ["score", "paths", "matrices"])


def _seq_to_bytes(seq):
    """Convert a sequence into bytes."""
    if isinstance(seq, Sequence):
        return seq._bytes
    elif isinstance(seq, str):
        return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    else:
        raise ValueError("Sequence must be a string or a `Sequence` object.")


def _moves(i, j, alnmat, seq1, seq2, mthmis, gap):
    """Calculate scores of moving in three directions."""
    return (
        # substitution (diagonal)
        alnmat[i - 1, j - 1] + mthmis[int(seq1[i - 1] != seq2[j - 1])],
        # TODO: alternative methods:
        # (a == b) * match + (a != b) * mismatch
        # match if a == b else mismatch
        # insertion (left to right)
        alnmat[i, j - 1] + gap,
        # deletion (upper to lower)
        alnmat[i - 1, j] + gap,
    )


### Create substitution matrix


def _submat_from_mm(seq1, seq2, match, mismatch):
    """Pre-compute a match/mismatch array to facilitate lookup."""
    match, mismatch = float(match), float(mismatch)
    return np.where(seq1[:, None] == seq2[None, :], match, mismatch)


def _submat_from_sm(seq1, seq2, sub_score):
    """Pre-compute a substitution array to facilitate lookup."""
    seq1 = sub_score._char_hash[seq1]
    seq2 = sub_score._char_hash[seq2]
    try:
        submat = sub_score._data[seq1][:, seq2]
    except IndexError:
        raise ValueError(
            "Sequences contain characters that are not present in the provided "
            "substitution matrix."
        )
    return np.ascontiguousarray(submat)


def _alloc_matrices(m, n, affine, dtype=np.float64):
    """Allocate alignment matrix(ces).

    Parameters
    ----------
    m : int
        Length of sequence 1.
    n : int
        Length of sequence 2.
    affine : bool
        Affine (True) or linear (False) gap penalty.
    dtype : type, optional
        Data type (np.float32 or np.float64)

    Returns
    -------
    scomat : ndarray of shape (m, n)
        Main matrix.
    insmat : ndarray of shape (m, n)
        Insertion matrix.
    delmat : ndarray of shape (m, n)
        Deletion matrix.

    """
    # Note: The array should be C-contiguous to facilitate row-wise iteration.
    # NumPy's default array is already C-contiguous. This is also enforced by the
    # unit test.
    shape = (m + 1, n + 1)
    scomat = np.empty(shape, dtype=dtype)
    if affine:
        insmat = np.empty(shape, dtype=dtype)
        delmat = np.empty(shape, dtype=dtype)
    else:
        insmat = None
        delmat = None
    return scomat, insmat, delmat


def _init_matrices(scomat, insmat, delmat, gap_open, gap_extend, local, free_ends):
    """Initialize alignment matrix(ces) by populating first column and row.

    Parameters
    ----------
    scomat : ndarray of shape (m, n)
        Main matrix.
    insmat : ndarray of shape (m, n)
        Insertion matrix.
    delmat : ndarray of shape (m, n)
        Deletion matrix.
    gap_open : float
        Gap opening penalty.
    gap_extend : float
        Gap extension penalty.
    local : bool
        Local or global alignment.
    free_ends : bool
        If end gaps are free from penalty.

    """
    m1, n1 = scomat.shape

    # initialize main scoring matrix
    scomat[0, 0] = 0
    if local or free_ends:
        scomat[1:m1, 0] = 0
        scomat[0, 1:n1] = 0
    else:
        series = np.arange(1, max(m1, n1), dtype=scomat.dtype)
        series *= -gap_extend
        if gap_open:
            series -= gap_open
        scomat[1:m1, 0] = series[: m1 - 1]
        scomat[0, 1:n1] = series[: n1 - 1]

    # initialize insertion and deletion matrices
    if gap_open:
        insmat[1:m1, 0] = -np.inf
        delmat[0, 1:n1] = -np.inf
        # series -= gap_open
        # insmat[1:m1, 0] = series[:m1 - 1]
        # delmat[0, 1:n1] = series[:n1 - 1]


def _fill_nan(matrices):
    """ "Fill empty cells of affine matrices with NaN before returning."""
    _, insmat, delmat = matrices
    insmat[0, :] = np.nan
    delmat[:, 0] = np.nan


def _one_stop(scomat, local, free_ends):
    """Locate one stop with optimal alignment score.

    Parameters
    ----------
    scomat : ndarray of shape (m, n)
        Main matrix.
    local : bool
        Local or global alignment.
    free_ends : bool
        If end gaps are free from penalty.

    Returns
    -------
    float
        Optimal alignment score.
    ndarray of int of shape (1, 2)
        Coordinates of one alignment stop.

    Notes
    -----
    When there is a tie, the smallest index (col, row) is chosen.

    """
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1

    # local alignment: maximum cell in the matrix
    if local:
        i, j = np.divmod(scomat.argmax(), n + 1)

    # semi-global alignment: maximum cell in the last column and row
    elif free_ends:
        i = scomat[:m, n].argmax()  # last column (ends with deletion)
        j = scomat[m, :].argmax()  # last row (ends with insertion)
        if scomat[i, n] >= scomat[m, j]:
            j = n
        else:
            i = m

    # global alignment: bottom right cell
    else:
        i, j = m, n

    return scomat[i, j], np.array([[i, j]])


def _all_stops(scomat, local, free_ends, eps=1e-7):
    """Locate all stops with optimal alignment score.

    Parameters
    ----------
    scomat : ndarray of shape (m, n)
        Main matrix.
    local : bool
        Local or global alignment.
    free_ends : bool
        If end gaps are free from penalty.
    eps : float, optional
        Absolute tolerance. Default is 1e-7.

    Returns
    -------
    float
        Optimal alignment score.
    ndarray of int of shape (k, 2)
        Coordinates of all (k) alignment stops.

    Notes
    -----
    Coordinates (col, row) are sorted in ascending order.

    """
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1

    # local alignment
    if local:
        max_ = scomat.max()
        if eps:
            tests = np.isclose(scomat, max_, rtol=0, atol=eps)
        else:
            tests = scomat == max_
        return max_, np.argwhere(tests)

    # semi-global alignment
    elif free_ends:
        max_ = np.max([scomat[:, n].max(), scomat[m, :].max()])
        if eps:
            test1 = np.isclose(scomat[:m, n], max_, rtol=0, atol=eps)
            test2 = np.isclose(scomat[m, :], max_, rtol=0, atol=eps)
        else:
            test1 = scomat[:m, n] == max_
            test2 = scomat[m, :] == max_
        ii = np.argwhere(test1)
        jj = np.argwhere(test2)
        return max_, np.vstack(
            (
                np.column_stack((ii, np.full(ii.shape[0], n))),
                np.column_stack((np.full(jj.shape[0], m), jj)),
            )
        )

    # global alignment
    else:
        return scomat[m, n], np.array([[m, n]])


# Encodings of moving directions during traceback. Columns are:
# 0. Row offset (i)
#   - Can be calculated with (state & 1) ^ 1
# 1. Column offset (j)
#   - Can be calculated with (state >> 1) ^ 1
# 2. Gap state
#   0. Substitution (no gap)
#   1. Insertion (gap in seq1)
#   2. Deletion (gap in seq2)
#   3. Invalid state
# 3. Matrix index
#   0. Main matrix
#   1. Insertion matrix (affine)
#   2. Deletion matrix (affine)
#   3. Insertion matrix 2 (2-piece affine)
#   4. Deletion matrix 2 (2-piece affine)
MOVES = np.array(
    [
        [1, 1, 0, 0],  # substitution
        [0, 1, 1, 0],  # insertion
        [1, 0, 2, 0],  # deletion
        [0, 1, 1, 1],  # extend insertion
        [1, 0, 2, 2],  # extend deletion
        [0, 0, 3, 1],  # jump to insertion matrix
        [0, 0, 3, 2],  # jump to deletion matrix
        [0, 1, 1, 3],  # extend insertion 2
        [1, 0, 2, 4],  # extend deletion 2
        [0, 0, 3, 3],  # jump to insertion matrix 2
        [0, 0, 3, 4],  # jump to deletion matrix 2
    ],
    dtype=int,
)


def _trailing_gaps(i, j, m, n, lengths, states):
    """Fill trailing gaps before traceback starts."""
    state = 3

    # bottom row: ends with insertions (gaps in seq1)
    if i == m and j < n:
        state = 1
        L = n - j

    # right-most column: ends with deletions (gaps in seq2)
    elif j == n and i < m:
        state = 2
        L = m - i

    # add a gap segment to path
    if state != 3:
        states.append(state)
        lengths.append(L)

    return state


def _leading_gaps(i, j, m, n, lengths, states):
    """Fill leading gaps after traceback ends."""
    state = 3

    # top row: starts with insertions (gaps in seq1)
    if i == 0 and j > 0:
        state = 1
        L = j

    # left-most column: starts with deletions (gaps in seq2)
    elif j == 0 and i > 0:
        state = 2
        L = i

    # add a gap segment to path
    if state != 3:
        # Note: The traceback algorithm terminates when reaching either edge, thus
        # guaranteeing that `state` cannot be equal to the previous state in the
        # path.
        states.append(state)
        lengths.append(L)


def _trace_one_linear_eq(lengths, states, i, j, matrices, gap, local, eps=1e-7):
    """Traceback across matrix body with linear gap penalty."""
    scomat = matrices[0]
    pos = len(lengths) - 1
    prev = states[pos] if pos >= 0 else 3

    # will stop when reaching either edge of the matrix
    while i and j:
        score = scomat[i, j]
        if local and score == 0:
            break

        # deletion (vertical; gap in seq2)
        if score == scomat[i - 1, j] - gap:
            state = 2
            i -= 1
        # insertion (horizontal; gap in seq1)
        elif score == scomat[i, j - 1] - gap:
            state = 1
            j -= 1
        # substitution (diagonal; no gap)
        else:
            state = 0
            i -= 1
            j -= 1

        # extend existing segment or create new segment
        if state == prev:
            lengths[pos] += 1
        else:
            lengths.append(1)
            states.append(state)
            prev = state
            pos += 1

    return i, j


def _trace_one_linear_tol(lengths, states, i, j, matrices, gap, local, eps=1e-7):
    """Traceback across matrix body with linear gap penalty."""
    scomat = matrices[0]
    pos = len(lengths) - 1
    prev = states[pos] if pos >= 0 else 3

    # will stop when reaching either edge of the matrix
    while i and j:
        score = scomat[i, j]
        if local and abs(score) <= eps:
            break
        gap_score = gap + score

        # deletion (vertical; gap in seq2)
        if abs(scomat[i - 1, j] - gap_score) <= eps:
            state = 2
            i -= 1
        # insertion (horizontal; gap in seq1)
        elif abs(scomat[i, j - 1] - gap_score) <= eps:
            state = 1
            j -= 1
        # substitution (diagonal; no gap)
        else:
            state = 0
            i -= 1
            j -= 1

        # extend existing segment or create new segment
        if state == prev:
            lengths[pos] += 1
        else:
            lengths.append(1)
            states.append(state)
            prev = state
            pos += 1

    return i, j


def _trace_one_affine_eq(lengths, states, i, j, matrices, gap, local, eps=1e-7):
    """Traceback across matrix body with linear gap penalty."""
    scomat, insmat, delmat = matrices[:3]
    pos = len(lengths) - 1
    prev = states[pos] if pos >= 0 else 3
    mat = 0

    while i and j:
        score = scomat[i, j]
        if local and score == 0:
            break
        state = 3

        # main matrix
        if mat == 0:
            if score == delmat[i, j]:  # jump to deletion matrix
                mat = 2
            elif score == insmat[i, j]:  # jump to insertion matrix
                mat = 1
            else:  # substitution (diagonal; no gap)
                state = 0
                i -= 1
                j -= 1

        # deletion matrix (vertical; gap in seq2)
        elif mat == 2:
            # open a new gap (jump back to main matrix) or extend an existing gap
            # (stay in the current matrix).
            if delmat[i, j] == scomat[i - 1, j] - gap:
                mat = 0
            state = 2
            i -= 1

        # insertion matrix (horizontal; gap in seq1)
        elif mat == 1:
            # same as above
            if insmat[i, j] == scomat[i, j - 1] - gap:
                mat = 0
            state = 1
            j -= 1

        # extend existing segment or create new segment
        if state != 3:
            if state == prev:
                lengths[pos] += 1
            else:
                lengths.append(1)
                states.append(state)
                prev = state
                pos += 1

    return i, j


def _trace_one_affine_tol(lengths, states, i, j, matrices, gap, local, eps=1e-7):
    """Traceback across matrix body with linear gap penalty."""
    scomat, insmat, delmat = matrices[:3]
    pos = len(lengths) - 1
    prev = states[pos] if pos >= 0 else 3
    mat = 0

    while i and j:
        score = scomat[i, j]
        if local and abs(score) <= eps:
            break
        state = 3

        # main matrix
        if mat == 0:
            if abs(delmat[i, j] - score) <= eps:  # jump to deletion matrix
                mat = 2
            elif abs(insmat[i, j] - score) <= eps:  # jump to insertion matrix
                mat = 1
            else:  # substitution (diagonal; no gap)
                state = 0
                i -= 1
                j -= 1

        # deletion matrix (vertical; gap in seq2)
        elif mat == 2:
            # open a new gap (jump back to main matrix) or extend an existing gap
            # (stay in the current matrix).
            if abs(scomat[i - 1, j] - gap - delmat[i, j]) <= eps:
                mat = 0
            state = 2
            i -= 1

        # insertion matrix (horizontal; gap in seq1)
        elif mat == 1:
            # same as above
            if abs(scomat[i, j - 1] - gap - insmat[i, j]) <= eps:
                mat = 0
            state = 1
            j -= 1

        # extend existing segment or create new segment
        if state != 3:
            if state == prev:
                lengths[pos] += 1
            else:
                lengths.append(1)
                states.append(state)
                prev = state
                pos += 1

    return i, j


def _make_path_obj(i, j, lengths, states):
    """Create a pairwise alignment path object."""
    return PairAlignPath(
        np.array(lengths[::-1], dtype=np.int64),
        np.array(states[::-1], dtype=np.uint8),
        np.array([i, j], dtype=np.int64),
    )


def _traceback_one(i, j, matrices, gap_open, gap_extend, local):
    """Traceback and return one optimal alignment path."""
    scomat = matrices[0]
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1
    gap = gap_open + gap_extend

    lengths = []
    states = []

    # fill trailing gaps
    if not local:
        _trailing_gaps(i, j, m, n, lengths, states)

    # traceback matrix body
    if gap_open:
        i, j = _trace_one_affine_tol(lengths, states, i, j, matrices, gap, local)
    else:
        i, j = _trace_one_linear_tol(lengths, states, i, j, matrices, gap, local)

    # fill leading gaps
    if not local:
        _leading_gaps(i, j, m, n, lengths, states)
        i, j = 0, 0

    return _make_path_obj(i, j, lengths, states)


def _linear_moves_eq(i, j, matrices, idx, submat, gap_oe, gap_extend, eps=1e-7):
    """Identify move direction(s) at a cell with linear gap penalty."""
    scomat = matrices[0]
    score = scomat[i, j]
    moves = []
    if score == scomat[i - 1, j - 1] + submat[i - 1, j - 1]:  # subsitution
        moves.append(0)
    if score == scomat[i, j - 1] - gap_extend:  # insertion
        moves.append(1)
    if score == scomat[i - 1, j] - gap_extend:  # deletion
        moves.append(2)
    return moves


def _linear_moves_tol(i, j, matrices, idx, submat, gap_oe, gap_extend, eps=1e-7):
    """Identify move direction(s) at a cell with linear gap penalty."""
    scomat = matrices[0]
    score = scomat[i, j]
    moves = []
    if abs(scomat[i - 1, j - 1] + submat[i - 1, j - 1] - score) <= eps:  # subsitution
        moves.append(0)
    if abs(scomat[i, j - 1] - gap_extend - score) <= eps:  # insertion
        moves.append(1)
    if abs(scomat[i - 1, j] - gap_extend - score) <= eps:  # deletion
        moves.append(2)
    return moves


def _affine_moves_eq(i, j, matrices, idx, submat, gap_oe, gap_extend, eps=1e-7):
    """Identify move direction(s) at a cell with affine gap penalty."""
    scomat, insmat, delmat = matrices[:3]
    moves = []

    # main matrix
    if idx == 0:
        score = scomat[i, j]
        if score == scomat[i - 1, j - 1] + submat[i - 1, j - 1]:  # subsitution
            moves.append(0)
        if score == insmat[i, j]:  # jump to insertion matrix
            moves.append(5)
        if score == delmat[i, j]:  # jump to deletion matrix
            moves.append(6)
    # insertion matrix
    elif idx == 1:
        score = insmat[i, j]
        if score == scomat[i, j - 1] - gap_oe:  # open insertion
            moves.append(1)
        if score == insmat[i, j - 1] - gap_extend:  # extend insertion
            moves.append(3)
    # deletion matrix
    else:
        score = delmat[i, j]
        if score == scomat[i - 1, j] - gap_oe:  # open deletion
            moves.append(2)
        if score == delmat[i - 1, j] - gap_extend:  # extend deletion
            moves.append(4)

    return moves


def _affine_moves_tol(i, j, matrices, idx, submat, gap_oe, gap_extend, eps=1e-7):
    """Identify move direction(s) at a cell with affine gap penalty."""
    scomat, insmat, delmat = matrices[:3]
    moves = []

    # main matrix
    if idx == 0:
        score = scomat[i, j]
        if (
            abs(scomat[i - 1, j - 1] + submat[i - 1, j - 1] - score) <= eps
        ):  # subsitution
            moves.append(0)
        if abs(insmat[i, j] - score) <= eps:  # jump to insertion matrix
            moves.append(5)
        if abs(delmat[i, j] - score) <= eps:  # jump to deletion matrix
            moves.append(6)
    # insertion matrix
    elif idx == 1:
        score = insmat[i, j]
        if abs(scomat[i, j - 1] - gap_oe - score) <= eps:  # open insertion
            moves.append(1)
        if abs(insmat[i, j - 1] - gap_extend - score) <= eps:  # extend insertion
            moves.append(3)
    # deletion matrix
    else:
        score = delmat[i, j]
        if abs(scomat[i - 1, j] - gap_oe - score) <= eps:  # open deletion
            moves.append(2)
        if abs(delmat[i - 1, j] - gap_extend - score) <= eps:  # extend deletion
            moves.append(4)

    return moves


def _traceback_all(
    stops, matrices, submat, gap_open, gap_extend, local, max_paths=None, eps=1e-7
):
    """Traceback and return all optimal alignment paths.

    max_paths must be >= 1 or None

    """
    scomat = matrices[0]
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1
    gap_oe = gap_open + gap_extend

    # Dispatch function that decides the moving direction at each cell.
    moves_func = _affine_moves_tol if gap_open else _linear_moves_tol

    # Results (optimal alignment paths)
    paths = []

    # Iterate over all starting points.
    for i, j in stops:
        lengths = []  # segment lengths
        states = []  # segment states (0 - substitution, 1 - insertion, 2 - deletion)

        prev_state = 3  # previous segment status (start with an impossible code 3)
        mat_idx = 0  # matrix index (start from main matrix (0))

        # Pre-populate the path with trailing gaps, if applicable.
        if not local:
            prev = _trailing_gaps(i, j, m, n, lengths, states)

        # Create a stack to store branching alignment paths and to enable depth-first
        # search (DFS).
        stack = [(i, j, mat_idx, lengths, states, prev_state)]

        # Note: The stack's maximum size is m + n, which is equivalent to the maximum
        # size of a full path, times (b - 1), where b is the number possible branches
        # per cell (3 for linear and affine). Regardless of the total number of full
        # paths, which could easily explode, the stack size is linear and manageable.

        while stack:
            i, j, mat_idx, lengths, states, prev_state = stack.pop()

            # Check whether a full path has been completed. There are two scenarios:
            finished = False

            # 1) Local alignment: Cell value is 0. Path is already completed.
            if local and scomat[i, j] == 0:  # TODO: floating point
                finished = True

            # 2) Global alignment: Reached the top or left-most edge of the matrix.
            # Path will extend straight left or up to the top-left cell.
            elif i == 0 or j == 0:
                if not local:
                    _leading_gaps(i, j, m, n, lengths, states)
                    i, j = 0, 0
                finished = True

            # Create a path object, and halt if the path number limit has been reached.
            if finished:
                paths.append(_make_path_obj(i, j, lengths, states))
                if max_paths and len(paths) == max_paths:
                    break
                else:
                    continue

            # Otherwise, the current cell must be within the main body of the matrix.
            # Check all possible moving directions and get ones that match the current
            # optimal alignment score.
            moves = moves_func(i, j, matrices, mat_idx, submat, gap_oe, gap_extend)

            # This is impossible. Raise error for debugging purpose.
            n_moves = len(moves)
            if n_moves == 0:
                raise ValueError("Traceback cannot proceed.")

            n_moves_1 = n_moves - 1
            for k in range(n_moves):  # TODO: avoid copy
                row = MOVES[moves[k]]

                # If path branches at this cell (i.e., 1 path becomes 2 or 3), we need
                # to create copies of the path. Otherwise, we can use the same path to
                # save memory and runtime.
                if k < n_moves_1:
                    lengths_ = lengths.copy()
                    states_ = states.copy()
                else:
                    lengths_ = lengths
                    states_ = states

                # Deal with gap state. 3 means jumping from main matrix into another
                # matrix without advancing the path, therefore reset the state to the
                # previous state. Otherwise (0, 1, 2), check if state is identical to
                # the previous state. If so, extend the current segment. Else, create
                # a new segment.
                state = row[2]
                if state == 3:
                    state = prev_state
                elif state == prev_state:
                    lengths_[-1] += 1
                else:
                    lengths_.append(1)
                    states_.append(state)

                stack.append((i - row[0], j - row[1], row[3], lengths_, states_, state))

    return paths
