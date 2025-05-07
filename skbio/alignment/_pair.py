# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numbers import Real

import numpy as np

from skbio._base import SkbioObject
from skbio.util._decorator import classonlymethod
from skbio.sequence import Sequence, SubstitutionMatrix


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


def _make_submat(seq1, seq2, match, mismatch):
    """Pre-compute a match/mismatch array to facilitate lookup."""
    # This is an alternative to mthmis[int(seq1[i, j] != seq2[i, j])]
    return np.where(seq1[:, None] == seq2[None, :], match, mismatch)


class PairAligner(SkbioObject):
    r"""Specify parameters for pairwise sequence alignment.

    Parameters
    ----------
    mode : str, optional
        Alignment mode. Can be "global" or "local".
    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Score of a match, mismatch or substitution.
    gap_score : float or tuple of (float, float)
        Score of a gap. The value is usually negative, representing a penalty from
        the alignment score.
    free_ends : bool, optional
        If True (default), gaps at sequence ends are not penalized. Otherwise False.

    See Also
    --------
    PairAlignPath

    """

    def __init__(
        self,
        mode="global",
        sub_score=(1, -1),
        gap_score=-2,
        free_ends=True,
    ):
        # alignment mode
        if mode not in ("global", "local"):
            raise ValueError("`mode` must be either 'global', 'local'.")
        self._local = mode == "local"

        # substitution matrix or match/mismatch scores
        if isinstance(sub_score, str):
            sub_score = SubstitutionMatrix.by_name(sub_score)
        if isinstance(sub_score, SubstitutionMatrix):
            self._submat = sub_score._data
            self._subidx = sub_score._char_hash
            self._mthmis = np.empty(0)
        else:
            self._submat = np.empty((0, 0))
            self._subidx = None
            self._mthmis = np.array(sub_score, dtype=float)

        # affine or linear gap penalty
        if isinstance(gap_score, Real):
            self._gap_open, self._gap_extend = 0, gap_score
            self._affine = False
        else:
            self._gap_open, self._gap_extend = gap_score
            self._affine = True

        # terminal gap policy (for global alignment)
        self._free_ends = free_ends

        # input sequences (integer codes) and lengths
        self._seq1 = None
        self._seq2 = None
        self._len1 = None
        self._len2 = None

        # alignment matrices
        self._alnmat = None
        self._insmat = None
        self._delmat = None

        self._score = None
        self._starts = None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return (
            f"Mode: {self.mode},"
            f"sub_score: {self.sub_score},"
            f"gap_score: {self.gap_score},"
        )

    @property
    def mode(self):
        return "local" if self._local else "global"

    @property
    def sub_score(self):
        """Substitution scoring policy."""
        # TODO: polish
        return tuple(map(float, self._mthmis))

    @property
    def gap_score(self):
        """Gap penalizing policy."""
        return (self._gap_open, self._gap_extend) if self._affine else self._gap_extend

    @property
    def free_ends(self):
        """Terminal gap penalizing policy."""
        return self._free_ends

    @property
    def score(self):
        """Alignment score."""
        return float(self._score)

    @property
    def matrix(self):
        """Alignment matrix."""
        return self._alnmat

    def align(self, seq1, seq2):
        r"""Compute an alignment matrix of two sequences."""
        self._seq1 = _seq_to_bytes(seq1)
        self._seq2 = _seq_to_bytes(seq2)
        self._len1 = self._seq1.size
        self._len2 = self._seq2.size

        # initialize alignment matrix
        if self._affine:
            self._alloc_affine_matrix()
            if self._local:
                self._init_local_affine_matrix()
            else:
                self._init_global_affine_matrix()
        else:
            self._alloc_linear_matrix()
            if self._local:
                self._init_local_linear_matrix()
            else:
                self._init_global_linear_matrix()

        # fill alignment matrix (compute-intensive)
        if self._affine:
            self._fill_affine_matrix()
        else:
            self._fill_linear_matrix()

        # find optimal alignment score
        if self._affine:
            pass
        else:
            if self._local:
                pass
            else:
                self._score_global_linear_endcost()

    @property
    def path(self):
        r"""Return the first optimal alignment path."""
        if self._affine:
            pass
        else:
            if self._local:
                pass
            else:
                return self._trace_global_linear_one()

    @property
    def paths(self):
        r"""Iterate over all optimal alignment paths."""
        if self._affine:
            pass
        else:
            if self._local:
                pass
            else:
                return self._trace_global_linear_iter()

    ### Step 1: Allocate alignment matrix ###

    def _alloc_linear_matrix(self):
        """Allocate alignment matrix with linear gap penalty."""
        # Note: The array should be C-contiguous to facilitate row-wise iteration.
        # NumPy's default array is already C-contiguous. This is also enforced by the
        # unit test.
        self._alnmat = np.empty((self._len1 + 1, self._len2 + 1))
        self._insmat = None
        self._delmat = None

    def _alloc_affine_matrix(self):
        """Allocate alignment matrices with affine gap penalty."""
        shape_ = (self._len1 + 1, self._len2 + 1)
        self._alnmat = np.empty(shape_)
        self._insmat = np.empty(shape_)
        self._delmat = np.empty(shape_)

    ### Step 2: Initialize alignment matrix ###

    def _init_global_linear_matrix(self):
        """Initialize global alignment matrix with linear gap penalty."""
        m1, n1 = self._len1 + 1, self._len2 + 1
        gap = self._gap_extend
        alnmat = self._alnmat
        alnmat[0, 0] = 0
        alnmat[1:, 0] = np.arange(1, m1) * gap
        alnmat[0, 1:] = np.arange(1, n1) * gap

    def _init_local_linear_matrix(self):
        """Initialize local alignment matrix with linear gap penalty."""
        alnmat = self._alnmat
        alnmat[0, 0] = 0
        alnmat[1:, 0] = 0
        alnmat[0, 1:] = 0

    def _init_global_affine_matrix(self):
        """Initiate global alignment matrices with affine gap penalty."""
        pass

    def _init_local_affine_matrix(self):
        """Initiate local alignment matrices with affine gap penalty."""
        pass

    ### Step 3: Fill alignment matrix (compute-intensive) ###

    def _fill_linear_matrix(self):
        """Fill alignment matrix with linear gap penalty."""
        alnmat = self._alnmat
        m1, n1 = self._len1 + 1, self._len2 + 1
        args = (alnmat, self._seq1, self._seq2, self._mthmis, self._gap_extend)
        # TODO: optimize
        base = 0 if self._local else -np.inf
        for i in range(1, m1):
            for j in range(1, n1):
                alnmat[i, j] = max(base, *_moves(i, j, *args))

    def _fill_affine_matrix(self):
        """Fill alignment matrix with affine gap penalty."""
        pass

    ### Step 4: Find optimal score and traceback starting points ###

    def _score_global_linear_endcost(self):
        """Score global alignment matrix with linear gap penalty."""
        i, j = self._len1, self._len2
        self._score = self._alnmat[i, j]
        self._starts = [(i, j)]

    ### Step 5: Traceback alignment matrix to find optimal paths ###

    def _trace_global_linear_one(self):
        """Traceback global matrix with linear gap penalty and return one path."""
        seq1, seq2 = self._seq1, self._seq2
        m, n = self._len1, self._len2
        alnmat = self._alnmat
        mthmis = self._mthmis
        gap = self._gap_extend

        pos = m + n
        path = np.empty(pos, dtype=np.uint8)

        i, j = self._starts[0]

        while True:
            if i > 0:
                if j > 0:
                    curr = alnmat[i, j]
                    moves = _moves(i, j, alnmat, seq1, seq2, mthmis, gap)
                    pos -= 1
                    if moves[2] == curr:
                        path[pos] = 2
                        i -= 1
                    elif moves[1] == curr:
                        path[pos] = 1
                        j -= 1
                    elif moves[0] == curr:  # can be omitted
                        path[pos] = 0
                        i -= 1
                        j -= 1
                    else:
                        raise ValueError("Floating-point error encountered.")

                else:  # left-most column
                    new_pos = pos - i
                    path[new_pos:pos] = 2
                    pos = new_pos
                    break
            else:
                if j > 0:  # top row
                    new_pos = pos - j
                    path[new_pos:pos] = 1
                    pos = new_pos
                    break
                else:  # top-left cell
                    break

        return path[pos:]

    def _trace_global_linear_iter(self):
        """Traceback global matrix with linear gap penalty and return an iterator."""
        seq1, seq2 = self._seq1, self._seq2
        m, n = self._len1, self._len2
        alnmat = self._alnmat
        mthmis = self._mthmis
        gap = self._gap_extend
        mn = m + n

        # Pre-allocation alignment path. It is a vector of 0 (substitution), 1
        # (insertion), and 2 (deletion). The maximum size of is the sum of the two
        # sequences. It will be filled in reverse order.
        path = np.empty(mn, dtype=np.uint8)

        # Iterate over all starting points.
        for i, j in self._starts:
            # Stack to store branching alignment paths and to enable depth-first search
            # (DFS). Each tuple contains: (row index, column index, path, and the index
            # in the path).
            stack = [(i, j, path, mn)]
            while stack:
                i, j, path, pos = stack.pop()
                if i > 0:
                    # Current cell is within the main body of the matrix. Will check
                    # all three possible directions and save the optimal ones.
                    if j > 0:
                        # Whether path branches at this point (i.e., 1 path becomes 2
                        # or 3). In True, we need to create copies of the path.
                        branching = False

                        moves = _moves(i, j, alnmat, seq1, seq2, mthmis, gap)
                        curr = alnmat[i, j]

                        for k in range(3):
                            if moves[k] == curr:
                                if not branching:
                                    new_path = path
                                    branching = True
                                else:
                                    new_path = path.copy()
                                new_pos = pos - 1
                                new_path[new_pos] = k

                                # Use bitwise operators to calculate moves.
                                # 0: i-1, j-1, 1: j-1, 2: i-1
                                # TODO: This can be further optimized by pre-creating
                                # a lookup table
                                stack.append(
                                    (
                                        i - ((k & 1) ^ 1),
                                        j - ((k >> 1) ^ 1),
                                        new_path,
                                        new_pos,
                                    )
                                )

                    # Reached the left-most column. Will move straight up to the top-
                    # left cell and create leading deletions.
                    else:
                        new_pos = pos - i
                        path[new_pos:pos] = 2
                        stack.append((0, 0, path, new_pos))
                        continue

                else:
                    # Reached the top row. Will move straight left and make leading
                    # insertions.
                    if j > 0:
                        new_pos = pos - j
                        path[new_pos:pos] = 1
                        stack.append((0, 0, path, new_pos))
                        continue

                    # Reached top-left cell. Will finish the current path and yield.
                    else:
                        yield path[pos:].copy()
                        continue
