# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from math import prod
from typing import TYPE_CHECKING

import numpy as np

from scipy.cluster.hierarchy import linkage
from skbio.sequence import Sequence, GrammaredSequence
from skbio.tree import TreeNode
from ._path import AlignPath
from ._pair import (
    pair_align,
    _init_matrices,
    _one_stop,
    _trailing_gaps,
    _leading_gaps,
    _encode_path,
)
from ._utils import encode_sequences, prep_gapcost, _get_seqids
from skbio.tree._utils import _tree_to_lnkmat
from ._cutils import (
    _fill_linear_matrix,
    _fill_affine_matrices,
    _trace_one_linear,
    _trace_one_affine,
    _fill_linear_rows,
    _fill_affine_rows,
    _trace_linear_rows,
    _trace_affine_rows,
)


if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Iterable
    from skbio.sequence import SubstitutionMatrix
    from ._utils import SequenceLike


def multi_align(
    sequences: Iterable[SequenceLike],
    /,
    sub_score: tuple[float, float] | SubstitutionMatrix | str = (1.0, -1.0),
    gap_cost: float | tuple[float, float] = 2.0,
    free_ends: bool = True,
    guide_tree: TreeNode | None = None,
    ids: Iterable[str] | None = None,
    method: str = "full",
    atol: float = 1e-5,
) -> AlignPath:
    r"""Align multiple sequences by progressive profile merging.

    Return a full-coverage alignment using a supplied guide tree or an automatically
    constructed UPGMA tree. Both linear and affine gap penalties are supported.

    .. versionadded:: 0.7.4

    Parameters
    ----------
    sequences : iterable of sequence-like
        Two or more nonempty, ungapped sequences of the same type. Supports the
        sequence types accepted by :func:`pair_align`, including scikit-bio sequences,
        strings, and sequences of tokens or numbers. Decode byte strings first.
        Gaps in grammared
        sequences, and '-' or '.' in other string-like inputs, are not allowed.
        Duplicate sequences are retained as separate rows.
    sub_score : tuple of (float, float), SubstitutionMatrix, or str, optional
        Match/mismatch scores, a substitution matrix, or its name. Scores must be
        finite, and the matrix must be symmetric. Default is (1.0, -1.0).
    gap_cost : float or tuple of (float, float), optional
        Nonnegative, finite gap penalties. A scalar ``g`` gives a cost of ``g * k``
        for a run of ``k`` newly inserted profile columns. A tuple ``(o, e)`` gives
        ``o + e * k``: the first column costs ``o + e``. Default is 2.0.
    free_ends : bool, optional
        If True (default), insertions before or after an entire child profile are
        free. All residues and terminal columns remain in the output. Boundaries of
        individual sequences within a profile do not determine this exemption.
    guide_tree : TreeNode, optional
        Rooted binary guide tree. Its tip names must match ``ids`` exactly. The node
        supplied is treated as the root; branch lengths are ignored. The tree is not
        modified. If provided, pairwise distance calculation is skipped entirely.
        Otherwise, all sequence pairs are aligned and their score-based distances
        are passed in condensed form to SciPy average linkage (UPGMA).
    ids : iterable of str, optional
        Unique string identifiers in input order to match tips in the provided guide
        tree. Override sequence metadata if provided. Otherwise, use metadata ``'id'``
        values if present in every sequence and unique, or use ``['0', '1', ...]`` if
        none is present. Partial, duplicate, or non-string metadata IDs raise an error.
    method : {'full', 'rolling'}, optional
        Profile DP backend. 'full' (default) reuses the pairwise fill and traceback
        kernels, retaining one linear or three affine score matrices. 'rolling'
        retains two score rows per state and one byte of traceback per cell.
        Both use the same recurrence, endpoint selection, and traceback tolerance.
        A per-call workspace grows as needed and is reused across guide pairwise
        alignments and profile merges. Guide alignments always use full matrices.
    atol : float, optional
        Nonnegative finite absolute tolerance for traceback score comparisons,
        following :func:`pair_align`. Default is 1e-5. Set to zero for exact
        comparisons. Positive tolerance can select a slightly lower-scoring path;
        the DP maximum itself is unaffected. Also used for guide pairwise alignments.

    Returns
    -------
    AlignPath
        One alignment path in original input sequence order, starting at position
        zero and consuming every sequence in full.

    Raises
    ------
    ValueError
        If inputs are empty or gapped, scores or costs are invalid, identifiers do
        not match, or a guide tree is not binary with exactly the required tips.
    ValueError
        If a pairwise score cannot be converted to a finite nonnegative distance.
        A supplied guide tree permits alignment without this distance requirement.
    TypeError
        If ``free_ends`` is not Boolean or ``guide_tree`` is not a TreeNode.

    See Also
    --------
    pair_align
    align_score
    AlignPath
    skbio.tree.upgma

    Notes
    -----
    This algorithm follows the profile-averaging modification of Feng and Doolittle
    [1]_, [2]_, with a fixed guide tree and configurable linear/affine profile-gap
    costs. It does not reproduce their complete sequence-order refinement procedure.

    Guide merges are stored as pairs of cluster indices. Automatic guides use
    SciPy linkage order; supplied trees use postorder with their child order
    preserved. Branch lengths are not needed for alignment.

    At each internal tree node, two child alignments are merged without changing
    their existing residue relationships. For profiles :math:`A` and :math:`B` with
    :math:`r` and :math:`s` rows,
    the score of matching columns :math:`i` and :math:`j` is

    .. math::

        C_{ij} = \frac{1}{rs}\sum_{a\in A}\sum_{b\in B}\widetilde M(a_i,b_j).

    Here, :math:`\widetilde M` is the substitution matrix extended with a private gap
    symbol that scores zero against every residue and itself. Equivalently, column
    frequencies, **including gaps in their denominators**, give
    :math:`C_{ij}=f_{A,i}^T M f_{B,j}` using only nongap entries. Frequencies are not
    renormalized after excluding gaps. Each input row has equal weight, including
    duplicates; merged frequencies are weighted by child row counts.

    A merge path consumes both child columns (D), only an A column while inserting
    gaps into B (X), or only a B column while inserting gaps into A (Y). Its objective
    is the sum of D scores minus the costs of maximal consecutive X or Y runs. A
    new run is charged **once per profile**, without multiplying by row count or
    adjusting for residue occupancy. Existing gaps incur no additional cost in a D
    move. A D move ends a run even if individual rows contain gaps there; switching
    directly between X and Y is allowed and opens a new run.

    Both backends use the pairwise affine states: an overall best score H and
    insertion/deletion scores I/D. A gap opens from H at cost :math:`o+e` or extends
    at cost :math:`e`. Leading boundary scores are initialized in Python; free
    trailing gaps are appended after selecting the best last-row/last-column stop.
    Thus a completely nonoverlapping alignment can have score zero. Stop ties
    prefer the smallest (row, column); traceback prefers deletion, insertion, then
    diagonal, and gap extension before opening, as in :func:`pair_align`.

    Profile-column substitution scores require a floating-point matrix of size
    :math:`mn` in both backends. Beyond this shared cost, 'full' retains one or
    three floating-point matrices, whereas 'rolling' uses one byte per cell and
    linear-sized score buffers. Flat NumPy work buffers double in capacity when
    necessary; active prefixes are reshaped into contiguous views. Retained child
    profiles and returned paths own their data and never alias this workspace.

    Each merge optimizes this objective, but the overall progressive alignment is
    heuristic. It does **not** optimize the induced-pair sum-of-pairs (SP) score:
    :func:`align_score` penalizes resulting residue-gap runs separately for every
    sequence pair, including gaps treated as neutral in later profile merges. Its
    terminal boundaries also refer to individual pairs. Use it for final evaluation,
    not as the internal merge objective. Exact affine-SP merging is a different,
    more demanding optimization problem [3]_.

    For automatic guide construction, one global pairwise alignment per pair gives
    score :math:`S_{ab}`, aligned length :math:`L` (including terminal columns), and
    charged gap cost :math:`K`. With original residue counts :math:`c_a,c_b`, the
    distance is

    .. math::

        S_{\mathrm{rand}} = c_a^T M c_b/L - K,\quad
        S_{\mathrm{max}} = (S_{aa}+S_{bb})/2,\quad
        d = -\ln\left(\frac{S_{ab}-S_{\mathrm{rand}}}
                            {S_{\mathrm{max}}-S_{\mathrm{rand}}}\right).

    Self-scores are optimal pairwise self-alignment scores. :math:`K` uses the selected
    gap and terminal settings. The self-score average is not a guaranteed upper
    bound for arbitrary scoring schemes. This extends the constant-gap random-score
    formula in [2]_ and is not a guaranteed mathematical metric. Every sequence
    pair is evaluated, including duplicates. Identical sequences can also have an
    undefined normalization (for example, identical homopolymers).
    Nonpositive numerator/denominator or a ratio above one
    (beyond eight machine epsilons of the scoring dtype) raises an error, rather
    than inventing a distance for unrelated sequences. Pairwise alignments are
    processed one at a time using shared score buffers; self-alignments omit
    traceback. DP and profile arithmetic use the substitution-matrix dtype; FD
    normalization uses float64 to limit additional rounding near the random baseline.
    The eight-epsilon allowance is a small boundary-roundoff policy, not a bound on
    accumulated alignment-score error.
    Two inputs bypass distances and use pair_align directly, regardless of `method`.
    Tied pairwise alignments can produce different distance statistics; tied guides
    and different floating-point arithmetic can consequently change the final MSA.

    All-pairs dynamic programming takes approximately :math:`O(n^2 L^2)` time for
    :math:`n` sequences of comparable length :math:`L`. Each profile merge requires
    quadratic time and memory in child alignment lengths, in addition to computing
    column scores. This method
    targets moderate research and educational use, not very large sequence sets.

    References
    ----------
    .. [1] Feng, D.-F. and Doolittle, R. F. (1987). Progressive sequence alignment as
       a prerequisite to correct phylogenetic trees. J. Mol. Evol. 25, 351-360.
       doi:10.1007/BF02603120. See the Note Added in Proof, p. 359.
    .. [2] Feng, D.-F. and Doolittle, R. F. (1996). Progressive alignment of amino acid
       sequences and construction of phylogenetic trees from them. Methods Enzymol.
       266, 368-382. doi:10.1016/S0076-6879(96)66023-6.
    .. [3] Wheeler, T. J. and Kececioglu, J. D. (2007). Multiple alignment by aligning
       alignments. Bioinformatics 23, i559-i568. doi:10.1093/bioinformatics/btm226.

    Examples
    --------
    Align three sequences and materialize their aligned strings:

    >>> from skbio.alignment import multi_align, align_score
    >>> seqs = ['ACGT', 'AGT', 'ACGT']
    >>> path = multi_align(seqs, free_ends=False)
    >>> path.to_aligned(seqs)
    ['ACGT', 'A-GT', 'ACGT']
    >>> align_score((path, seqs), free_ends=False)
    6.0

    A supplied tree controls merge order and avoids automatic distances:

    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(['((a,b),c);'])
    >>> path = multi_align(seqs, ids=['a', 'b', 'c'], guide_tree=tree,
    ...                    gap_cost=(2, 1), free_ends=False)
    >>> path.shape[0]
    3

    Sequence metadata IDs, including those read from FASTA headers, are inferred:

    >>> from skbio import DNA
    >>> named = [DNA(seq, metadata={'id': name})
    ...          for seq, name in zip(seqs, ['a', 'b', 'c'])]
    >>> multi_align(named, guide_tree=tree, free_ends=False).to_aligned(named)
    ['ACGT', 'A-GT', 'ACGT']

    """
    sequences = list(sequences)
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required.")
    if method not in ("full", "rolling"):
        raise ValueError("`method` must be 'full' or 'rolling'.")
    if not np.isscalar(atol) or not np.isfinite(atol) or atol < 0:
        raise ValueError("`atol` must be finite and nonnegative.")
    if not isinstance(free_ends, (bool, np.bool_)):
        raise TypeError("`free_ends` must be Boolean.")
    for seq in sequences:
        if isinstance(seq, bytes):
            raise TypeError("Byte strings must be decoded before alignment.")
        if isinstance(seq, GrammaredSequence):
            gapped = seq.has_gaps()
        elif isinstance(seq, (Sequence, str)):
            chars = str(seq)
            gapped = "-" in chars or "." in chars  # TODO: is this correct?
        else:
            gapped = False
        if gapped:
            raise ValueError("Input sequences must be ungapped.")

    encoded, matrix, _ = encode_sequences(sequences, sub_score)
    # Check symmetry. It's not that MSA absolutely require it. But its behavior will be
    # unpredictable with asymmetric substitution scores.
    if not np.isfinite(matrix).all() or not np.array_equal(matrix, matrix.T):
        raise ValueError("Substitution scores must be finite and symmetric.")

    atol = matrix.dtype.type(atol)
    if not np.isfinite(atol):
        raise ValueError("`atol` is too large for the scoring dtype.")

    # Prepare gap penalties
    gap_open, gap_extend = prep_gapcost(gap_cost, matrix.dtype.type)
    if not np.isfinite([gap_open, gap_extend]).all() or min(gap_open, gap_extend) < 0:
        raise ValueError("Gap costs must be finite and nonnegative.")

    # Prepare custom guide tree
    # TODO: Skip if tree is not supplied
    ids = _get_seqids(sequences, ids)

    # Fall back to pairwise alignment
    # TODO: Merge into main workflow.
    if len(sequences) == 2:
        path = pair_align(
            *sequences,
            sub_score=sub_score,
            gap_cost=gap_cost,
            free_ends=free_ends,
            atol=atol,
        ).paths[0]
        return AlignPath.from_bits(path.to_bits())

    # Shrink substitution matrix to observed characters only. This accelerates the
    # calculation without changing the result.
    # TODO: investigate
    alphabet, inverse = np.unique(np.concatenate(encoded), return_inverse=True)
    encoded = np.split(inverse, np.cumsum([len(x) for x in encoded])[:-1])
    matrix = np.ascontiguousarray(matrix[np.ix_(alphabet, alphabet)])

    workspace = _ProfileWorkspace()
    n = len(sequences)

    # Decide merging order. If a guide tree is provided, convert it into a linkage
    # matrix (only first two columns are needed). Otherwise, perform pairwise
    # alignments, calculate a distance matrix using the Feng-Doolittle metric, then
    # calculate a guide tree using UPGMA and retain the linkage matrix.
    # `merges` is an index array of (n_seqs - 1, 2)
    if guide_tree is None:
        distances = _multi_distances(
            encoded, matrix, gap_open, gap_extend, free_ends, ids, workspace, atol
        )
        merges = linkage(distances, method="average")[:, :2].astype(np.intp)
    else:
        merges = _tree_to_lnkmat(guide_tree, ids)

    profiles = {}
    eye = np.eye(len(matrix), dtype=matrix.dtype)
    for parent, children in enumerate(merges, n):
        # Materialize leaves when first used and release children after each merge.
        for child in children:
            if child < n:
                profiles[child] = (
                    eye[encoded[child]],
                    np.zeros((1, len(encoded[child])), dtype=bool),
                    [child],
                )
        a, b = (profiles.pop(child) for child in children)
        profiles[parent] = _merge_profiles(
            a, b, matrix, gap_open, gap_extend, free_ends, workspace, method, atol
        )
    _, bits, order = profiles[2 * n - 2]
    return AlignPath.from_bits(bits[np.argsort(order)])


def _multi_distances(
    encoded, matrix, gap_open, gap_extend, free_ends, ids, workspace, atol
):
    """Stream pairwise alignments into SciPy's condensed distance ordering."""
    n = len(encoded)
    # DP uses the matrix dtype. FD arithmetic uses float64 to avoid additional
    # rounding in composition sums and subtraction near the random baseline;
    # widening cannot recover precision already lost during alignment.
    counts = np.array(
        [np.bincount(x, minlength=len(matrix)) for x in encoded], dtype=np.float64
    )
    weighted = counts @ matrix.astype(np.float64)
    fd_open, fd_extend = float(gap_open), float(gap_extend)
    self_scores = np.empty(n, dtype=np.float64)
    data = np.empty(n * (n - 1) // 2, dtype=np.float64)
    # Small allowance at the dimensionless boundary one, not a general DP error
    # bound or a tolerance for declaring alignment scores tied.
    allowance = 8 * np.finfo(matrix.dtype).eps
    args = (gap_open, gap_extend, free_ends, workspace, atol)
    # Upper-triangle traversal needs future sequences' self-scores in advance.
    # Recreate query tables in the pair loop instead of retaining all of them.
    for i, seq in enumerate(encoded):
        _, self_scores[i] = _align_pair(matrix[seq], seq, *args, traceback=False)
    position = 0
    for i in range(n - 1):
        seq = encoded[i]
        query = matrix[seq]
        for j in range(i + 1, n):
            moves, score = _align_pair(query, encoded[j], *args)
            path = _encode_path(moves, 0, len(seq), 0, len(encoded[j]))
            s_max = (self_scores[i] + self_scores[j]) / 2
            try:
                distance = _fd_dist(
                    score,
                    path,
                    weighted[i] @ counts[j],
                    s_max,
                    fd_open,
                    fd_extend,
                    free_ends,
                    allowance,
                )
            except ValueError as e:
                raise ValueError(
                    f"Cannot calculate a distance between {ids[i]!r} and {ids[j]!r}: "
                    f"{e}; supply a guide tree instead."
                ) from e
            # Consecutive pairs already follow SciPy's condensed ordering.
            data[position] = distance
            position += 1
    return data


def _fd_dist(score, path, expected, s_max, gap_open, gap_extend, free_ends, allowance):
    """Convert an alignment score to a Feng-Doolittle distance.

    `expected` is the unnormalized composition sum c_a.T @ M @ c_b; `s_max`
    is the average optimal self-score. The path consumes both sequences in full.
    Existing validation guarantees finite scoring parameters and a nonempty path.
    """
    bits = path.to_bits(expand=False).astype(bool)
    lengths = path.lengths
    charged = bits.any(axis=0)
    if free_ends:
        # A path segment is a maximal run with a fixed gap state. Terminal gaps
        # have consumed either none or all of the residues in the gapped row.
        positions = np.cumsum(~bits * lengths, axis=1)
        terminal = ((positions == 0) | (positions == positions[:, -1:])) & bits
        charged &= ~terminal.any(axis=0)
    cost = gap_open * charged.sum() + gap_extend * lengths[charged].sum()
    s_rand = expected / lengths.sum() - cost
    numerator = score - s_rand
    denominator = s_max - s_rand
    if not np.isfinite([numerator, denominator]).all():
        raise ValueError("nonfinite score normalization")
    if denominator <= 0:
        raise ValueError("self-score normalization is not positive")
    if numerator <= 0:
        raise ValueError("alignment score does not exceed the random baseline")
    ratio = numerator / denominator
    if ratio > 1 + allowance:
        raise ValueError("alignment score exceeds the self-score normalization")
    return -np.log(min(ratio, 1.0))


def _merge_profiles(
    a, b, matrix, gap_open, gap_extend, free_ends, workspace=None, method="full", atol=0
):
    """Merge residue counts and gap masks; original row order travels with them."""
    counts_a, bits_a, order_a = a
    counts_b, bits_b, order_b = b
    if workspace is None:
        workspace = _ProfileWorkspace()
    scores = workspace.get("scores", (len(counts_a), len(counts_b)), matrix.dtype)
    np.matmul(
        (counts_a / len(order_a)) @ matrix, (counts_b / len(order_b)).T, out=scores
    )
    indices, _ = _align_profiles(
        scores, gap_open, gap_extend, free_ends, workspace, method, atol
    )
    width = indices.shape[1]
    counts = np.zeros((width, len(matrix)), dtype=matrix.dtype)
    bits = []
    for old_counts, old_bits, take in (
        (counts_a, bits_a, indices[0]),
        (counts_b, bits_b, indices[1]),
    ):
        present = take >= 0
        counts[present] += old_counts[take[present]]
        new_bits = np.ones((len(old_bits), width), dtype=bool)
        new_bits[:, present] = old_bits[:, take[present]]
        bits.append(new_bits)
    return counts, np.concatenate(bits), order_a + order_b


class _ProfileWorkspace:
    """Per-call flat buffers shared by pairwise alignment and profile merging."""

    def __init__(self):
        self.buffers = {}

    def get(self, name, shape, dtype):
        size = prod(shape)
        old = self.buffers.get(name)
        if old is None or old.size < size or old.dtype != dtype:
            capacity = size if old is None else max(size, 2 * old.size)
            self.buffers[name] = np.empty(capacity, dtype=dtype)
        return self.buffers[name][:size].reshape(shape)


def _align_pair(
    query, target, gap_open, gap_extend, free_ends, workspace, atol, traceback=True
):
    """Align an encoded query/target with prepared costs in the query dtype.

    Return dense moves borrowed from the workspace (or None) and a Python-float
    score. Consume moves before the next call; score-only calls skip traceback.
    """
    m, n = len(query), len(target)
    affine = gap_open != 0
    matrices = tuple(
        workspace.get(f"dp{k}", (m + 1, n + 1), query.dtype)
        for k in range(3 if affine else 1)
    )
    # Reshaping changes row stride, so smaller matrices also need fresh boundaries.
    _init_matrices(matrices, gap_open, gap_extend, False, free_ends, free_ends)
    if affine:
        _fill_affine_matrices(*matrices, query, target, gap_open, gap_extend, False)
    else:
        _fill_linear_matrix(matrices[0], query, target, gap_extend, False)
    score, stops = _one_stop(matrices[0], False, free_ends, free_ends)
    if not traceback:
        return None, float(score)
    i, j = stops[0]
    path = workspace.get("path", (m + n,), np.uint8)
    pos, _, _ = _trailing_gaps(path, m + n, i, j, m, n, True, True)
    if affine:
        pos, i, j = _trace_one_affine(
            path, pos, i, j, *matrices, gap_extend, False, atol
        )
    else:
        pos, i, j = _trace_one_linear(
            path, pos, i, j, matrices[0], gap_extend, False, atol
        )
    pos, _, _ = _leading_gaps(path, pos, i, j, True, True)
    return path[pos:], float(score)


def _align_pair_roll(scores, gap_open, gap_extend, free_ends, workspace, atol):
    """Align a dense column-score table using rolling rows and stored traceback.

    Like `_align_pair`, return workspace-borrowed dense moves and a Python-float
    score. The rolling kernels consume column scores directly, without a target
    index array. Costs and tolerance are already in the scoring dtype.
    """
    m, n = scores.shape
    dtype = scores.dtype
    affine = gap_open != 0
    matrices = tuple(
        workspace.get(f"dp{k}", (2, n + 1), dtype) for k in range(3 if affine else 1)
    )
    _init_matrices(matrices, gap_open, gap_extend, False, free_ends, free_ends)
    edge = workspace.get("edge", (m + 1,), dtype)
    if free_ends:
        edge[:] = 0
    else:
        edge[0] = 0
        edge[1:] = np.arange(1, m + 1, dtype=dtype)
        edge[1:] *= -gap_extend
        if gap_open:
            edge[1:] -= gap_open
    trace = workspace.get("trace", (m + 1, n + 1), np.uint8)
    if affine:
        _fill_affine_rows(*matrices, edge, trace, scores, gap_open, gap_extend, atol)
    else:
        _fill_linear_rows(matrices[0], edge, trace, scores, gap_extend, atol)
    last = matrices[0][m % 2]
    i, j = m, n
    if free_ends:
        # Same endpoint ordering as _one_stop: exclude (m,n) from column scan.
        i, j = int(edge[:m].argmax()), int(last.argmax())
        if edge[i] >= last[j]:
            j = n
        else:
            i = m
    score = edge[i] if j == n else last[j]
    path = workspace.get("path", (m + n,), np.uint8)
    pos, _, _ = _trailing_gaps(path, m + n, i, j, m, n, True, True)
    if affine:
        pos, i, j = _trace_affine_rows(path, pos, i, j, trace)
    else:
        pos, i, j = _trace_linear_rows(path, pos, i, j, trace)
    pos, _, _ = _leading_gaps(path, pos, i, j, True, True)
    moves = path[pos:]
    return moves, float(score)


def _align_profiles(
    scores, gap_open, gap_extend, free_ends, workspace=None, method="full", atol=0
):
    """Align profile columns with costs already prepared in the scoring dtype.

    Return independent column indices and the DP maximum. With positive `atol`,
    the traced path can score below that maximum, just as in pair_align.
    """
    if workspace is None:
        workspace = _ProfileWorkspace()
    n = scores.shape[1]
    if method == "full":
        # A column-score table is a pairwise query table with target indices 0..n-1.
        target = workspace.get("target", (n,), np.intp)
        target[:] = np.arange(n)
        moves, score = _align_pair(
            scores, target, gap_open, gap_extend, free_ends, workspace, atol
        )
    else:
        moves, score = _align_pair_roll(
            scores, gap_open, gap_extend, free_ends, workspace, atol
        )
    present = np.array([moves != 1, moves != 2])
    indices = np.cumsum(present, axis=1) - 1
    indices[~present] = -1
    return indices, float(score)
