# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from typing import Any, NamedTuple, TYPE_CHECKING

import numpy as np

from skbio.sequence import Sequence, GrammaredSequence
from skbio.alignment import AlignPath
from skbio.tree import TreeNode
from skbio.stats.distance import DistanceMatrix
from skbio.tree._utils import _tree_to_lnkmat
from skbio.util._array import ArrayWorkspace
from ._pair import (
    _init_matrices,
    _one_stop,
    _trailing_gaps,
    _leading_gaps,
    _encode_path,
)
from ._utils import encode_sequences, prep_gapcost, _get_seqids, _check_atol
from ._cutils import (
    _fill_matrix_linear,
    _fill_matrix_affine,
    _trace_one_linear,
    _trace_one_affine,
    _fill_matrix_linear_mn,
    _fill_matrix_affine_mn,
)


if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Iterable
    from skbio.sequence import SubstitutionMatrix
    from ._utils import SequenceLike


class MultiAlignResult(NamedTuple):
    path: AlignPath
    tree: TreeNode | None = None
    distmat: DistanceMatrix | None = None


def multi_align(
    sequences: Iterable[SequenceLike],
    /,
    sub_score: tuple[float, float] | SubstitutionMatrix | str = (1.0, -1.0),
    gap_cost: float | tuple[float, float] = 2.0,
    free_ends: bool = True,
    guide_tree: TreeNode | None = None,
    ids: Iterable[str] | None = None,
    atol: float = 1e-5,
    keep_tree: bool = False,
    keep_distmat: bool = False,
) -> MultiAlignResult:
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
    atol : float, optional
        Nonnegative finite absolute tolerance for traceback score comparisons,
        following :func:`pair_align`. Default is 1e-5. Set to zero for exact
        comparisons. Positive tolerance can select a slightly lower-scoring path;
        the DP maximum itself is unaffected. Also used for guide pairwise alignments.
    keep_tree : bool, optional
        If True, include the guide tree in the returned object. Default is False.
    keep_distmat : bool, optional
        If True, and if the guide tree is not provided, include the constructed
        distance matrix in the returned object. Default is False.

    Returns
    -------
    path : AlignPath
        One alignment path in original input sequence order, starting at position
        zero and consuming every sequence in full.
    tree : TreeNode or None
        Constructed or provided guide tree determining the merging order of sequences
        (if ``keep_tree`` is True).
    distmat : DistanceMatrix or None
        Distance matrix constructed based on pairwise alignments and used to compute
        the guide tree (if ``keep_distmat`` is True).

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
    >>> from skbio.alignment import multi_align, align_score

    Align three DNA sequences using default parameters.

    >>> from skbio.sequence import DNA
    >>> seqs = [DNA('CATTAACGT'),
    ...         DNA('CGTTACGGT'),
    ...         DNA('AGTTAACGG')]
    >>> path = multi_align(seqs).path
    >>> path
    <AlignPath, sequences: 3, positions: 11, segments: 7>

    Print the aligned sequences.

    >>> for seq in path.to_aligned(seqs):
    ...     print(seq)
    CA-TTAACGT-
    -CGTTA-CGGT
    -AGTTAACGG-

    The quality of the alignment can be evaluated using the `align_score` function,
    which calculates the sum-of-pairs (SP) score. It has the same default parameter
    settings as `multi_align` does.

    >>> from skbio.alignment import align_score
    >>> align_score((path, seqs))
    7.0

    Under the hood, the function performs pairwise alignments, calculates a distance
    matrix, then infers a guide tree which determines the merging order. The tree and
    distance matrix can be retained for diagnostic and educational purposes.

    >>> path, tree, dm = multi_align(seqs, keep_tree=True, keep_distmat=True)
    >>> print(tree.ascii_art())
              /-1
    ---------|
             |          /-0
              \--------|
                        \-2

    >>> print(dm)
    3x3 distance matrix
    IDs:
    '0', '1', '2'
    Data:
    [[ 0.          0.70444674  0.40215932]
     [ 0.70444674  0.          0.40215932]
     [ 0.40215932  0.40215932  0.        ]]

    One can supply a custom guide tree to skip the costly automatic pairwise alignment
    and tree building process. An accurate tree may improve the alignment quality.

    >>> from skbio.tree import TreeNode
    >>> tree = TreeNode.read(['((1,2),0);'])
    >>> path = multi_align(seqs, guide_tree=tree).path
    >>> for seq in path.to_aligned(seqs):
    ...     print(seq)
    CATTAACGT-
    CGTTA-CGGT
    AGTTAACGG-

    >>> align_score((path, seqs))
    9.0

    By default, sequences match taxa (tip names) of the tree by incremental indices
    '0', '1', '2'... Alternatively, explicit sequence IDs can be defined using the
    ``'id'`` key of sequence metadata or supplied by the ``ids`` parameter of this
    function.

    >>> for seq, id_ in zip(seqs, 'abc'):
    ...     seq.metadata['id'] = id_
    >>> tree = TreeNode.read(['((b,c),a);'])
    >>> res = multi_align(seqs, guide_tree=tree)

    One can customize the alignment parameters, including substitution scores, gap
    penalties, and terminal gap policy. Refer to :func:`pair_align` for details of
    the parameters.

    >>> params = dict(sub_score=(2, -3), gap_cost=(2, 5), free_ends=False)
    >>> path = multi_align(seqs, **params).path
    >>> for seq in path.to_aligned(seqs):
    ...     print(seq)
    CATTAACGT
    CGTTACGGT
    AGTTAACGG

    Supply the same parameters when calculating the alignment score.

    >>> align_score((path, seqs), **params)
    4.0

    The entire process of reading a multi-FASTA file of original sequences, performing
    multiple sequence alignment, and writing the aligned sequences into a multi-FASTA
    file is:

    >>> from skbio.io import read as sk_read  # doctest: +SKIP
    >>> from skbio.alignment import TabularMSA  # doctest: +SKIP
    >>> it = sk_read('input.fa', format='fasta', constructor=DNA)  # doctest: +SKIP
    >>> seqs = list(it)  # doctest: +SKIP
    >>> path = multi_align(seqs, **params).path  # doctest: +SKIP
    >>> msa = TabularMSA.from_path_seqs(path, seqs)  # doctest: +SKIP
    >>> msa.write('output.fa')  # doctest: +SKIP

    """
    seqs = list(sequences)
    if (n := len(seqs)) < 2:
        raise ValueError("At least two sequences are required.")
    for seq in seqs:
        if isinstance(seq, bytes):
            raise TypeError("Byte strings must be decoded before alignment.")
        if isinstance(seq, GrammaredSequence):
            gapped = seq.has_gaps()
        elif isinstance(seq, (Sequence, str)):
            chars = str(seq)
            gapped = "-" in chars or "." in chars
        else:
            gapped = False
        if gapped:
            raise ValueError("Input sequences must be ungapped.")
    ids = _get_seqids(seqs, ids)

    # Multiple alignment doesn't support custom terminal penalties.
    if not isinstance(free_ends, (bool, np.bool_)):
        raise TypeError("`free_ends` must be Boolean.")

    seqs, submat, _ = encode_sequences(seqs, sub_score)
    dtype = submat.dtype.type

    # Check symmetry. It's not that MSA absolutely require it. But its behavior will be
    # unpredictable with asymmetric substitution scores.
    if not np.isfinite(submat).all() or not np.array_equal(submat, submat.T):
        raise ValueError("Substitution scores must be finite and symmetric.")

    gap_o, gap_e = prep_gapcost(gap_cost, dtype=dtype)
    if not np.isfinite([gap_o, gap_e]).all() or min(gap_o, gap_e) < 0:
        raise ValueError("Gap costs must be finite and non-negative.")
    atol = _check_atol(atol, dtype=dtype)

    # Shrink substitution matrix to observed characters only. This accelerates the
    # calculation without changing the result. For example, if DNA sequences contain
    # only ACGT, the matrix will be (4, 4).
    alphabet, inv = np.unique(np.concatenate(seqs), return_inverse=True)
    seqs = np.split(inv, np.cumsum([len(x) for x in seqs])[:-1])
    submat = np.ascontiguousarray(submat[np.ix_(alphabet, alphabet)])

    # Re-usable array buffers that are capable of auto-growth
    works = ArrayWorkspace()

    # Decide merging order. If a guide tree is provided, convert it into a linkage
    # matrix (only the first two columns are needed). Otherwise, perform pairwise
    # alignments, compute a distance matrix using the Feng-Doolittle score distance
    # metric, then compute a guide tree using UPGMA and retain the linkage matrix.
    # `merges` is an index array of (n_seqs - 1, 2)
    if guide_tree is None:
        from scipy.cluster.hierarchy import linkage

        dm = _score_dists(seqs, submat, gap_o, gap_e, free_ends, works, atol)
        lm = linkage(dm, method="average")
        merges = lm[:, :2].astype(np.intp)
    else:
        merges = _tree_to_lnkmat(guide_tree, ids)

    eye = np.eye(len(submat), dtype=dtype)

    # Perform iterative profile alignment to merge all sequences. A profile consists of
    # a count matrix, a gap mask, and a list of member sequences.
    profiles = {}
    for parent, children in enumerate(merges, start=n):
        for child in children:
            if child < n:  # is a tip in linkage matrix
                seq = seqs[child]
                profiles[child] = (
                    eye[seq],
                    np.zeros((1, len(seq)), dtype=bool),
                    [child],
                )
        a, b = (profiles.pop(child) for child in children)
        profiles[parent] = _merge_align(
            a, b, submat, gap_o, gap_e, free_ends, works, atol
        )
    _, bits, order = profiles[2 * n - 2]

    # Reorder sequences to match input order before outputting
    path = AlignPath.from_bits(bits[np.argsort(order)])

    # Prepare extra outputs
    if not keep_tree:
        tree = None
    elif guide_tree is None:
        tree = TreeNode.from_linkage_matrix(lm, ids)
    else:
        tree = guide_tree
    if keep_distmat and guide_tree is None:
        dm = DistanceMatrix(dm, ids)
    else:
        dm = None

    return MultiAlignResult(path, tree, dm)


def multi_align_nucl(
    sequences: Iterable[SequenceLike],
    /,
    **kwargs: Any,
) -> MultiAlignResult:
    r"""Align multiple nucleotide sequences.

    This is a convenience wrapper of ``multi_align`` for nucleotide sequence alignment.
    It is preloaded with a scoring scheme consistent with BLASTN's defaults [1]_: match
    score = 2, mismatch score = -3, gap opening penalty = 5, gap extension penalty = 2.
    All parameters remain customizable. Refer to :func:`multi_align` for full
    documentation.

    See Also
    --------
    multi_align
    multi_align_prot

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK279684/

    Examples
    --------
    >>> from skbio.sequence import DNA
    >>> from skbio.alignment import multi_align_nucl
    >>> seqs = [DNA("CAGCTATATATCGCTACG"),
    ...         DNA("CTGCTTATATCCCTAGG"),
    ...         DNA("AAGCTATACATCCTTCACG")]
    >>> path = multi_align_nucl(seqs).path
    >>> for seq in path.to_aligned(seqs):
    ...     print(seq)
    CAGCTATATATCGCT-ACG
    CTGCT-TATATCCCT-AGG
    AAGCTATACATCCTTCACG

    """
    params: dict[str, Any] = dict(sub_score=(2.0, -3.0), gap_cost=(5.0, 2.0))
    params.update(kwargs)
    return multi_align(sequences, **params)


def multi_align_prot(
    sequences: Iterable[SequenceLike],
    /,
    **kwargs: Any,
) -> MultiAlignResult:
    r"""Align multiple protein sequences.

    This is a convenience wrapper of ``multi_align`` for protein sequence alignment.
    It is preloaded with a scoring scheme consistent with BLASTP's defaults [1]_:
    substitution matrix = BLOSUM62, gap opening penalty = 11, gap extension penalty
    = 1. All parameters remain customizable. Refer to :func:`multi_align` for full
    documentation.

    See Also
    --------
    multi_align
    multi_align_nucl

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK279684/

    Examples
    --------
    >>> from skbio.sequence import Protein
    >>> from skbio.alignment import multi_align_prot
    >>> seqs = [Protein("MKTAVLGHDPQRSIF"),
    ...         Protein("MKTSVLGHDPKRAIF"),
    ...         Protein("MRAAAVLNYDPPQSVF"),
    ...         Protein("MKTGAVLGHEDPQRTIF"),
    ...         Protein("MSTGVLGYDPQRSIL")]
    >>> path = multi_align_prot(seqs).path
    >>> for seq in path.to_aligned(seqs):
    ...     print(seq)
    MKTA-VLGH-DPQRSIF
    MKTS-VLGH-DPKRAIF
    MRAAAVLNY-DPPQSVF
    MKTGAVLGHEDPQRTIF
    MSTG-VLGY-DPQRSIL

    """
    params: dict[str, Any] = dict(sub_score="BLOSUM62", gap_cost=(11.0, 1.0))
    params.update(kwargs)
    return multi_align(sequences, **params)


def _score_dists(encoded, submat, gap_o, gap_e, free, works, atol):
    """Calculate alignment score distances between all sequences.

    Returns
    -------
    ndarray of shape (n * (n - 1) // 2,)
        Condensed (upper triangle) distance matrix between sequences.

    """
    n = len(encoded)

    # Force data type to be float64 (even if upstream code uses float32) in order to
    # retain precision in the Feng-Doolittle metric calculation.
    # TODO: Revisit this decision.
    dtype = np.float64
    gap_o64, gap_e64 = dtype(gap_o), dtype(gap_e)

    # Count occurrence per character per sequence
    counts = np.array(
        [np.bincount(x, minlength=len(submat)) for x in encoded], dtype=dtype
    )  # (n_sequences, n_symbols)

    # Intermediate to facilitate total substitution score calculation
    weight = counts @ submat.astype(dtype)  # (n_sequences, n_symbols)

    # Compute self-alignment scores
    params = (gap_o, gap_e, free, works, atol)
    s_self = [_align_pair(submat[seq], seq, *params, trace=False)[1] for seq in encoded]

    # Fill a condensed distance matrix
    dm = np.empty(n * (n - 1) // 2, dtype=dtype)
    pos = 0
    for i in range(n - 1):
        query = submat[encoded[i]]
        for j in range(i + 1, n):
            s_max = (s_self[i] + s_self[j]) / 2

            # Compute pairwise alignment
            moves, score = _align_pair(query, encoded[j], *params)
            path = _encode_path(moves)

            # Calculate total substitution score
            # S_subs = counts1.T @ submat @ counts2
            lens = path.lengths
            s_subs = weight[i] @ counts[j] / lens.sum()

            # Calculate total gap penalty
            gaps = path.states.ravel().astype(bool)
            if free:
                gaps[0] = gaps[-1] = False
            s_gaps = gap_o64 * gaps.sum() + gap_e64 * lens[gaps].sum()

            # Calculate score distance
            dm[pos] = _score_dist(score, s_subs, s_gaps, s_max)
            pos += 1

    return dm


def _score_dist(score, s_subs, s_gaps, s_max):
    """Calculate alignment score distance between a two sequences.

    The original metric (Eqs. 1-3 of Feng & Doolittle (1996)) is defined as:

        D = -ln S_eff, where S_eff = ((S - S_rand) / (S_max - S_rand))

    However, this equation may be undefined in certain edge cases. For example,
    aligning a homopolymer (e.g., "AAA") to itself will have both numerator and
    denominator = 0. Unrealistic substitution and gap score settings can also
    result in a non-positive numerator/denominator.

    Therefore, two modifications were introduced to this metric:

    1. When denominator = 0, directly return 0. The homopolymer case falls into this
       scenario. A distance of 0 between two identical sequences is justified.

    2. Clip S_eff to [eps, 1]. S_eff = 1 suggests that S = S_max, implicating maximum
       similarity between the two sequences. Thus D = 0 is justified. eps is a small,
       fixed floor. Here, we set it to 1e-6. Thus D_max = -ln 1e-6 ~= 13.82.

    """
    s_rand = s_subs - s_gaps
    numer = score - s_rand
    denom = s_max - s_rand
    if denom == 0:
        return 0.0
    s_eff = numer / denom
    s_eff = max(min(s_eff, 1.0), 1e-6)
    return -np.log(s_eff)


def _merge_align(aln1, aln2, submat, gap_o, gap_e, free, works, atol=0):
    """Merge two alignments.

    Parameters
    ----------
    aln1, aln2 : 3-tuple of (counts, gaps, order)
        The two alignments to merge. Elements are:

        counts : ndarray of shape (n_columns, n_symbols)
            Count per site per character in the alphabet. For one sequence, it is the
            one-hot encoding by alphabet. For multiple sequences, it is the sum of
            character frequencies in the current alignment.
        gaps : ndarray of shape (n_sequences, n_columns)
            Gap positions in the alignment.
        order: list of int
            Indices of merged sequences in merging order.

    Returns
    -------
    3-tuple of (counts, gaps, order)
        Output alignment. See above.

    Notes
    -----
    Output retains original row order.

    """
    counts1, gaps1, order1 = aln1
    counts2, gaps2, order2 = aln2
    dtype = submat.dtype

    n1, L1 = gaps1.shape
    n2, L2 = gaps2.shape

    # Pre-calculate an (n_columns_1, n_columns_2) score matrix before feeding into DP.
    # This is fully vectorized, but consumes extra memory. Workspace is pre-allocated
    # to reduce allocation overhead.
    # scores = (counts1 / n1) @ submat @ (counts2 / n2).T
    scores = works.get("scores", (L1, L2), dtype)
    weight = works.get("weight", (L1, submat.shape[0]), dtype)
    np.matmul(counts1, submat, out=weight)
    np.matmul(weight, counts2.T, out=scores)
    scores /= n1 * n2

    # Perform pairwise alignment using DP and return a dense path.
    path, _ = _align_pair(scores, None, gap_o, gap_e, free, works, atol)
    L = path.size

    # (2, n_columns) array of pre-merging column indices in the merged alignment. Gaps
    # are -1. This format is consistent with `AlignPath.to_indices`.
    # Here n_columns is the number of columns in the merged alignment.
    mask = path != np.array([1, 2])[:, None]
    indices = np.cumsum(mask, axis=1) - 1
    indices[~mask] = -1

    # Sum counts and concatenate gaps of the two alignments.
    counts = np.zeros((L, submat.shape[0]), dtype=dtype)
    gaps = np.ones((n1 + n2, L), dtype=bool)
    i = 0
    for counts_, bits_, n_, take in (
        (counts1, gaps1, n1, indices[0]),
        (counts2, gaps2, n2, indices[1]),
    ):
        mask = take >= 0
        cols = take[mask]
        j = i + n_
        counts[mask] += counts_[cols]
        gaps[i:j, mask] = bits_[:, cols]
        i = j

    return counts, gaps, order1 + order2


def _align_pair(query, target, gap_o, gap_e, free, works, atol, trace=True):
    """Align an encoded query/target with prepared costs in the query dtype.

    Returns
    -------
    path : ndarray of uint8 of (n_columns,)
        Dense alignment path (without run-length encoding).
    score : float
        Optimal alignment score.

    """
    if whole := target is None:
        m, n = query.shape
        scores = (query,)
    else:
        m, n = len(query), len(target)
        scores = (query, target)

    if affine := gap_o != 0:
        gaps = (gap_o, gap_e)
        fill_f = _fill_matrix_affine_mn if whole else _fill_matrix_affine
        trace_f = _trace_one_affine
    else:
        gaps = (gap_e,)
        fill_f = _fill_matrix_linear_mn if whole else _fill_matrix_linear
        trace_f = _trace_one_linear

    matrices = tuple(
        works.get(f"dp{k}", (m + 1, n + 1), query.dtype)
        for k in range(3 if affine else 1)
    )
    _init_matrices(matrices, gap_o, gap_e, False, free, free)
    fill_f(*matrices, *scores, *gaps, False)

    score, stops = _one_stop(matrices[0], False, free, free)
    if not trace:
        return None, float(score)

    i, j = stops[0]
    path = works.get("path", (m + n,), np.uint8)
    pos, _, _ = _trailing_gaps(path, m + n, i, j, m, n, True, True)
    pos, i, j = trace_f(path, pos, i, j, *matrices, gap_e, False, atol)
    pos, _, _ = _leading_gaps(path, pos, i, j, True, True)

    return path[pos:], float(score)
