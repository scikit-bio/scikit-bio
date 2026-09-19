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
from ._utils import encode_sequences, prep_gapcost, _get_seqids, _prep_atol
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
    r"""Perform progressive alignment of multiple sequences.

    .. versionadded:: 0.7.4

    A versatile workflow for multiple sequence alignment (MSA) by integrating the
    pairwise alignment engine (:func:`pair_align`), alignment distance calculation,
    and UPGMA tree construction following the progressive alignment approach.

    Parameters
    ----------
    sequences : iterable of Sequence, str, or sequence of scalar
        Sequences to be aligned. Must be non-empty and ungapped. At least two
        sequences must be provided.
    sub_score : tuple of (float, float), SubstitutionMatrix, or str, optional
        Score of a substitution. May be two numbers (match, mismatch), a substitution
        matrix, or its name. See :func:`pair_align` for details. Default is
        (1.0, -1.0).
    gap_cost : float or tuple of (float, float), optional
        Penalty of a gap. May be one (linear) or two numbers (affine). See
        :func:`pair_align` for details. Default is 2.0.
    free_ends : bool, optional
        If True (default), gaps at the sequence terminals are free from penalization.
    guide_tree : TreeNode, optional
        A guide tree determining the merging order of sequences. Must be strictly
        bifurcating. Tip names must match sequence IDs. If not provided, the function
        will align all sequence pairs, calculate score-based distances, and compute a
        guide tree using UPGMA. A provided guide tree will skip these procedures.
    ids : iterable of str, optional
        Unique identifiers in input order to match tips in the provided guide tree.
        tree. If not provided, the function will use sequence metadata ``'id'`` if
        present in every sequence and unique, or use ``['0', '1', ...]`` if none is
        present.
    atol : float, optional
        Absolute tolerance in comparing scores of alternative alignment paths. See
        :func:`pair_align` for details. Default is 1e-5.
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
        Distance matrix constructed based on pairwise alignment scores and used to
        compute the guide tree (if ``keep_distmat`` is True).

    Raises
    ------
    ValueError
        If inputs are empty or gapped, scores or costs are invalid, identifiers do
        not match, or a guide tree is not binary with exactly the required tips.
    TypeError
        If ``free_ends`` is not Boolean or ``guide_tree`` is not a TreeNode.

    See Also
    --------
    pair_align
    align_score
    AlignPath

    Notes
    -----
    This function implements the classic progressive alignment method for multiple
    sequence alignment, originally introduced in [1]_, with later improvements
    described in [2]_ and [3]_. Compared with the historical method, this
    implementation represents a refined form of progressive alignment that is often
    described in educational materials. Specifically, the algorithm consists of the
    following steps, detailed in reverse order:

    **Procedures**

    An alignment of all sequences is constructed by iteratively merging sub-alignments
    containing one or more sequences. This function adopts the **profile alignment**
    approach [3]_, which aligns two sub-alignments ("profiles") using the same dynamic
    programming (DP) algorithm for pairwise sequence alignment (:func:`pair_align`).
    Refer to the latter's documentation for settings and considerations.

    The substitution score :math:`S` for two matching columns between profiles is
    calculated as the average substitution score :math:`s` across all pairs of
    characters from the two columns:

    .. math::
        S = \frac{1}{mn}\sum_{x\in A}\sum_{y\in B}s(x,y)

    where :math:`x` and :math:`y` are characters in the two columns :math:`A` and
    :math:`B`, which contain :math:`m` and :math:`n` rows (sequences), respectively.

    Under the "once a gap, always a gap" rule [1]_, existing gaps within each profile
    are treated as neutral characters and assigned a substitution score of 0 with any
    character. Gap penalties are calculated only for gaps introduced during the DP
    alignment.

    The order of merging is determined by a **guide tree**. The program traverses the
    tree in postorder and merges the two child sub-alignments at each internal node. If
    the guide tree is not explicitly supplied, the program computes one using the UPGMA
    method (see :func:`~skbio.tree.upgma`) from a distance matrix containing all
    pairwise sequence distances. Each distance is calculated by aligning the two
    sequences :math:`a` and :math:`b` using DP (see :func:`pair_align`) and
    normalizing the alignment score :math:`S_{a,b}` as follows [1]_:

    .. math::
        D = -\ln S_{\mathrm{eff}} = -\ln\left(\frac{S_{a,b} - S_{\mathrm{rand}}}
        {S_{\mathrm{iden}} - S_{\mathrm{rand}}}\right)

    where :math:`S_{\mathrm{iden}} = (S_{a,a} + S_{b,b})/2` is the average score of the
    two sequences aligned to themselves. :math:`S_{\mathrm{rand}}` is the *random
    score*, calculated as [2]_:

    .. math::
        S_{\mathrm{rand}} = \frac{1}{L}\sum_{x\in a}\sum_{y\in b}s(x,y)N_a(x)N_b(y) - G

    where :math:`L` is the length of the alignment, :math:`G` is the total gap penalty,
    and :math:`N` is the number of occurrences of a character in the corresponding
    source sequence.

    There are two notes regarding this calculation. First, a factor of 100 is omitted
    from the effective score :math:`S_{\mathrm{eff}}` compared with the original work.
    Second, to guard against edge cases that would produce undefined, infinite, or
    negative distances (e.g., when aligning two homopolymers), :math:`S_{\mathrm{eff}}`
    is clipped to the range [1e-6, 1] in this implementation.

    **Solution quality**

    The **sum-of-pairs** (SP) score is the optimality criterion for multiple sequence
    alignment. This metric can be calculated by applying the :func:`align_score`
    function to the resulting alignment. It should be noted that progressive alignment
    is a heuristic algorithm and the resulting alignment is not guaranteed to be
    optimal.

    **Computational efficiency**

    The algorithm is dominated by the all-vs-all pairwise alignment step, which takes
    *O*\(*n*:sup:`2` *L*:sup:`2`) time for *n* sequences of comparable length *L*. Peak
    memory usage is *O*\(*n*:sup:`2` + *L*:sup:`2`) for storing the distance matrix and
    each DP matrix (the algorithm reuses memory for DP matrices). When a guide tree is
    supplied, time reduces to *O*\(*nL*:sup:`2` + *n*:sup:`2` *L*), and memory to
    *O*\(*L*:sup:`2` + *nL*).

    **Terminal gap policy**

    The function defaults to ``free_ends=True`` which prevents terminal gaps from
    being penalized. This setting is broadly applicable to homologous sequences with
    incomplete coverage, different domain boundaries, or terminal extensions. However,
    it can favor short overlaps between weakly related sequences. When sequences are
    expected to span the same homologous region with defined boundaries (e.g., the
    coding sequence of a gene), setting ``free_ends=False`` is often preferable.

    References
    ----------
    .. [1] Feng, D. F., & Doolittle, R. F. (1987). Progressive sequence alignment as a
       prerequisite to correct phylogenetic trees. Journal of Molecular Evolution,
       25(4), 351-360.
    .. [2] Feng, D. F., & Doolittle, R. F. (1996). [21] Progressive alignment of amino
       acid sequences and construction of phylogenetic trees from them. In Methods in
       enzymology (Vol. 266, pp. 368-382). Academic Press.
    .. [3] Corpet, F. (1988). Multiple sequence alignment with hierarchical clustering.
       Nucleic Acids Research, 16(22), 10881-10890.

    Examples
    --------
    >>> from skbio.alignment import multi_align

    Align three DNA sequences using default parameters and obtain an alignment path.

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

    Or convert the alignment path and sequences into a :class:`TabularMSA` object.

    >>> from skbio.alignment import TabularMSA
    >>> msa = TabularMSA.from_path_seqs(path, seqs)
    >>> msa
    TabularMSA[DNA]
    ----------------------
    Stats:
        sequence count: 3
        position count: 11
    ----------------------
    CA-TTAACGT-
    -CGTTA-CGGT
    -AGTTAACGG-

    The quality of the alignment can be evaluated using the :func:`align_score`
    function, which calculates the sum-of-pairs (SP) score, and has the same default
    parameter settings as ``multi_align`` does.

    >>> from skbio.alignment import align_score
    >>> align_score((path, seqs))
    7.0

    Under the hood, this function performs pairwise alignments, calculates a distance
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
    >>> seqs = list(
    ...     sk_read('input.fa', format='fasta', constructor=DNA)
    ... )  # doctest: +SKIP
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
    atol = _prep_atol(atol, dtype=dtype)

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
        if not isinstance(guide_tree, TreeNode):
            raise TypeError("`guide_tree` must be a TreeNode.")
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
