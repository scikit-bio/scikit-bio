# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import NamedTuple, Optional, Union, Tuple, List, TYPE_CHECKING

import numpy as np

from skbio.alignment import PairAlignPath
from ._utils import encode_sequences, prep_gapcost
from ._cutils import (
    _fill_linear_matrix,
    _fill_affine_matrices,
    _trace_one_linear,
    _trace_one_affine,
)

if TYPE_CHECKING:  # pragma: no cover
    from skbio.sequence import SubstitutionMatrix
    from ._utils import SequenceLike


class PairAlignResult(NamedTuple):
    score: float
    paths: Optional[List["PairAlignPath"]] = None
    matrices: Optional[Tuple[np.ndarray, ...]] = None


def pair_align(
    seq1: "SequenceLike",
    seq2: "SequenceLike",
    /,
    mode: str = "global",
    sub_score: Union[Tuple[float, float], "SubstitutionMatrix", str] = (1.0, -1.0),
    gap_cost: Union[float, Tuple[float, float]] = 2.0,
    free_ends: Union[bool, Tuple[bool, bool], Tuple[bool, bool, bool, bool]] = True,
    trim_ends: bool = False,
    max_paths: Optional[int] = 1,
    atol: float = 1e-5,
    keep_matrices: bool = False,
) -> PairAlignResult:
    r"""Perform pairwise alignment of two sequences.

    A versatile, efficient and generalizable function for pairwise alignment using
    dynamic programming. It supports:

    * Global, local and semi-global alignments (with all four ends customizable).
    * Nucleotide, protein, and un-grammared sequences, plain strings (ASCII and
      Unicode), words/tokens, and numbers.
    * Match/mismatch scores or substitution matrix.
    * Linear and affine gap penalties.
    * Integer, decimal and infinite scores.
    * Returning one, multiple or all optimal alignment paths.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    seq1 : :class:`~skbio.sequence.Sequence`, str, or sequence of scalar
        The first sequence to be aligned.

    seq2 : :class:`~skbio.sequence.Sequence`, str, or sequence of scalar
        The second sequence to be aligned.

    mode : {'global', 'local'}, optional
        Mode of alignment. Options are:

        - "global" (default): Global alignment, as in the Needleman-Wunsch algorithm.
          The two sequences are aligned end-to-end. See also ``free_ends`` below.
        - "local": Local alignment, as in the Smith-Waterman algorithm. It identifies
          similar regions between two sequences.

    sub_score : tuple of (float, float), SubstitutionMatrix, or str, optional
        Score of a substitution. May be one of the following:

        - Tuple of two numbers: Match score (same symbol) and mismatch score (different
          symbols).
        - :class:`~skbio.sequence.SubstitutionMatrix`: A matrix of substitution scores
          between all symbols in the alphabet.
        - String: Name of the substitution matrix that can be recognized by
          ``SubstitutionMatrix.by_name``, such as "NUC.4.4" or "BLOSUM62".

        Default is (1.0, -1.0) (match, mismatch).

    gap_cost : float or tuple of (float, float), optional
        Penalty of a gap. Values are usually positive, representing subtractions from
        the alignment score. May be one of the following:

        - One number: Linear gap penalty. Each gap position is penalized by this value
          (*g*). A contiguous gap of length *k* has a total penalty of *g* * *k*.
        - Tuple of two numbers: Affine gap penalty. The two numbers (*o*, *e*)
          represent gap opening penalty and gap extension penalty, respectively. A
          contiguous gap of length *k* has a total penalty of *o* + *e* * *k*.
          See also notes below.

        Default is 2.0 (linear gap penalty).

        .. note::
           Infinites are valid values for ``sub_score`` and ``gap_cost``. For example,
           ``gap_cost=np.inf`` disables gaps. ``sub_score=(..., -np.inf)`` disables
           mismatches.

    free_ends : bool, 2-tuple of bool, or 4-tuple of bool, optional
        Whether gaps at the sequence terminals are free from penalization. Relevant
        when ``mode`` is "global". Alignment with free ends is known as semi-global
        (or "glocal") alignment. Values can be:

        - False: Penalize terminal gaps using the same method defined by ``gap_cost``.
          This behavior is the true global alignment.
        - True: Do not penalize any terminal gap. This behavior is known as "overlap",
          as it identifies the maximum overlap between two sequences.
        - Tuple of two Booleans: Whether terminal gaps of seq1 and seq2 are free,
          respectively. For example, (True, False) is the typical semi-global
          alignment, useful for aligning a short seq1 within a long seq2.
        - Tuple of four Booleans: Whether leading gap of seq1, trailing gap of seq1,
          leading gap of seq2, and trailing gap of seq2 are free, respectively. For
          example, (False, True, True, False) joins the tail of seq1 with the head
          of seq2.

        Default is True (overlap).

    trim_ends : bool, optional
        If True, penalty-free terminal gaps defined by ``free_ends`` will be excluded
        from returned paths. Useful for locating a short query in a long reference.
        Default is False.

    max_paths : int, optional
        Maximum number of alignment paths to return. Default is 1, which is generated
        through a performance-oriented traceback algorithm. A value larger than 1 will
        trigger a less efficient traceback algorithm to enumerate up to this number of
        paths. Setting it as None will return all paths. However, be cautious that the
        total number of paths may be extremely large. Setting it as 0 will disable
        traceback and return no path (None).

        .. note::
           When ``mode="local"``, it is possible that there is no path with a score >
           0. In this case, an empty list will be returned. ``mode="global"``
           guarantees to return at least one path, even if it's completely misaligned.

    atol : float, optional
        Absolute tolerance in comparing scores of alternative alignment paths. This is
        to ensure floating-point arithmetic safety when ``sub_score`` or ``gap_cost``
        involve decimal numbers. Default is 1e-5. Setting it to 0 or None will slightly
        increase performance, and is usually safe when there are only integers (e.g.,
        2.0) or half integers (e.g., 2.5) involved.

        .. note::
           Relative tolerance is not involved in the calculation.

    keep_matrices : bool, optional
        Whether to include the alignment matrix(ces) in the returned value. They are
        typically for diagnostic or educational purposes. Default is False, which lets
        the memory space free up after the function finishes.

    Returns
    -------
    score : float
        Optimal alignment score.

    paths : list of :class:`~skbio.alignment.PairAlignPath`, optional
        Alignment paths. Up to ``max_paths`` paths will be returned. Note that all
        paths are optimal and share the same alignment score.

    matrices : tuple of ndarray of shape (m + 1, n + 1), optional
        Alignment matrices generated during the computation (if ``keep_matrices`` is
        True). *m* and *n* are the lengths of seq1 and seq2, respectively.

        - For linear gap penalty, one main matrix will be returned.
        - For affine gap penalty, the main matrix, plus an insertion matrix (gap in
          seq1) and a deletion matrix (gap in seq2) will be returned.

    See Also
    --------
    align_score
    skbio.alignment.PairAlignPath

    Notes
    -----
    This function implements the classic dynamic programming (DP) method for pairwise
    sequence alignment. It is commonly known as the Needleman-Wunsch algorithm [1]_ for
    global alignment, or the Smith-Waterman algorithm [2]_ for local alignment. These
    two algorithms use linear gap penalty. When affine gap penalty is specified, the
    underlying method is the Gotoh algorithm [3]_, with later modifications [4]_. The
    algorithms are exact and the output alignments are guaranteed to be optimal.

    The algorithms are quadratic (*O*\(*mn*\)) in both time and space, where *m* and
    *n* are the lengths of the two sequences, respectively. Affine gap penalty costs
    3x as much memory and 2-3x as much runtime than linear gap penalty.

    **Parameters**

    The default parameters ``sub_score=(1, -1), gap_cost=2`` is a simple scoring scheme
    that quickly captures sequence similarity. However, in bioinformatics applications,
    one usually wants to replace them with more realistic parameters according to the
    specific task. For reference, below are the default parameters of some common
    bioinformatics programs:

    - NCBI MegaBLAST (nucleotide): ``sub_score=(1, -2), gap_cost=2.5``
    - NCBI BLASTN (nucleotide): ``sub_score=(2, -3), gap_cost=(5, 2)``
      (see also :func:`pair_align_nucl`)
    - NCBI BLASTP (protein): ``sub_score="BLOSUM62", gap_cost=(11, 1)``
      (see also :func:`pair_align_prot`)
    - EMBOSS Needle/Water (protein): ``sub_score="BLOSUM62", gap_cost=(9.5, 0.5)``
      (for nucleotide, replace ``sub_score`` with ``"NUC.4.4"``)
    - EBI FASTA (protein): ``sub_score="BLOSUM50", gap_cost=(8, 2)``

    .. note::
       These are scoring schemes only. ``pair_align`` does not reproduce the output
       of the programs mentioned.

    The flexibility of these parameters, especially the support for infinity, enables
    some common tasks in text analysis:

    - Edit (Levenshtein) distance: ``sub_score=(0, -1), gap_cost=1`` (negate the score
      to get the distance)
    - Longest common subsequence (LCS):
      ``mode="local", sub_score=(1, -np.inf), gap_cost=0``
    - Longest common substring:
      ``mode="local", sub_score=(1, -np.inf), gap_cost=np.inf``

    **Input sequences**

    The format of input sequences is flexible. The function is most efficient when they
    are subclasses of :class:`~skbio.sequence.GrammaredSequence`, such as ``DNA`` and
    ``Protein``, or any user-customized class. They have a finite alphabet, and permit
    auto-replacing degenerate codes with wildcard (such as nucleotides to "N" and amino
    acids to "X"), as long as the substitution matrix has the wildcard character.

    For non-grammared sequences, such as plain strings, bytes, and integer arrays, the
    the function will attempt to encode them into ASCII codes (*n* = 128), which permit
    efficient indexing. If ASCII encoding is not possible (for example lists of words,
    large or negative integers, or floats), unique symbols will be extracted from the
    sequences and indexed.

    **Performance**

    The underlying dynamic programming kernel is a plain dual loop structure, without
    further vectorization or parallelization techniques.

    This function does not discriminate between seq1 and seq2. Nevertheless, aligning a
    shorter seq1 (often referred to as "query") against a longer seq2 (often referred to
    as "target" or "reference") is usually more efficient than the other way around.

    The algorithm defaults to float32, which is sufficient for most use cases. It also
    supports float64, with a higher memory cost (2x) and moderately increased runtime.
    If float64 is what you need, you may supply a substitution matrix of float64 type,
    using e.g., ``SubstitutionMatrix.identity("ACGT", 1, -2, dtype="float64")``, which
    will be respected by the function without casting throughout calculation.

    **Affine gap penalty**

    Under the affine gap penalty mode, the penalty of a contiguous gap of length
    :math:`k` is calculated as:

    .. math::

       G(k) = o + e \times k \tag{1}

    where :math:`o` is the gap opening penalty and :math:`e` is the gap extension
    penalty.

    It should be noted that, discrepancy exists among literature and implementations
    regarding whether gap extension penalty should apply to the first position of a
    gap. scikit-bio's equation is consistent with multiple common alignment tools,
    such as BLAST [5]_, Minimap2, SeqAn3, and WFA2-lib.

    Meanwhile, multiple other tools, such as EMBOSS, parasail, Biopython and Biotite,
    use the following equation instead:

    .. math::

       G(k) = o + e \times (k - 1) \tag{2}

    Therefore, if you intend to reproduce the behavior of a software tool of the
    latter category using scikit-bio, you will need to subtract :math:`e` from
    :math:`o` when adopting its parameters. For example, EMBOSS' default parameters
    ``o=10, e=0.5`` will become ``o=9.5, e=0.5`` in scikit-bio. Vice versa.

    .. versionchanged:: 0.7.0
        Previous alignment algorithms in scikit-bio used Eq. 2. These functions were
        deprecated in 0.5.x and will be removed in 0.6.x. Future functions will
        uniformly use Eq. 1.

    Related: A simple but less common scenario is constant gap penalty, where a
    contiguous gap has a constant penalty :math:`g` regardless of its length. This can
    be specified as ``gap_cost=(g, 0)``.

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

    .. [5] https://www.ncbi.nlm.nih.gov/books/NBK62051/

    Examples
    --------
    >>> from skbio.sequence import DNA, Protein
    >>> from skbio.alignment import pair_align

    Align two DNA sequences using default parameters (global alignment with free end
    gaps, match/mismatch scores, linear gap penalty, returning just one path):

    >>> seq1 = DNA('ACTACCAGATTACTTACGGATCAGGTACTTGCCAACAA')
    >>> seq2 = DNA('CGAAACTACTAGATTACGGATCTTACTTTCCAGCAAGG')
    >>> res = pair_align(seq1, seq2)

    The result is a named tuple consisting of score, path(s), and matrices (off by
    default). The score represents the "goodness" of the alignment.

    >>> res.score
    12.0

    The alignment path can be represented by a
    :wiki:`CIGAR string <Sequence_alignment#CIGAR_Format>` . See
    :class:`PairAlignPath` for details.

    >>> path = res.paths[0]
    >>> path
    <PairAlignPath, positions: 44, segments: 7, CIGAR: '4I13M4D6M2D13M2I'>

    Extract aligned sequences:

    >>> aln = path.to_aligned((seq1, seq2))
    >>> aln
    ['----ACTACCAGATTACTTACGGATCAGGTACTTGCCAACAA--',
     'CGAAACTACTAGATTAC----GGATCT--TACTTTCCAGCAAGG']

    Alternatively, you can do ``TabularMSA.from_path_seqs(path, (seq1, seq2))`` to
    create a dedicated sequence alignment structure. See :class:`TabularMSA` for
    details.

    Align two protein sequences using local alignment, substitution matrix, and affine
    gap penalty (the default parameters of BLASTP):

    >>> seq1 = Protein('MKRTLKGHFVQWC')
    >>> seq2 = Protein('MQMLKTHYAQTRN')
    >>> score, paths, _ = pair_align(seq1, seq2, mode='local',
    ...                              sub_score='BLOSUM62', gap_cost=(11, 1))
    >>> score
    23.0

    >>> paths[0].to_aligned((seq1, seq2))
    ['LKGHFVQ', 'LKTHYAQ']

    Search a sequencing read against a reference genome using the default parameters of
    MegaBLAST and return locations of all hits.

    >>> query = "ACCGT"
    >>> target = "AAACGCTACCGTCCGTAGACCGTGACCGTGCGAAGC"
    >>> score, paths, _ = pair_align(query, target, mode='global',
    ...                              sub_score=(1, -2), gap_cost=2.5,
    ...                              free_ends=(True, False), trim_ends=True,
    ...                              max_paths=None)
    >>> len(paths)
    3

    >>> for x in paths:
    ...     print(x.ranges[1])
    [ 7 12]
    [18 23]
    [24 29]

    Note: ``ranges[1]`` stores the start and stop positions of the aligned region in
    the target (``ranges[0]`` stores those of the query). The numbers are BED-like
    (i.e., 0-based, half-open).

    Calculate the edit distance between two sentences (lists of words).

    >>> text1 = 'The quick brown fox jumps over the lazy dog'.split()
    >>> text2 = 'The slow brown wolf jumps over a lazy dog'.split()
    >>> res = pair_align(text1, text2, mode='global',
    ...                  sub_score=(0, -1), gap_cost=1,
    ...                  free_ends=False, max_paths=0)
    >>> -int(res.score)
    3

    """
    # This function implements the classic dynamic programming method for pairwise
    # sequence alignment. While the most time-consuming step (matrix filling) is done
    # in Cython, this function retains as many steps as possible in Python, thereby
    # "outsourcing" computation and optimization to the upstream library (NumPy).

    # Determine alignment mode (False - global, True - local)
    local = _prep_mode(mode)

    # Determine end gap policy. These are ternary flags (0 - penalize, 1 - free and
    # keep, 2 - free and trim).
    lead1, trail1, lead2, trail2 = _prep_free_ends(local, free_ends, trim_ends)

    # Prepare sequences and substitution matrix.
    # If `sub_score` consists of match/mismatch scores, an identity matrix will be
    # created or retrieved from cache. The two sequences are converted into indices in
    # the substitution matrix, which facilitate lookup and memory locality.
    (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), sub_score)

    # Profile seq1 (query), which will be aligned against seq2 (target).
    # The profile is essentially a position-specific scoring matrix (PSSM). This design
    # is useful for future implementation of profile search (like PSI-BLAST) and one-
    # vs-many search.
    # Query and target must be C-contiguous for efficient loops. submat, seq1 and seq2
    # are guaranteed to be C-contiguous by upstream code. submat[seq1] should be C-
    # contiguous too due to NumPy's advanced indexing. But `ascontiguousarray` is still
    # called here to be safe.
    query, target = np.ascontiguousarray(submat[seq1]), seq2
    dtype = query.dtype.type

    # Prepare affine or linear gap penalties.
    gap_open, gap_extend = prep_gapcost(gap_cost, dtype=dtype)
    affine = gap_open != 0

    # Allocate alignment matrices.
    # There is one matrix for linear gap penalty or three matrices for affine gap
    # penalty. Each matrix is (m + 1) by (n + 1). They can become quite large and
    # challenge the memory capacity. To overcome this, future implementation of
    # Hirschberg's algorithm with linear space is desired.
    matrices = _alloc_matrices(query.shape[0], target.size, affine, dtype=dtype)

    # Initialize alignment matrices.
    _init_matrices(matrices, gap_open, gap_extend, local, lead1, lead2)

    # Fill alignment matrices (quadratic; compute-intensive).
    if affine:
        _fill_affine_matrices(*matrices, query, target, gap_open, gap_extend, local)
    else:
        _fill_linear_matrix(matrices[0], query, target, gap_extend, local)

    # Get optimal alignment score and corresponding stop(s).
    if max_paths == 1 or max_paths == 0:
        score, stops = _one_stop(matrices[0], local, trail1, trail2)
    else:
        score, stops = _all_stops(matrices[0], local, trail1, trail2)

    # Traceback from each stop to reconstruct optimal alignment path(s).
    args = (gap_open, gap_extend, local, lead1, lead2, trail1, trail2, atol)
    if max_paths == 0:  # traceback is disabled
        paths = None
    elif local and abs(score) <= atol:  # no path is found
        paths = []
    elif max_paths == 1:
        i, j = stops[0]
        paths = [_traceback_one(i, j, matrices, *args)]
    else:
        paths = _traceback_all(stops, max_paths, matrices, query, target, *args)

    # Discard or keep matrices.
    if not keep_matrices:
        matrices = None
    elif affine:
        _fill_nan(matrices)

    return PairAlignResult(float(score), paths, matrices)


def pair_align_nucl(
    seq1: "SequenceLike",
    seq2: "SequenceLike",
    /,
    **kwargs,
) -> PairAlignResult:
    r"""Align two nucleotide sequences.

    This is a convenience wrapper of ``pair_align`` for nucleotide sequence alignment.
    It is preloaded with a scoring scheme consistent with BLASTN's defaults [1]_: match
    score = 2, mismatch score = -3, gap opening penalty = 5, gap extension penalty = 2.
    All parameters remain customizable. Refer to :func:`pair_align` for full
    documentation.

    See Also
    --------
    pair_align
    pair_align_prot

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK279684/

    Examples
    --------
    >>> from skbio.sequence import DNA
    >>> from skbio.alignment import pair_align_nucl
    >>> seq1 = DNA('GATCGTC')
    >>> seq2 = DNA('ATCGCTC')
    >>> res = pair_align_nucl(seq1, seq2)
    >>> res.score
    5.0

    >>> res.paths[0]
    <PairAlignPath, positions: 8, segments: 4, CIGAR: '1D4M1I2M'>

    >>> res.paths[0].to_aligned((seq1, seq2))
    ['GATCG-TC', '-ATCGCTC']

    """
    params = dict(sub_score=(2.0, -3.0), gap_cost=(5.0, 2.0))
    params.update(kwargs)
    return pair_align(seq1, seq2, **params)


def pair_align_prot(
    seq1: "SequenceLike",
    seq2: "SequenceLike",
    /,
    **kwargs,
) -> PairAlignResult:
    r"""Align two protein sequences.

    This is a convenience wrapper of ``pair_align`` for protein sequence alignment.
    It is preloaded with a scoring scheme consistent with BLASTP's defaults [1]_:
    substitution matrix = BLOSUM62, gap opening penalty = 11, gap extension penalty
    = 1. All parameters remain customizable. Refer to :func:`pair_align` for full
    documentation.

    See Also
    --------
    pair_align
    pair_align_nucl

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK279684/

    Examples
    --------
    >>> from skbio.sequence import Protein
    >>> from skbio.alignment import pair_align_prot
    >>> seq1 = Protein('PKKKRKV')
    >>> seq2 = Protein('PAAKRVKLD')
    >>> res = pair_align_prot(seq1, seq2)
    >>> res.score
    11.0

    >>> res.paths[0]
    <PairAlignPath, positions: 9, segments: 2, CIGAR: '7M2I'>

    >>> res.paths[0].to_aligned((seq1, seq2))
    ['PKKKRKV--', 'PAAKRVKLD']

    """
    params = dict(sub_score="BLOSUM62", gap_cost=(11.0, 1.0))
    params.update(kwargs)
    return pair_align(seq1, seq2, **params)


def _prep_mode(mode):
    """Prepare pairwise alignment mode.

    Parameters
    ----------
    mode : str
        Alignment mode.

    Returns
    -------
    bool
        Whether alignment mode is local.

    """
    if mode == "global":
        return False
    elif mode == "local":
        return True
    else:
        raise ValueError(f"Invalid mode: {mode}.")


def _prep_free_ends(local, free_ends, trim_ends):
    """Prepare end gap policy for pairwise alignment.

    Parameters
    ----------
    local : bool
        Local or global alignment.
    free_ends : bool, (bool, bool), or (bool, bool, bool, bool)
        Free ends policy.
    trim_ends : bool
        End trimming policy.

    Returns
    -------
    4-tuple of {0, 1, 2}
        Ternary flags. Positions:
            0. Leading gaps of sequence 1.
            1. Trailing gaps of sequence 1.
            2. Leading gaps of sequence 2.
            3. Trailing gaps of sequence 2.
        Codes:
            0. Penalize.
            1. Free and keep.
            2. Free and trim.

    """
    # not relevant
    if local:
        return 2, 2, 2, 2
    trim = bool(trim_ends)

    # all four ends
    if np.isscalar(free_ends):
        free = bool(free_ends)
        free += free and trim
        return free, free, free, free

    # both ends for seq1 and seq2
    elif (n := len(free_ends)) == 2:
        free1, free2 = (x + (x and trim) for x in map(bool, free_ends))
        return free1, free1, free2, free2

    # leading and trailing ends for seq1 and seq2
    elif n == 4:
        return tuple(x + (x and trim) for x in map(bool, free_ends))
    else:
        raise ValueError("`free_ends` must be one, two or four Booleans.")


def _alloc_matrices(m, n, affine, dtype=np.float32):
    """Allocate alignment matrix(ces).

    Parameters
    ----------
    m, n: int
        Lengths of sequences 1 and 2, respectively.
    affine : bool
        Affine (True) or linear (False) gap penalty.
    dtype : type, optional
        Data type (np.float32 or np.float64).

    Returns
    -------
    tuple of ndarray of shape (m + 1, n + 1)
        Alignment matrix(ces).

    Notes
    -----
    The arrays should be C-contiguous to facilitate row-major operations. NumPy's
    default array is already C-contiguous. This is also enforced by the unit test.

    """
    shape = (m + 1, n + 1)
    scomat = np.empty(shape, dtype=dtype)
    if not affine:
        return (scomat,)
    insmat = np.empty(shape, dtype=dtype)
    delmat = np.empty(shape, dtype=dtype)
    return scomat, insmat, delmat


def _init_matrices(matrices, gap_open, gap_extend, local, lead1, lead2):
    """Initialize alignment matrix(ces) by populating first column and row.

    Parameters
    ----------
    matrices : tuple of ndarray of shape (m + 1, n + 1)
        Alignment matrices.
    gap_open : float
        Gap opening penalty.
    gap_extend : float
        Gap extension penalty.
    local : bool
        Local or global alignment.
    lead1, lead2 : bool
        Whether leading gaps of sequence 1 and 2 are free.

    Notes
    -----
    The initial values in the insertion and deletion matrices are -inf, which follows
    conventions and should be correct. However, one should note from Flouri et al.
    (2015) [1]_ that the original Gotoh algorithm (1982) made a mistake in it by using
    gap_open + k * gap_extend. Flouri et al. suggested that the minimum values must be
    2 * gap_open + k * gap_extend.

    References
    ----------
    .. [1] Flouri, T., Kobert, K., Rognes, T., & Stamatakis, A. (2015). Are all global
       alignment algorithms and implementations correct?. bioRxiv, 031500.

    """
    scomat = matrices[0]
    m1, n1 = scomat.shape

    # initialize main scoring matrix
    scomat[0, 0] = 0
    if local:
        scomat[1:m1, 0] = 0
        scomat[0, 1:n1] = 0
    else:
        series = np.arange(1, max(m1, n1), dtype=scomat.dtype)
        series *= -gap_extend
        if gap_open:
            series -= gap_open
        scomat[1:m1, 0] = 0 if lead2 else series[: m1 - 1]
        scomat[0, 1:n1] = 0 if lead1 else series[: n1 - 1]

    # initialize insertion and deletion matrices
    if gap_open:
        insmat = matrices[1]
        delmat = matrices[2]
        insmat[1:m1, 0] = -np.inf
        delmat[0, 1:n1] = -np.inf
        # Flouri's minimum value (see Notes):
        # series -= gap_open
        # insmat[1:m1, 0] = series[:m1 - 1]
        # delmat[0, 1:n1] = series[:n1 - 1]


def _fill_nan(matrices):
    """Fill empty cells of affine matrices with NaN before returning.

    These are not useful during alignment, but can make things less confusing for
    diagnostic or educational purposes.

    This function supports arbitrary number of affine pieces.

    """
    n = len(matrices)
    # odd numbers are insertion matrices
    for i in range(1, n, 2):
        matrices[i][0, :] = np.nan
    # even numbers are deletion matrices
    for i in range(2, n, 2):
        matrices[i][:, 0] = np.nan


def _one_stop(scomat, local, trail1, trail2):
    """Locate one stop with optimal alignment score.

    Parameters
    ----------
    scomat : ndarray of shape (m + 1, n + 1)
        Main matrix.
    local : bool
        Local or global alignment.
    trail1, trail2 : bool
        Whether trailing gaps of sequence 1 and 2 are free.

    Returns
    -------
    float
        Optimal alignment score.
    ndarray of int of shape (1, 2)
        Coordinates of one alignment stop.

    Notes
    -----
    When there is a tie, the smallest index (row, column) is chosen.

    """
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1

    # local alignment: maximum cell in the matrix
    if local:
        i, j = np.unravel_index(scomat.argmax(), scomat.shape)

    # semi-global alignment
    # free trailing gaps for both: maximum cell in the last column and row
    elif trail1 and trail2:
        i = scomat[:m, n].argmax()  # last column (ends with deletion)
        j = scomat[m, :].argmax()  # last row (ends with insertion)
        if scomat[i, n] >= scomat[m, j]:
            j = n
        else:
            i = m

    # free trailing gaps in seq1: maximum in last row
    elif trail1:
        i, j = m, scomat[m, :].argmax()

    # free trailing gaps in seq2: maximum in last column
    elif trail2:
        i, j = scomat[:, n].argmax(), n

    # global alignment: bottom right cell
    else:
        i, j = m, n

    return scomat[i, j], np.array([[i, j]])


def _all_stops(scomat, local, trail1, trail2, eps=1e-5):
    """Locate all stops with optimal alignment score.

    Parameters
    ----------
    scomat : ndarray of shape (m + 1, n + 1)
        Main matrix.
    local : bool
        Local or global alignment.
    trail1, trail2 : bool
        Whether trailing gaps of sequence 1 and 2 are free.
    eps : float, optional
        Absolute tolerance. Default is 1e-5.

    Returns
    -------
    float
        Optimal alignment score.
    ndarray of int of shape (k, 2)
        Coordinates of all (k) alignment stops.

    Notes
    -----
    Coordinates (row, column) are sorted in ascending order.

    """
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1

    # local alignment
    if local:
        best = scomat.max()
        if eps:
            test = np.isclose(scomat, best, rtol=0, atol=eps)
        else:
            test = scomat == best
        return best, np.argwhere(test)

    # semi-global alignment, free trailing gaps for both
    elif trail1 and trail2:
        col_n, row_m = scomat[:m, n], scomat[m, :]
        best = np.max([col_n.max(), row_m.max()])
        if eps:
            test1 = np.isclose(col_n, best, rtol=0, atol=eps)
            test2 = np.isclose(row_m, best, rtol=0, atol=eps)
        else:
            test1 = col_n == best
            test2 = row_m == best
        ii = np.argwhere(test1)
        jj = np.argwhere(test2)
        return best, np.vstack(
            (
                np.column_stack((ii, np.full(ii.shape[0], n))),
                np.column_stack((np.full(jj.shape[0], m), jj)),
            )
        )

    # free trailing gaps in seq1
    elif trail1:
        row_m = scomat[m, :]
        best = row_m.max()
        jj = np.argwhere(
            np.isclose(row_m, best, rtol=0, atol=eps) if eps else row_m == best
        )
        return best, np.column_stack((np.full(jj.shape[0], m), jj))

    # free trailing gaps in seq2
    elif trail2:
        col_n = scomat[:, n]
        best = col_n.max()
        ii = np.argwhere(
            np.isclose(col_n, best, rtol=0, atol=eps) if eps else col_n == best
        )
        return best, np.column_stack((ii, np.full(ii.shape[0], n)))

    # global alignment
    else:
        return scomat[m, n], np.array([[m, n]])


def _encode_path(path, i0, i1, j0, j1):
    """Perform run-length encoding (RLE) on a dense alignment path.

    Parameters
    ----------
    path : ndarray of uint8 of shape (n_positions,)
        Dense alignment path.
    i0, i1 : int
        Start and stop positions in sequences 1, respectively.
    j0, j1 : int
        Start and stop positions in sequences 2, respectively.

    Returns
    -------
    PairAlignPath
        Encoded alignment path.

    See Also
    --------
    skbio.alignment.AlignPath.from_bits

    """
    if L := path.size:
        segs = np.append(0, np.flatnonzero(path[:-1] != path[1:]) + 1)
        lens = np.append(segs[1:] - segs[:-1], L - segs[-1])
        ints = path[segs]
    else:
        lens = np.array([], dtype=np.intp)
        ints = path
    ranges = np.array([[i0, i1], [j0, j1]], dtype=np.intp)
    return PairAlignPath(lens, ints, ranges=ranges)


def _trailing_gaps(path, pos, i, j, m, n, fill1, fill2):
    """Fill trailing gaps before traceback starts.

    Parameters
    ----------
    path : ndarray of uint8 of shape (n_positions,)
        Dense alignment path.
    pos : int
        Current start position of the path.
    i, j: int
        Current row and column indices in the matrix.
    m, n : int
        Last row and column indices in the matrix (lengths of seq1 and seq2).
    fill1, fill2 : bool
        Fill trailing gaps in sequences 1 and 2, respectively.

    Returns
    -------
    int
        Updated start position of the path.
    int, int
        Stop positions of the path in sequences 1 and 2, respectively.

    """
    # bottom row: ends with insertions (gaps in seq1)
    if fill1 and i == m and j < n:
        pos -= n - j
        path[pos:] = 1
        return pos, m, n

    # right-most column: ends with deletions (gaps in seq2)
    elif fill2 and j == n and i < m:
        pos -= m - i
        path[pos:] = 2
        return pos, m, n

    else:
        return pos, i, j


def _leading_gaps(path, pos, i, j, fill1, fill2):
    """Fill leading gaps after traceback ends.

    Parameters
    ----------
    path : ndarray of uint8 of shape (n_positions,)
        Dense alignment path.
    pos : int
        Current start position of the path.
    i, j: int
        Current row and column indices in the matrix.
    fill1, fill2 : bool
        Fill leading gaps in sequences 1 and 2, respectively.

    Returns
    -------
    int
        Updated start position of the path.
    int, int
        Start positions of the path in sequences 1 and 2, respectively.

    """
    stop = pos

    # top row: starts with insertions (gaps in seq1)
    if fill1 and i == 0 and j > 0:
        pos -= j
        path[pos:stop] = 1
        return pos, 0, 0

    # left-most column: starts with deletions (gaps in seq2)
    elif fill2 and j == 0 and i > 0:
        pos -= i
        path[pos:stop] = 2
        return pos, 0, 0

    else:
        return pos, i, j


def _traceback_one(
    i,
    j,
    matrices,
    gap_open,
    gap_extend,
    local=False,
    lead1=0,
    lead2=0,
    trail1=0,
    trail2=0,
    eps=1e-5,
):
    """Traceback and return one optimal alignment path.

    Parameters
    ----------
    i, j : int
        Stop positions in sequences 1 and 2, respectively.
    matrices : tuple of ndarray of float of shape (m + 1, n + 1)
        Alignment matrices.
    gap_open : float
        Gap opening penalty.
    gap_extend : float
        Gap extension penalty.
    local : bool, optional
        Global (False) or local (True) alignment.
    lead1, lead2, trail1, trail2 : int, optional
        If =2, trim leading / trailing gaps of sequences 1 / 2, respectively.
    eps : float, optional
        Absolute tolerance.

    Returns
    -------
    PairAlignPath
        Optimal alignment path.

    See Also
    --------
    _traceback_all

    Notes
    -----
    This function pre-allocates memory space for the path, which has a maximum length
    of m + n (assuming none of the characters are aligned). Each element represents the
    gap status of one position in the alignment:

        0: no gap
        1: gap in seq1
        2: gap in seq2

    The path will be filled in reverse order (from stop to start, i.e., trace*back*).

    After the path is completed, run-length encoding (RLE) will be performed to
    compress it into a compact form.

    Alternatively, one can encode the path while constructing it, which reduces memory
    consumption. However, that will complicate Cythonization (it will need C++ vector).
    There is no obvious runtime difference between the two methods.

    """
    scomat = matrices[0]
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1

    # current start position of the path, i.e., the index right *after* the next
    # position to be filled.
    pos = m + n

    # allocate space for path
    path = np.empty(pos, dtype=np.uint8)

    # fill trailing gaps (from edge to bottom-right cell).
    pos, i1, j1 = _trailing_gaps(path, pos, i, j, m, n, trail1 < 2, trail2 < 2)

    # Traceback matrix body. This is a time-consuming step. However, it is O(m + n),
    # thus is faster than the matrix filling step.
    if gap_open:
        pos, i, j = _trace_one_affine(
            path, pos, i, j, *matrices, gap_extend, local, eps
        )
    else:
        pos, i, j = _trace_one_linear(path, pos, i, j, scomat, gap_extend, local, eps)

    # fill leading gaps (from top-left cell to edge).
    pos, i0, j0 = _leading_gaps(path, pos, i, j, lead1 < 2, lead2 < 2)

    # encode path
    return _encode_path(path[pos:], i0, i1, j0, j1)


"""Encodings of moving directions during traceback. Columns are:

    0. Row offset (i)
        - Can be calculated with (state & 1) ^ 1

    1. Column offset (j)
        - Can be calculated with (state >> 1) ^ 1

    2. Gap state
        0. Substitution (no gap)
        1. Insertion (gap in seq1)
        2. Deletion (gap in seq2)
        3. Invalid state

    3. Matrix index
        0. Main matrix
        1. Insertion matrix (affine)
        2. Deletion matrix (affine)
        3. Insertion matrix 2 (2-piece affine)
        4. Deletion matrix 2 (2-piece affine)

"""
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
    ]
)


def _linear_moves(i, j, scomat, query, target, gap, eps):
    """Identify move direction(s) at a cell with linear gap penalty.

    Parameters
    ----------
    i, j : int
        Current row and column indices in the matrix, respectively.
    scomat : ndarray of shape (m + 1, n + 1)
        Main matrix.
    query : ndarray of float of shape (m, n_symbols)
        Query profile.
    target : ndarray of int of shape (n,)
        Target sequence.
    gap : float
        Gap penalty.
    eps : float
        Absolute tolerance.

    Returns
    -------
    list of int
        Move directions.

    """
    score = scomat[i, j]
    moves = []
    if (
        abs(scomat[i - 1, j - 1] + query[i - 1, target[j - 1]] - score) <= eps
    ):  # substitution
        moves.append(0)
    if abs(scomat[i, j - 1] - gap - score) <= eps:  # insertion
        moves.append(1)
    if abs(scomat[i - 1, j] - gap - score) <= eps:  # deletion
        moves.append(2)
    return moves


def _affine_moves(i, j, mat, scomat, insmat, delmat, query, target, gap_oe, gap_e, eps):
    """Identify move direction(s) at a cell with affine gap penalty.

    Parameters
    ----------
    i, j : int
        Current row and column indices in the matrix, respectively.
    mat : int
        Current matrix index.
    scomat : ndarray of shape (m + 1, n + 1)
        Main matrix.
    insmat : ndarray of shape (m + 1, n + 1)
        Insertion matrix.
    delmat : ndarray of shape (m + 1, n + 1)
        Deletion matrix.
    query : ndarray of float of shape (m, n_symbols)
        Query profile.
    target : ndarray of int of shape (n,)
        Target sequence.
    gap_oe : float
        Gap opening + extension penalty.
    gap_e : float
        Gap extension penalty.
    eps : float
        Absolute tolerance.

    Returns
    -------
    list of int
        Move directions.

    """
    moves = []

    # main matrix
    if mat == 0:
        score = scomat[i, j]
        if (
            abs(scomat[i - 1, j - 1] + query[i - 1, target[j - 1]] - score) <= eps
        ):  # substitution
            moves.append(0)
        if abs(insmat[i, j] - score) <= eps:  # jump to insertion matrix
            moves.append(5)
        if abs(delmat[i, j] - score) <= eps:  # jump to deletion matrix
            moves.append(6)
    # insertion matrix
    elif mat == 1:
        score = insmat[i, j]
        if abs(scomat[i, j - 1] - gap_oe - score) <= eps:  # open insertion
            moves.append(1)
        if abs(insmat[i, j - 1] - gap_e - score) <= eps:  # extend insertion
            moves.append(3)
    # deletion matrix
    else:
        score = delmat[i, j]
        if abs(scomat[i - 1, j] - gap_oe - score) <= eps:  # open deletion
            moves.append(2)
        if abs(delmat[i - 1, j] - gap_e - score) <= eps:  # extend deletion
            moves.append(4)

    return moves


def _traceback_all(
    stops,
    max_paths,
    matrices,
    query,
    target,
    gap_open,
    gap_extend,
    local=False,
    lead1=0,
    lead2=0,
    trail1=0,
    trail2=0,
    eps=1e-5,
):
    """Traceback and return all optimal alignment paths.

    Parameters
    ----------
    stops : sequence of (int, int)
        Stop positions in both sequences.
    max_paths : int or None
        Maximum number of paths to return, or None to return all paths. Cannot be 0.
    matrices : tuple of ndarray of float of shape (m + 1, n + 1)
        Alignment matrices.
    query : ndarray of float of shape (m, n_symbols)
        Query profile.
    target : ndarray of int of shape (n,)
        Target sequence.
    gap_open : float
        Gap opening penalty.
    gap_extend : float
        Gap extension penalty.
    local : bool, optional
        Global (False) or local (True) alignment.
    lead1, lead2, trail1, trail2 : int, optional
        If =2, trim leading / trailing gaps of sequences 1 / 2, respectively.
    eps : float, optional
        Absolute tolerance.

    Returns
    -------
    list of PairAlignPath
        Optimal alignment paths.

    See Also
    --------
    _traceback_one

    Notes
    -----
    This function enumerates all optimal alignment paths or until a given number of
    paths has been reached.

    A stack is used to store branching paths and to enable depth-first search (DFS).
    For linear gap penalty, the stack's maximum size is m + n, which is equivalent to
    the maximum size of a full path, times (b - 1), where b is the number of possible
    branches (3). For affine gap penalty, the space is appr. (m + n) * 3. Regardless
    of the total number of full paths, which could easily explode, the stack size is
    linear and manageable.

    This function is purely Python and is not as efficient as it could be. Cythonizing
    it is possible but complicated (requiring C++ vector and class). This can be done
    in the future.

    """
    lead1, lead2, trail1, trail2 = lead1 < 2, lead2 < 2, trail1 < 2, trail2 < 2

    scomat = matrices[0]
    m = scomat.shape[0] - 1
    n = scomat.shape[1] - 1
    max_len = m + n

    if gap_open:
        gap_oe = gap_open + gap_extend
        insmat = matrices[1]
        delmat = matrices[2]

    # results (optimal alignment paths)
    paths = []

    # allocate the first path
    path = np.empty(max_len, dtype=np.uint8)

    # iterate over all starting points
    for i, j in stops:
        pos = max_len

        # fill trailing gaps
        pos, i1, j1 = _trailing_gaps(path, pos, i, j, m, n, trail1, trail2)

        # matrix index (start from main matrix (0))
        mat = 0

        # perform DFS to enumerate all paths
        stack = [(path, pos, i, j, mat)]
        while stack:
            path, pos, i, j, mat = stack.pop()

            # Check whether a full path has been completed. There are two scenarios:
            finished = False

            # 1) Local alignment: Path terminates where cell = 0.
            # Note: The matrix filling functions guarantee that in local alignment,
            # cell values cannot be negative. Otherwise this can terminate the loop
            # prematurely.
            if local and abs(scomat[i, j]) <= eps:
                i0, j0 = i, j
                finished = True

            # 2) Global alignment: Reached the top or left-most edge of the matrix.
            # Path will extend straight left or up to the top-left cell.
            elif i == 0 or j == 0:
                pos, i0, j0 = _leading_gaps(path, pos, i, j, lead1, lead2)
                finished = True

            # Create a path object, and halt if the path number limit has been reached.
            if finished:
                paths.append(_encode_path(path[pos:], i0, i1, j0, j1))
                if max_paths and len(paths) == max_paths:
                    break
                else:
                    continue

            # Otherwise, the current cell must be within the main body of the matrix.
            # Check all possible moving directions and get ones that match the current
            # optimal alignment score.
            if gap_open:
                moves = _affine_moves(
                    i,
                    j,
                    mat,
                    scomat,
                    insmat,
                    delmat,
                    query,
                    target,
                    gap_oe,
                    gap_extend,
                    eps,
                )
            else:
                moves = _linear_moves(i, j, scomat, query, target, gap_extend, eps)

            # This is impossible. Raise error for debugging purpose.
            # if not moves:
            #     raise ValueError("Traceback cannot proceed.")

            last = len(moves) - 1
            for k, move in enumerate(moves):
                di, dj, state, mat_ = MOVES[move]

                # If path branches at this cell (i.e., 1 path becomes 2 or 3), we need
                # to create a copy of the path. Otherwise, we can use the same path to
                # save memory and runtime.
                if k < last:
                    path_ = path.copy()
                else:
                    path_ = path

                # Deal with gap state (3 means jumping from main matrix into another
                # matrix without advancing the path).
                if state < 3:
                    pos_ = pos - 1
                    path_[pos_] = state
                else:
                    pos_ = pos

                stack.append((path_, pos_, i - di, j - dj, mat_))

    return paths
