# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._utils import encode_alignment, prepare_gapcost
from ._cutils import _trim_end_gaps, _multi_align_score


def align_score(alignment, sub_score, gap_cost, free_ends=True, gap_chars="-."):
    r"""Calculate the alignment score of two or more aligned sequences.

    For two sequences, their pairwise alignment score will be calculated. For three or
    more sequences, the sum-of-pairs (SP) alignment score will be returned.

    .. versionadded:: 0.6.4

    Parameters
    ----------
    alignment : TabularMSA, iterable, or (AlignPath, iterable)
        Aligned sequences. Can be any of the following:

        - :class:`~skbio.alignment.TabularMSA`.
        - List of *aligned* sequences as raw strings or ``Sequence`` objects.
        - Tuple of :class:`~skbio.alignment.AlignPath` and the corresponding list of
          *original* (unaligned) sequences.

    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Score of a match, mismatch or substitution. It can be one of the following:

        - Tuple of two numbers: Match score and mismatch score.
        - :class:`~skbio.sequence.SubstitutionMatrix`: A matrix of substitution scores.
        - String: Name of the substitution matrix that can be recognized by
          ``SubstitutionMatrix.by_name``.

    gap_cost : float or tuple of (float, float)
        Penalty of a gap. The value is usually positive, representing a subtraction
        from the alignment score. It can be one of the following:

        - One number: Linear gap penalty. Each gap position is penalized by this value
          (*g*). A contiguous gap of length *n* has a total penalty of *g* * *n*.
        - Tuple of two numbers: Affine gap penalty. The two numbers (*o*, *e*)
          represent gap opening penalty and gap extension penalty, respectively. A
          contiguous gap of length *n* has a total penalty of *o* + *e* * *n*.
          See also notes below.

    free_ends : bool, optional
        Whether gaps at the sequence terminals are free from penalization. It can be:

        - True (default): Do not penalize terminal gaps. This behavior is known as the
          semi-global (or "glocal") alignment.
        - False: Penalize terminal gaps using the same method defined by ``gap_cost``.
          This behavior is the true global alignment.

    gap_chars : iterable of 1-length str, optional
        Character(s) that represent gaps. Only relevant when ``alignment`` is
        a list of aligned sequences.

    Returns
    -------
    float
        Alignment score.

    Raises
    ------
    ValueError
        If there are less than two sequences in the alignment.
    ValueError
        If the alignment has zero length.
    ValueError
        If any sequence in the alignment contains only gaps.
    ValueError
        If any sequence contains characters not present in the designated
        substitution matrix.

    See Also
    --------
    skbio.sequence.SubstitutionMatrix

    Notes
    -----
    Under the affine gap penalty mode, which is the most common situation, the penalty
    of a contiguous gap of length :math:`n` is calculated as:

    .. math::

       G(n) = o + e \times n \tag{1}

    where :math:`o` is the gap opening penalty and :math:`e` is the gap extension
    penalty.

    It should be noted that, discrepancy exists among literature and implementations
    regarding whether gap extension penalty should apply to the first position of a
    gap. scikit-bio's equation is consistent with multiple common alignment tools,
    such as BLAST [1]_, Minimap2, SeqAn3, and WFA2-lib.

    Meanwhile, multiple other tools, such as EMBOSS, parasail, Biopython and Biotite,
    use the following equation instead:

    .. math::

       G(n) = o + e \times (n - 1) \tag{2}

    Therefore, if you intend to reproduce the behavior of a software tool of the
    latter category using scikit-bio, you will need to subtract :math:`e` from
    :math:`o` when adopting its parameters. For example, EMBOSS' default parameters
    ``o=10, e=0.5`` will become ``o=9.5, e=0.5`` in scikit-bio. Vice versa.

    .. versionchanged:: 0.6.4
        Previous alignment algorithms in scikit-bio used Eq. 2. These functions were
        deprecated in 0.5.x and will be removed in 0.6.x. Future functions will
        uniformly use Eq. 1.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK62051/

    Examples
    --------
    >>> from skbio.sequence import DNA, Protein
    >>> from skbio.alignment import TabularMSA, align_score

    Calculate the score of a pair of aligned DNA sequences, with match score = 2,
    mismatch score = -3, gap opening penalty = 5, and gap extension penalty = 2 (the
    default BLASTN parameters).

    >>> seq1 = DNA("CGGTCGTAACGCGTA---CA")
    >>> seq2 = DNA("CAG--GTAAG-CATACCTCA")
    >>> align_score([seq1, seq2], (2, -3), (5, 2))
    -14.0

    Calculate the score of a multiple alignment of protein sequences, using the
    BLOSUM62 substitution matrix, with gap opening and extension penalties being 11
    and 1 (the default BLASTP parameters). Note that terminal gaps are not penalized
    by default unless ``free_ends`` is set to False.

    >>> msa = TabularMSA([Protein("MKQ-PSV"),
    ...                   Protein("MKIDTS-"),
    ...                   Protein("MVIDPSS")])
    >>> align_score(msa, "BLOSUM62", (11, 1))
    11.0

    """
    # process input alignment
    seqs, submat, bits, lens = encode_alignment(alignment, sub_score, gap_chars)
    if (n := seqs.shape[0]) < 2:
        raise ValueError("Alignment must contain at least two sequences.")

    # determine gap penalty method
    gap_open, gap_extend = prepare_gapcost(gap_cost, dtype=submat.dtype.type)

    # identify terminal gaps
    starts, stops = np.empty(n, dtype=int), np.empty(n, dtype=int)
    _trim_end_gaps(bits, starts, stops)
    if not (starts + stops).all():
        raise ValueError("The alignment contains gap-only sequence(s).")

    # TODO: Add support for different treatments of 5' and 3' terminal gaps.
    return _multi_align_score(
        seqs, bits, lens, starts, stops, submat, gap_open, gap_extend, free_ends
    )
