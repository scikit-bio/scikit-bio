# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Union, Tuple, TYPE_CHECKING

import numpy as np

from ._utils import encode_alignment, prep_gapcost
from ._cutils import _trim_end_gaps, _multi_align_score

if TYPE_CHECKING:  # pragma: no cover
    from ._utils import AlignmentLike
    from skbio.sequence import SubstitutionMatrix


def align_score(
    alignment: "AlignmentLike",
    sub_score: Union[Tuple[float, float], "SubstitutionMatrix", str] = (1.0, -1.0),
    gap_cost: Union[float, Tuple[float, float]] = 2.0,
    free_ends: bool = True,
    gap_chars: str = "-.",
) -> float:
    r"""Calculate the alignment score of two or more aligned sequences.

    For two sequences, their pairwise alignment score will be calculated. For three or
    more sequences, the sum-of-pairs (SP) alignment score will be returned.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    alignment : TabularMSA, iterable, or (AlignPath, iterable)
        Aligned sequences. Can be any of the following:

        - :class:`~skbio.alignment.TabularMSA`.
        - List of *aligned* sequences as raw strings or ``Sequence`` objects.
        - Tuple of :class:`~skbio.alignment.AlignPath` and the corresponding list of
          *original* (unaligned) sequences.

    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Score of a substitution. May be two numbers (match, mismatch), a substitution
        matrix, or its name. See :func:`pair_align` for details. Default is
        (1.0, -1.0).

    gap_cost : float or tuple of (float, float)
        Penalty of a gap. May be one (linear) or two numbers (affine). See
        :func:`pair_align` for details. Default is 2.0.

    free_ends : bool, optional
        If True (default), gaps at the sequence terminals are free from penalization.

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
    pair_align
    skbio.sequence.SubstitutionMatrix

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
    gap_open, gap_extend = prep_gapcost(gap_cost, dtype=submat.dtype.type)

    # identify terminal gaps
    starts, stops = np.empty(n, dtype=int), np.empty(n, dtype=int)
    _trim_end_gaps(bits, starts, stops)
    if not (starts + stops).all():
        raise ValueError("The alignment contains gap-only sequence(s).")

    return _multi_align_score(
        seqs, bits, lens, starts, stops, submat, gap_open, gap_extend, free_ends
    )
