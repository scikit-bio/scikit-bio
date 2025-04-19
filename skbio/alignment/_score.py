# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numbers import Real

import numpy as np

from skbio.sequence import Sequence, SubstitutionMatrix
from skbio.alignment import TabularMSA, AlignPath
from ._c_score import _trim_terminal_gaps, _multi_align_score


def _seqs_to_bytes(seqs):
    """Convert sequences into bytes."""
    res = []
    res_append = res.append
    for seq in seqs:
        if isinstance(seq, Sequence):
            res_append(seq._bytes)
        elif isinstance(seq, str):
            res_append(np.frombuffer(seq.encode("ascii"), dtype=np.uint8))
        else:
            raise ValueError("Sequences must be strings or Sequence objects.")
    return res


def _parse_alignment(aln, gap_chars):
    """Parse input alignment in any supported format."""
    # 1. tabular alignment
    if isinstance(aln, TabularMSA):
        seqs = [x._bytes for x in aln]
        gap_codes = aln.dtype._gap_codes
        not_path = True

    # 2. alignment path and list of original sequences
    elif len(aln) == 2 and isinstance(aln[0], AlignPath):
        seqs = _seqs_to_bytes(aln[1])
        seqs, gaps, bits, lens = aln[0]._to_matrices(seqs)
        not_path = False

    # 3. list of aligned sequences
    else:
        seqs = _seqs_to_bytes(aln)
        gap_codes = [ord(x) for x in gap_chars]
        not_path = True

    if len(seqs) < 2:
        raise ValueError("Alignment must contain at least two sequences.")

    if not_path:
        try:
            seqs = np.vstack(seqs)
        except ValueError:
            raise ValueError("Sequence lengths do not match.")

    if seqs.shape[1] == 0:
        raise ValueError("The alignment has a length of zero.")

    if not_path:
        gaps = np.isin(seqs, gap_codes)
        bits, lens = _get_align_path(gaps)

    return seqs, gaps, bits, lens


def _get_align_path(bits):
    r"""Calculate the path of an alignment.

    This process is similar to ``AlignPath.from_tabular``, except that bits don't need
    to be packed into uint8's.

    """
    idx = np.append(0, np.where((bits[:, :-1] != bits[:, 1:]).any(axis=0))[0] + 1)
    lens = np.append(idx[1:] - idx[:-1], bits.shape[1] - idx[-1])
    bits = bits[:, idx]
    return bits, lens


def trim_terminal_gaps(alignment, gap_chars="-."):
    r"""Identify the terminal gap-free region of a multiple alignment.

    Parameters
    ----------
    alignment : TabularMSA, iterable of Sequence or str
        Aligned sequences.
    gap_chars : iterable of 1-length str, optional
        Character(s) that represent gaps. Only relevant when ``alignment`` is
        not a ``TabularMSA`` (which itself defines gap character(s)).

    Returns
    -------
    ndarray of int of shape (n_sequences,)
        Start position of terminal gap-free region of each sequence.
        (i.e., index of the first position within non-gap)
    ndarray of int of shape (n_sequences,)
        End position of terminal gap-free region of each sequence.
        (i.e., index of the first position after non-gap)

    Notes
    -----
    If any sequence contains only gaps, the returned start and end positions will both
    be zeros.

    """
    seqs, gaps, _, _ = _parse_alignment(alignment, gap_chars)
    n = seqs.shape[0]
    starts = np.empty(n, dtype=int)
    stops = np.empty(n, dtype=int)
    _trim_terminal_gaps(gaps, starts, stops)
    return starts, stops


def align_score(alignment, sub_score, gap_cost, terminal_gaps=False, gap_chars="-."):
    r"""Calculate the alignment score of two or more aligned sequences.

    For two sequences, their pairwise alignment score will be calculated. For three or
    more sequences, the sum-of-pairs (SP) alignment score will be returned.

    Parameters
    ----------
    alignment : TabularMSA, list of Sequence or str, or (AlignPath, list of Sequence
    or str)
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
          contiguous gap of length *n* has a total penalty of *o* + *e* * (*n* - 1).
          See also notes below.

    terminal_gaps : bool, optional
        Whether gaps at the terminals of the sequences should be penalized. It can be:

        - False (default): Do not penalize terminal gaps. This behavior is known as the
          semi-global (or "glocal") alignment.
        - True: Penalize terminal gaps using the same method defined by ``gap_cost``.
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

       G(n) = o + e \times (n - 1)

    where :math:`o` is the gap opening penalty and :math:`e` is the gap extension
    penalty.

    It should be noted that, discrepancy exists among literature and implementations
    regarding whether gap extension penalty should apply to the first position of a
    gap. scikit-bio's formula is consistent with multiple common alignment tools,
    including EMBOSS, parasail, Biopython and Biotite.

    Meanwhile, multiple other tools, such as BLAST ([1]_), Minimap2, SeqAn3, and
    WFA2-lib, use the following formula instead:

    .. math::

       G(n) = o + e \times n

    Therefore, if you intend to reproduce the behavior of a software tool of the
    latter category using scikit-bio, you will need to add :math:`e` to :math:`o`
    while adopting their parameters. For example, BLASTN's default parameters
    ``o=5, e=2`` will become ``o=7, e=2`` in scikit-bio. Vice versa.

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
    -8.0

    Calculate the score of a multiple alignment of protein sequences, using the
    BLOSUM62 substitution matrix, with gap opening and extension penalties being 11
    and 1 (the default BLASTP parameters). Note that terminal gaps are not penalized
    by default unless ``terminal_gaps`` is set to True.

    >>> msa = TabularMSA([Protein("MKQ-PSV"),
    ...                   Protein("MKIDTS-"),
    ...                   Protein("MVIDPSS")])
    >>> align_score(msa, "BLOSUM62", (11, 1))
    13.0

    """
    # process input alignment
    seqs, gaps, bits, lens = _parse_alignment(alignment, gap_chars)

    # substitution matrix or match/mismatch scores
    if isinstance(sub_score, str):
        sub_score = SubstitutionMatrix.by_name(sub_score)
    if isinstance(sub_score, SubstitutionMatrix):
        # convert sequences into indices in the matrix
        seqs = sub_score._char_hash[seqs]
        # TODO: add a `validate` flag to skip this check
        if (seqs[~gaps] >= sub_score.shape[0]).any():
            raise ValueError(
                "Sequences contain characters that are not present in the provided "
                "substitution matrix."
            )
        submat = sub_score._data
        match, mismatch = 0, 0
    else:
        submat = np.empty((0, 0))
        match, mismatch = sub_score

    # affine or linear gap penalty
    # TODO: Add support for dual affine gap penalty (four numbers).
    if isinstance(gap_cost, Real):
        gap_open, gap_extend = gap_cost, gap_cost
    else:
        gap_open, gap_extend = gap_cost

    # identify terminal gaps
    n = seqs.shape[0]
    starts = np.empty(n, dtype=int)
    stops = np.empty(n, dtype=int)
    _trim_terminal_gaps(bits, starts, stops)
    if not (starts + stops).all():
        raise ValueError("The alignment contains gap-only sequence(s).")

    # TODO: Add support for different treatments of 5' and 3' terminal gaps.
    return _multi_align_score(
        seqs,
        bits,
        lens,
        starts,
        stops,
        submat,
        match,
        mismatch,
        gap_open,
        gap_extend,
        terminal_gaps,
    )
