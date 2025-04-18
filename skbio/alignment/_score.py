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
from skbio.alignment import TabularMSA
from ._c_score import _trim_terminal_gaps, _multi_align_score


def _extract_seqs(alignment, gap_chars):
    r"""Extract individual sequences from an alignment.

    Various supported formats can be parsed.

    Returns
    -------
    ndarray of uint8 of shape (n_sequences, n_positions)
        Matrix of stacked sequences as ASCII codes.
    list of int
        Gap characters.

    """
    if isinstance(alignment, TabularMSA):
        seqs = [x._bytes for x in alignment]
        gap_chars = alignment.dtype.gap_chars
    else:
        seqs = [
            (x if isinstance(x, Sequence) else Sequence(x))._bytes for x in alignment
        ]
    if not seqs:
        raise ValueError("There is no sequence in the alignment.")
    try:
        seqs = np.vstack(seqs)
    except ValueError:
        raise ValueError("Sequence lengths do not match.")
    gap_chars = [ord(x) for x in gap_chars]
    return seqs, gap_chars


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
    seqs, gap_chars = _extract_seqs(alignment, gap_chars)
    bits = np.isin(seqs, gap_chars)
    n = seqs.shape[0]
    starts = np.empty(n, dtype=int)
    stops = np.empty(n, dtype=int)
    _trim_terminal_gaps(bits, starts, stops)
    return starts, stops


def align_score(alignment, sub_score, gap_cost, terminal_gaps=False, gap_chars="-."):
    r"""Calculate the alignment score of two or more aligned sequences.

    For two sequences, their pairwise alignment score will be calculated. For three or
    more sequences, the sum-of-pairs (SP) alignment score will be returned.

    Parameters
    ----------
    alignment : TabularMSA, iterable of Sequence or str
        Aligned sequences.

    sub_score : int, float, array_like of (2,), SubstitutionMatrix, or str
        Score of a match, mismatch or substitution. It can be one of the following:

        - Tuple of two numbers: Match score and mismatch score.
        - :class:`~skbio.sequence.SubstitutionMatrix`: A matrix of substitution scores.
        - String: Name of the substitution matrix that can be recognized by
          ``SubstitutionMatrix.by_name``.

    gap_cost : int, float, or array_like of (2,)
        Penalty of a gap. The value is usually positive, representing a subtraction
        from the alignment score. It can be one of the following:

        - One number: Linear gap penalty. Each gap position is penalized by this value
          (*g*). A contiguous gap of length *n* has a total penalty of *g* * *n*.
        - Tuple of two numbers: Affine gap penalty. The two numbers (*o*, *e*)
          represent gap open penalty and gap extension penalty, respectively. A
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
        not a ``TabularMSA`` object, which itself defines gap character(s).

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

    where :math:`o` is the gap open penalty and :math:`e` is the gap extension penalty.

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

    """
    # process sequences
    # TODO: Add support for (path, seqs) structure.
    seqs, gap_chars = _extract_seqs(alignment, gap_chars)
    if (n_seqs := seqs.shape[0]) == 1:
        raise ValueError("There is only one sequence in the alignment.")
    if not seqs.shape[1]:
        raise ValueError("The alignment has a length of 0.")

    # create a bit array representing gaps
    bits = np.isin(seqs, gap_chars)

    # substitution matrix or match/mismatch scores
    if isinstance(sub_score, str):
        sub_score = SubstitutionMatrix.by_name(sub_score)
    if isinstance(sub_score, SubstitutionMatrix):
        # convert sequences into indices in the matrix
        seqs = sub_score._char_hash[seqs]
        # TODO: add a `validate` flag to skip this check
        if (seqs[~bits] >= sub_score.shape[0]).any():
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

    # convert bit array into alignment path
    bits, lens = _get_align_path(bits)

    # identify terminal gaps
    starts = np.empty(n_seqs, dtype=int)
    stops = np.empty(n_seqs, dtype=int)
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
