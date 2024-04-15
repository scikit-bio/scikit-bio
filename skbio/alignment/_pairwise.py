# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
from itertools import product

import numpy as np

from skbio.alignment import TabularMSA
from skbio.alignment._ssw_wrapper import StripedSmithWaterman
from skbio.sequence import DNA, RNA, Protein
from skbio.sequence import GrammaredSequence
from skbio.sequence import SubstitutionMatrix
from skbio.util import EfficiencyWarning
from skbio.util._warning import _warn_deprecated


def local_pairwise_align_nucleotide(
    seq1,
    seq2,
    gap_open_penalty=5,
    gap_extend_penalty=2,
    match_score=2,
    mismatch_score=-3,
    substitution_matrix=None,
):
    """Locally align exactly two nucleotide seqs with Smith-Waterman.

    Parameters
    ----------
    seq1 : DNA or RNA
        The first unaligned sequence.
    seq2 : DNA or RNA
        The second unaligned sequence.
    gap_open_penalty : int or float, optional
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float, optional
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    match_score : int or float, optional
        The score to add for a match between a pair of bases (this is added
        to the previous best alignment score, so is typically positive).
    mismatch_score : int or float, optional
        The score to add for a mismatch between a pair of bases (this is
        added to the previous best alignment score, so is typically
        negative).
    substitution_matrix: 2D dict (or similar)
        Lookup for substitution scores (these values are added to the
        previous best alignment score). If provided, this overrides
        ``match_score`` and ``mismatch_score``.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_protein
    global_pairwise_align_nucelotide

    Notes
    -----
    Default ``match_score``, ``mismatch_score``, ``gap_open_penalty`` and
    ``gap_extend_penalty`` parameters are derived from the NCBI BLAST
    Server [1]_.

    References
    ----------
    .. [1] http://blast.ncbi.nlm.nih.gov/Blast.cgi

    """
    for seq in seq1, seq2:
        if not isinstance(seq, (DNA, RNA)):
            raise TypeError(
                "`seq1` and `seq2` must be DNA or RNA, not type %r" % type(seq).__name__
            )

    # use the substitution matrix provided by the user, or compute from
    # match_score and mismatch_score if a substitution matrix was not provided
    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.identity(
            "ACGTU", match_score, mismatch_score
        ).to_dict()

    return local_pairwise_align(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
    )


def local_pairwise_align_protein(
    seq1, seq2, gap_open_penalty=11, gap_extend_penalty=1, substitution_matrix=None
):
    """Locally align exactly two protein seqs with Smith-Waterman.

    Parameters
    ----------
    seq1 : Protein
        The first unaligned sequence.
    seq2 : Protein
        The second unaligned sequence.
    gap_open_penalty : int or float, optional
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float, optional
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    substitution_matrix: 2D dict (or similar), optional
        Lookup for substitution scores (these values are added to the
        previous best alignment score); default is BLOSUM 50.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_nucleotide
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_protein
    global_pairwise_align_nucelotide

    Notes
    -----
    Default ``gap_open_penalty`` and ``gap_extend_penalty`` parameters are
    derived from the NCBI BLAST Server [1]_.

    The BLOSUM (blocks substitution matrices) amino acid substitution matrices
    were originally defined in [2]_.

    References
    ----------
    .. [1] http://blast.ncbi.nlm.nih.gov/Blast.cgi
    .. [2] Amino acid substitution matrices from protein blocks.
       S Henikoff and J G Henikoff.
       Proc Natl Acad Sci U S A. Nov 15, 1992; 89(22): 10915-10919.

    """
    for seq in seq1, seq2:
        if not isinstance(seq, Protein):
            raise TypeError(
                "`seq1` and `seq2` must be Protein, not type %r" % type(seq).__name__
            )

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.by_name("BLOSUM50").to_dict()

    return local_pairwise_align(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
    )


def local_pairwise_align(
    seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
):
    """Locally align exactly two seqs with Smith-Waterman.

    Parameters
    ----------
    seq1 : GrammaredSequence
        The first unaligned sequence.
    seq2 : GrammaredSequence
        The second unaligned sequence.
    gap_open_penalty : int or float
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    substitution_matrix: 2D dict (or similar)
        Lookup for substitution scores (these values are added to the
        previous best alignment score).

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_protein
    global_pairwise_align_nucelotide

    Notes
    -----
    This algorithm was originally described in [1]_. The scikit-bio
    implementation was validated against the EMBOSS water web server [2]_.

    References
    ----------
    .. [1] Identification of common molecular subsequences.
       Smith TF, Waterman MS.
       J Mol Biol. 1981 Mar 25;147(1):195-7.
    .. [2] http://www.ebi.ac.uk/Tools/psa/emboss_water/

    """
    warn(
        "You're using skbio's python implementation of Smith-Waterman "
        "alignment. This will be very slow (e.g., thousands of times slower) "
        "than skbio.alignment.local_pairwise_align_ssw.",
        EfficiencyWarning,
    )

    for seq in seq1, seq2:
        if not isinstance(seq, GrammaredSequence):
            raise TypeError(
                "`seq1` and `seq2` must be %r subclasses, not type %r"
                % (GrammaredSequence.__name__, type(seq).__name__)
            )

    if type(seq1) is not type(seq2):
        raise TypeError(
            "`seq1` and `seq2` must be the same type: %r != %r"
            % (type(seq1).__name__, type(seq2).__name__)
        )

    seq1 = _coerce_alignment_input_type(seq1)
    seq2 = _coerce_alignment_input_type(seq2)

    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        substitution_matrix,
        new_alignment_score=0.0,
        init_matrices_f=_init_matrices_sw,
    )

    end_row_position, end_col_position = np.unravel_index(
        np.argmax(score_matrix), score_matrix.shape
    )

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = _traceback(
        traceback_matrix, score_matrix, seq1, seq2, end_row_position, end_col_position
    )
    start_end_positions = [
        (seq1_start_position, end_col_position - 1),
        (seq2_start_position, end_row_position - 1),
    ]

    msa = TabularMSA(aligned1 + aligned2)

    return msa, score, start_end_positions


def global_pairwise_align_nucleotide(
    seq1,
    seq2,
    gap_open_penalty=5,
    gap_extend_penalty=2,
    match_score=1,
    mismatch_score=-2,
    substitution_matrix=None,
    penalize_terminal_gaps=False,
):
    """Globally align nucleotide seqs or alignments with Needleman-Wunsch.

    Parameters
    ----------
    seq1 : DNA, RNA, or TabularMSA[DNA|RNA]
        The first unaligned sequence(s).
    seq2 : DNA, RNA, or TabularMSA[DNA|RNA]
        The second unaligned sequence(s).
    gap_open_penalty : int or float, optional
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float, optional
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    match_score : int or float, optional
        The score to add for a match between a pair of bases (this is added
        to the previous best alignment score, so is typically positive).
    mismatch_score : int or float, optional
        The score to add for a mismatch between a pair of bases (this is
        added to the previous best alignment score, so is typically
        negative).
    substitution_matrix: 2D dict (or similar)
        Lookup for substitution scores (these values are added to the
        previous best alignment score). If provided, this overrides
        ``match_score`` and ``mismatch_score``.
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_protein

    Notes
    -----
    Default ``match_score``, ``mismatch_score``, ``gap_open_penalty`` and
    ``gap_extend_penalty`` parameters are derived from the NCBI BLAST
    Server [1]_.

    This function can be use to align either a pair of sequences, a pair of
    alignments, or a sequence and an alignment.

    References
    ----------
    .. [1] http://blast.ncbi.nlm.nih.gov/Blast.cgi

    """
    for seq in seq1, seq2:
        if not isinstance(seq, (DNA, RNA, TabularMSA)):
            raise TypeError(
                "`seq1` and `seq2` must be DNA, RNA, or TabularMSA, not type "
                "%r" % type(seq).__name__
            )
        if isinstance(seq, TabularMSA) and not issubclass(seq.dtype, (DNA, RNA)):
            raise TypeError(
                "`seq1` and `seq2` must be TabularMSA with DNA or RNA dtype, "
                "not dtype %r" % seq.dtype.__name__
            )

    # use the substitution matrix provided by the user, or compute from
    # match_score and mismatch_score if a substitution matrix was not provided
    if substitution_matrix is None:
        substitution_matrix = make_identity_substitution_matrix(
            match_score, mismatch_score
        )

    return global_pairwise_align(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        substitution_matrix,
        penalize_terminal_gaps=penalize_terminal_gaps,
    )


def global_pairwise_align_protein(
    seq1,
    seq2,
    gap_open_penalty=11,
    gap_extend_penalty=1,
    substitution_matrix=None,
    penalize_terminal_gaps=False,
):
    """Globally align pair of protein seqs or alignments with Needleman-Wunsch.

    Parameters
    ----------
    seq1 : Protein or TabularMSA[Protein]
        The first unaligned sequence(s).
    seq2 : Protein or TabularMSA[Protein]
        The second unaligned sequence(s).
    gap_open_penalty : int or float, optional
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float, optional
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    substitution_matrix: 2D dict (or similar), optional
        Lookup for substitution scores (these values are added to the
        previous best alignment score); default is BLOSUM 50.
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_nucelotide

    Notes
    -----
    Default ``gap_open_penalty`` and ``gap_extend_penalty`` parameters are
    derived from the NCBI BLAST Server [1]_.

    The BLOSUM (blocks substitution matrices) amino acid substitution matrices
    were originally defined in [2]_.

    This function can be use to align either a pair of sequences, a pair of
    alignments, or a sequence and an alignment.

    References
    ----------
    .. [1] http://blast.ncbi.nlm.nih.gov/Blast.cgi
    .. [2] Amino acid substitution matrices from protein blocks.
       S Henikoff and J G Henikoff.
       Proc Natl Acad Sci U S A. Nov 15, 1992; 89(22): 10915-10919.

    """
    for seq in seq1, seq2:
        if not isinstance(seq, (Protein, TabularMSA)):
            raise TypeError(
                "`seq1` and `seq2` must be Protein or TabularMSA, not type %r"
                % type(seq).__name__
            )
        if isinstance(seq, TabularMSA) and not issubclass(seq.dtype, Protein):
            raise TypeError(
                "`seq1` and `seq2` must be TabularMSA with Protein dtype, "
                "not dtype %r" % seq.dtype.__name__
            )

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.by_name("BLOSUM50").to_dict()

    return global_pairwise_align(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        substitution_matrix,
        penalize_terminal_gaps=penalize_terminal_gaps,
    )


def global_pairwise_align(
    seq1,
    seq2,
    gap_open_penalty,
    gap_extend_penalty,
    substitution_matrix,
    penalize_terminal_gaps=False,
):
    """Globally align a pair of seqs or alignments with Needleman-Wunsch.

    Parameters
    ----------
    seq1 : GrammaredSequence or TabularMSA
        The first unaligned sequence(s).
    seq2 : GrammaredSequence or TabularMSA
        The second unaligned sequence(s).
    gap_open_penalty : int or float
        Penalty for opening a gap (this is substracted from previous best
        alignment score, so is typically positive).
    gap_extend_penalty : int or float
        Penalty for extending a gap (this is substracted from previous best
        alignment score, so is typically positive).
    substitution_matrix: 2D dict (or similar)
        Lookup for substitution scores (these values are added to the
        previous best alignment score).
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.alignment.local_pairwise_align_ssw
    global_pairwise_align_protein
    global_pairwise_align_nucelotide

    Notes
    -----
    This algorithm (in a slightly more basic form) was originally described
    in [1]_. The scikit-bio implementation was validated against the
    EMBOSS needle web server [2]_.

    This function can be use to align either a pair of sequences, a pair of
    alignments, or a sequence and an alignment.

    References
    ----------
    .. [1] A general method applicable to the search for similarities in
       the amino acid sequence of two proteins.
       Needleman SB, Wunsch CD.
       J Mol Biol. 1970 Mar;48(3):443-53.
    .. [2] http://www.ebi.ac.uk/Tools/psa/emboss_needle/

    """
    warn(
        "You're using skbio's python implementation of Needleman-Wunsch "
        "alignment. This is known to be very slow (e.g., thousands of times "
        "slower than a native C implementation). We'll be adding a faster "
        "version soon (see https://github.com/scikit-bio/scikit-bio/issues/"
        "254 to track progress on this).",
        EfficiencyWarning,
    )

    for seq in seq1, seq2:
        # We don't need to check the case where `seq` is a `TabularMSA` with a
        # dtype that isn't a subclass of `GrammaredSequence`, this is
        # guaranteed by `TabularMSA`.
        if not isinstance(seq, (GrammaredSequence, TabularMSA)):
            raise TypeError(
                "`seq1` and `seq2` must be GrammaredSequence subclasses or "
                "TabularMSA, not type %r" % type(seq).__name__
            )

    seq1 = _coerce_alignment_input_type(seq1)
    seq2 = _coerce_alignment_input_type(seq2)

    if seq1.dtype is not seq2.dtype:
        raise TypeError(
            "`seq1` and `seq2` must have the same dtype: %r != %r"
            % (seq1.dtype.__name__, seq2.dtype.__name__)
        )

    if penalize_terminal_gaps:
        init_matrices_f = _init_matrices_nw
    else:
        init_matrices_f = _init_matrices_nw_no_terminal_gap_penalty

    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        substitution_matrix,
        new_alignment_score=-np.inf,
        init_matrices_f=init_matrices_f,
        penalize_terminal_gaps=penalize_terminal_gaps,
    )

    end_row_position = traceback_matrix.shape[0] - 1
    end_col_position = traceback_matrix.shape[1] - 1

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = _traceback(
        traceback_matrix, score_matrix, seq1, seq2, end_row_position, end_col_position
    )
    start_end_positions = [
        (seq1_start_position, end_col_position - 1),
        (seq2_start_position, end_row_position - 1),
    ]

    msa = TabularMSA(aligned1 + aligned2)

    return msa, score, start_end_positions


def local_pairwise_align_ssw(sequence1, sequence2, **kwargs):
    """Align query and target sequences with Striped Smith-Waterman.

    Parameters
    ----------
    sequence1 : DNA, RNA, or Protein
        The first unaligned sequence
    sequence2 : DNA, RNA, or Protein
        The second unaligned sequence
    kwargs : dict
        Additional keyword arguments to pass to ``StripedSmithWaterman``.

    Returns
    -------
    tuple
        ``TabularMSA`` object containing the aligned sequences, alignment score
        (float), and start/end positions of each input sequence (iterable
        of two-item tuples). Note that start/end positions are indexes into the
        unaligned sequences.

    Warnings
    --------
    ``local_pairwise_align_ssw`` is deprecated as of ``0.5.8`` and will be removed in
    favor of more general-purpose and performant aligners. Additional details at
    :repo:`issues/1814`.

    Notes
    -----
    This is a wrapper for the SSW package [1]_.

    For a complete list of optional keyword-arguments that can be provided,
    see ``skbio.alignment.StripedSmithWaterman``.

    The following kwargs will not have any effect: `suppress_sequences`,
    `zero_index`, and `protein`

    If an alignment does not meet a provided filter, `None` will be returned.

    References
    ----------
    .. [1] Zhao, Mengyao, Wan-Ping Lee, Erik P. Garrison, & Gabor T.
       Marth. "SSW Library: An SIMD Smith-Waterman C/C++ Library for
       Applications". PLOS ONE (2013). Web. 11 July 2014.
       http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138

    See Also
    --------
    skbio.alignment.StripedSmithWaterman

    """
    # @deprecated
    _warn_deprecated(
        local_pairwise_align_ssw,
        "0.5.8",
        msg="It will be removed in favor of more general purpose and performant "
        "aligners. Additional details at "
        "https://github.com/scikit-bio/scikit-bio/issues/1814.",
    )

    for seq in sequence1, sequence2:
        if not isinstance(seq, (DNA, RNA, Protein)):
            raise TypeError(
                "`sequence1` and `sequence2` must be DNA, RNA, or Protein, "
                "not type %r" % type(seq).__name__
            )

    if type(sequence1) is not type(sequence2):
        raise TypeError(
            "`sequence1` and `sequence2` must be the same type: %r != %r"
            % (type(sequence1).__name__, type(sequence2).__name__)
        )

    # We need the sequences for `TabularMSA` to make sense, so don't let the
    # user suppress them.
    kwargs["suppress_sequences"] = False
    kwargs["zero_index"] = True

    kwargs["protein"] = False
    if isinstance(sequence1, Protein):
        kwargs["protein"] = True

    query = StripedSmithWaterman(str(sequence1), **kwargs)
    alignment = query(str(sequence2))

    # If there is no cigar, then it has failed a filter. Return None.
    if not alignment.cigar:
        return None

    start_end = None
    if alignment.query_begin != -1:
        start_end = [
            (alignment.query_begin, alignment.query_end),
            (alignment.target_begin, alignment.target_end_optimal),
        ]

    metadata1 = metadata2 = None
    if sequence1.has_metadata():
        metadata1 = sequence1.metadata
    if sequence2.has_metadata():
        metadata2 = sequence2.metadata

    constructor = type(sequence1)
    msa = TabularMSA(
        [
            constructor(
                alignment.aligned_query_sequence, metadata=metadata1, validate=False
            ),
            constructor(
                alignment.aligned_target_sequence, metadata=metadata2, validate=False
            ),
        ]
    )

    return msa, alignment.optimal_alignment_score, start_end


def make_identity_substitution_matrix(match_score, mismatch_score, alphabet="ACGTU"):
    """Generate substitution matrix where all matches are scored equally.

    Parameters
    ----------
    match_score : int, float
        The score that should be assigned for all matches. This value is
        typically positive.
    mismatch_score : int, float
        The score that should be assigned for all mismatches. This value is
        typically negative.
    alphabet : iterable of str, optional
        The characters that should be included in the substitution matrix.

    Returns
    -------
    dict of dicts
        All characters in alphabet are keys in both dictionaries, so that any
        pair of characters can be looked up to get their match or mismatch
        score.

    Warnings
    --------
    ``make_identity_substitution_matrix`` is deprecated as of ``0.4.0``. It has been
    replaced by a SubstitutionMatrix class. Additional details at :repo:`pull/1913`.

    """
    # @deprecated
    _warn_deprecated(
        make_identity_substitution_matrix,
        "0.4.0",
        msg="It has been "
        "replaced by the SubstitutionMatrix class. Additional "
        "details at "
        "https://github.com/scikit-bio/scikit-bio/pull/1913.",
    )

    result = {}
    for c1 in alphabet:
        row = {}
        for c2 in alphabet:
            if c1 == c2:
                row[c2] = match_score
            else:
                row[c2] = mismatch_score
        result[c1] = row
    return result


# Functions from here allow for generalized (global or local) alignment. I
# will likely want to put these in a single object to make the naming a little
# less clunky.


def _coerce_alignment_input_type(seq):
    if isinstance(seq, GrammaredSequence):
        return TabularMSA([seq])
    else:
        return seq


_traceback_encoding = {
    "match": 1,
    "vertical-gap": 2,
    "horizontal-gap": 3,
    "uninitialized": -1,
    "alignment-end": 0,
}


def _init_matrices_sw(aln1, aln2, gap_open_penalty, gap_extend_penalty):
    shape = (aln2.shape.position + 1, aln1.shape.position + 1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=int)
    traceback_matrix += _traceback_encoding["uninitialized"]
    traceback_matrix[0, :] = _traceback_encoding["alignment-end"]
    traceback_matrix[:, 0] = _traceback_encoding["alignment-end"]
    return score_matrix, traceback_matrix


def _init_matrices_nw(aln1, aln2, gap_open_penalty, gap_extend_penalty):
    shape = (aln2.shape.position + 1, aln1.shape.position + 1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=int)
    traceback_matrix += _traceback_encoding["uninitialized"]
    traceback_matrix[0, 0] = _traceback_encoding["alignment-end"]

    # cache some values for quicker access
    vgap = _traceback_encoding["vertical-gap"]
    hgap = _traceback_encoding["horizontal-gap"]

    for i in range(1, shape[0]):
        score_matrix[i, 0] = -gap_open_penalty - ((i - 1) * gap_extend_penalty)
        traceback_matrix[i, 0] = vgap

    for i in range(1, shape[1]):
        score_matrix[0, i] = -gap_open_penalty - ((i - 1) * gap_extend_penalty)
        traceback_matrix[0, i] = hgap

    return score_matrix, traceback_matrix


def _init_matrices_nw_no_terminal_gap_penalty(
    aln1, aln2, gap_open_penalty, gap_extend_penalty
):
    shape = (aln2.shape.position + 1, aln1.shape.position + 1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=int)
    traceback_matrix += _traceback_encoding["uninitialized"]
    traceback_matrix[0, 0] = _traceback_encoding["alignment-end"]

    # cache some values for quicker access
    vgap = _traceback_encoding["vertical-gap"]
    hgap = _traceback_encoding["horizontal-gap"]

    for i in range(1, shape[0]):
        traceback_matrix[i, 0] = vgap

    for i in range(1, shape[1]):
        traceback_matrix[0, i] = hgap

    return score_matrix, traceback_matrix


def _compute_substitution_score(
    aln1_chars, aln2_chars, substitution_matrix, gap_substitution_score, gap_chars
):
    substitution_score = 0
    for aln1_char, aln2_char in product(aln1_chars, aln2_chars):
        if aln1_char in gap_chars or aln2_char in gap_chars:
            substitution_score += gap_substitution_score
        else:
            try:
                substitution_score += substitution_matrix[aln1_char][aln2_char]
            except KeyError:
                offending_chars = [
                    c for c in (aln1_char, aln2_char) if c not in substitution_matrix
                ]
                raise ValueError(
                    "One of the sequences contains a character that is "
                    "not contained in the substitution matrix. Are you "
                    "using an appropriate substitution matrix for your "
                    "sequence type (e.g., a nucleotide substitution "
                    "matrix does not make sense for aligning protein "
                    "sequences)? Does your sequence contain invalid "
                    "characters? The offending character(s) is: "
                    " %s." % ", ".join(offending_chars)
                )
    substitution_score /= len(aln1_chars) * len(aln2_chars)
    return substitution_score


def _compute_score_and_traceback_matrices(
    aln1,
    aln2,
    gap_open_penalty,
    gap_extend_penalty,
    substitution_matrix,
    new_alignment_score=-np.inf,
    init_matrices_f=_init_matrices_nw,
    penalize_terminal_gaps=True,
    gap_substitution_score=0,
):
    """Return dynamic programming (score) and traceback matrices.

    A note on the ``penalize_terminal_gaps`` parameter. When this value is
    ``False``, this function is no longer true Smith-Waterman/Needleman-Wunsch
    scoring, but when ``True`` it can result in biologically irrelevant
    artifacts in Needleman-Wunsch (global) alignments. Specifically, if one
    sequence is longer than the other (e.g., if aligning a primer sequence to
    an amplification product, or searching for a gene in a genome) the shorter
    sequence will have a long gap inserted. The parameter is ``True`` by
    default (so that this function computes the score and traceback matrices as
    described by the original authors) but the global alignment wrappers pass
    ``False`` by default, so that the global alignment API returns the result
    that users are most likely to be looking for.

    """
    aln1_length = aln1.shape.position
    aln2_length = aln2.shape.position
    # cache some values for quicker/simpler access
    aend = _traceback_encoding["alignment-end"]
    match = _traceback_encoding["match"]
    vgap = _traceback_encoding["vertical-gap"]
    hgap = _traceback_encoding["horizontal-gap"]

    new_alignment_score = (new_alignment_score, aend)

    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    score_matrix, traceback_matrix = init_matrices_f(
        aln1, aln2, gap_open_penalty, gap_extend_penalty
    )

    # Iterate over the characters in aln2 (which corresponds to the vertical
    # sequence in the matrix)
    for aln2_pos, aln2_chars in enumerate(aln2.iter_positions(ignore_metadata=True), 1):
        aln2_chars = str(aln2_chars)

        # Iterate over the characters in aln1 (which corresponds to the
        # horizontal sequence in the matrix)
        for aln1_pos, aln1_chars in enumerate(
            aln1.iter_positions(ignore_metadata=True), 1
        ):
            aln1_chars = str(aln1_chars)

            # compute the score for a match/mismatch
            substitution_score = _compute_substitution_score(
                aln1_chars,
                aln2_chars,
                substitution_matrix,
                gap_substitution_score,
                aln1.dtype.gap_chars,
            )

            diag_score = (
                score_matrix[aln2_pos - 1, aln1_pos - 1] + substitution_score,
                match,
            )

            # compute the score for adding a gap in aln2 (vertical)
            if not penalize_terminal_gaps and (aln1_pos == aln1_length):
                # we've reached the end of aln1, so adding vertical gaps
                # (which become gaps in aln1) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                up_score = (score_matrix[aln2_pos - 1, aln1_pos], vgap)
            elif traceback_matrix[aln2_pos - 1, aln1_pos] == vgap:
                # gap extend, because the cell above was also a gap
                up_score = (
                    score_matrix[aln2_pos - 1, aln1_pos] - gap_extend_penalty,
                    vgap,
                )
            else:
                # gap open, because the cell above was not a gap
                up_score = (
                    score_matrix[aln2_pos - 1, aln1_pos] - gap_open_penalty,
                    vgap,
                )

            # compute the score for adding a gap in aln1 (horizontal)
            if not penalize_terminal_gaps and (aln2_pos == aln2_length):
                # we've reached the end of aln2, so adding horizontal gaps
                # (which become gaps in aln2) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                left_score = (score_matrix[aln2_pos, aln1_pos - 1], hgap)
            elif traceback_matrix[aln2_pos, aln1_pos - 1] == hgap:
                # gap extend, because the cell to the left was also a gap
                left_score = (
                    score_matrix[aln2_pos, aln1_pos - 1] - gap_extend_penalty,
                    hgap,
                )
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (
                    score_matrix[aln2_pos, aln1_pos - 1] - gap_open_penalty,
                    hgap,
                )

            # identify the largest score, and use that information to populate
            # the score and traceback matrices
            best_score = _first_largest(
                [new_alignment_score, left_score, diag_score, up_score]
            )
            score_matrix[aln2_pos, aln1_pos] = best_score[0]
            traceback_matrix[aln2_pos, aln1_pos] = best_score[1]

    return score_matrix, traceback_matrix


def _traceback(traceback_matrix, score_matrix, aln1, aln2, start_row, start_col):
    # cache some values for simpler reference
    aend = _traceback_encoding["alignment-end"]
    match = _traceback_encoding["match"]
    vgap = _traceback_encoding["vertical-gap"]
    hgap = _traceback_encoding["horizontal-gap"]
    gap_character = aln1.dtype.default_gap_char

    # initialize the result alignments
    aln1_sequence_count = aln1.shape.sequence
    aligned_seqs1 = [[] for e in range(aln1_sequence_count)]

    aln2_sequence_count = aln2.shape.sequence
    aligned_seqs2 = [[] for e in range(aln2_sequence_count)]

    current_row = start_row
    current_col = start_col

    best_score = score_matrix[current_row, current_col]
    current_value = None

    while current_value != aend:
        current_value = traceback_matrix[current_row, current_col]

        if current_value == match:
            for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
                aligned_seq.append(str(input_seq[current_col - 1]))
            for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
                aligned_seq.append(str(input_seq[current_row - 1]))
            current_row -= 1
            current_col -= 1
        elif current_value == vgap:
            for aligned_seq in aligned_seqs1:
                aligned_seq.append(gap_character)
            for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
                aligned_seq.append(str(input_seq[current_row - 1]))
            current_row -= 1
        elif current_value == hgap:
            for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
                aligned_seq.append(str(input_seq[current_col - 1]))
            for aligned_seq in aligned_seqs2:
                aligned_seq.append(gap_character)
            current_col -= 1
        elif current_value == aend:
            continue
        else:
            raise ValueError("Invalid value in traceback matrix: %s" % current_value)

    for i, (aligned_seq, original) in enumerate(zip(aligned_seqs1, aln1)):
        aligned_seq = "".join(aligned_seq)[::-1]
        constructor = aln1.dtype
        metadata = None
        if original.has_metadata():
            metadata = original.metadata
        aligned_seqs1[i] = constructor(aligned_seq, metadata=metadata, validate=False)

    for i, (aligned_seq, original) in enumerate(zip(aligned_seqs2, aln2)):
        aligned_seq = "".join(aligned_seq)[::-1]
        constructor = aln2.dtype
        metadata = None
        if original.has_metadata():
            metadata = original.metadata
        aligned_seqs2[i] = constructor(aligned_seq, metadata=metadata, validate=False)

    return aligned_seqs1, aligned_seqs2, best_score, current_col, current_row


def _first_largest(scores):
    """Similar to max, but returns the first element achieving the high score.

    If max receives a tuple, it will break a tie for the highest value
    of entry[i] with entry[i+1]. We don't want that here - to better match
    with the results of other tools, we want to be able to define which
    entry is returned in the case of a tie.
    """
    result = scores[0]
    for score, direction in scores[1:]:
        if score > result[0]:
            result = (score, direction)
    return result
