# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range, zip
from warnings import warn

import numpy as np

from skbio.core.warning import EfficiencyWarning
from .alignment import Alignment
from skbio import BiologicalSequence

# This is temporary: blosum50 does not exist in skbio yet as per
# issue 161. When the issue is resolved, this should be removed in favor
# of an import.
blosum50 = \
    {
        '*': {'*': 1, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5,
              'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5,
              'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5,
              'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5},
        'A': {'*': -5, 'A': 5, 'C': -1, 'B': -2, 'E': -1, 'D': -2, 'G': 0,
              'F': -3, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -2,
              'N': -1, 'Q': -1, 'P': -1, 'S': 1, 'R': -2, 'T': 0, 'W': -3,
              'V': 0, 'Y': -2, 'X': -1, 'Z': -1},
        'C': {'*': -5, 'A': -1, 'C': 13, 'B': -3, 'E': -3, 'D': -4,
              'G': -3, 'F': -2, 'I': -2, 'H': -3, 'K': -3, 'M': -2,
              'L': -2, 'N': -2, 'Q': -3, 'P': -4, 'S': -1, 'R': -4,
              'T': -1, 'W': -5, 'V': -1, 'Y': -3, 'X': -1, 'Z': -3},
        'B': {'*': -5, 'A': -2, 'C': -3, 'B': 6, 'E': 1, 'D': 6, 'G': -1,
              'F': -4, 'I': -4, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 5,
              'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -5, 'V': -3,
              'Y': -3, 'X': -1, 'Z': 1},
        'E': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 6, 'D': 2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': -1, 'R': 0, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5},
        'D': {'*': -5, 'A': -2, 'C': -4, 'B': 6, 'E': 2, 'D': 8, 'G': -1,
              'F': -5, 'I': -4, 'H': -1, 'K': -1, 'M': -4, 'L': -4, 'N': 2,
              'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -5, 'V': -4,
              'Y': -3, 'X': -1, 'Z': 1},
        'G': {'*': -5, 'A': 0, 'C': -3, 'B': -1, 'E': -3, 'D': -1, 'G': 8,
              'F': -4, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0,
              'Q': -2, 'P': -2, 'S': 0, 'R': -3, 'T': -2, 'W': -3, 'V': -4,
              'Y': -3, 'X': -1, 'Z': -2},
        'F': {'*': -5, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -5,
              'G': -4, 'F': 8, 'I': 0, 'H': -1, 'K': -4, 'M': 0, 'L': 1,
              'N': -4, 'Q': -4, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 1,
              'V': -1, 'Y': 4, 'X': -1, 'Z': -4},
        'I': {'*': -5, 'A': -1, 'C': -2, 'B': -4, 'E': -4, 'D': -4,
              'G': -4, 'F': 0, 'I': 5, 'H': -4, 'K': -3, 'M': 2, 'L': 2,
              'N': -3, 'Q': -3, 'P': -3, 'S': -3, 'R': -4, 'T': -1,
              'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -3},
        'H': {'*': -5, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2,
              'F': -1, 'I': -4, 'H': 10, 'K': 0, 'M': -1, 'L': -3, 'N': 1,
              'Q': 1, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -3, 'V': -4,
              'Y': 2, 'X': -1, 'Z': 0},
        'K': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 6, 'M': -2, 'L': -3, 'N': 0,
              'Q': 2, 'P': -1, 'S': 0, 'R': 3, 'T': -1, 'W': -3, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 1},
        'M': {'*': -5, 'A': -1, 'C': -2, 'B': -3, 'E': -2, 'D': -4,
              'G': -3, 'F': 0, 'I': 2, 'H': -1, 'K': -2, 'M': 7, 'L': 3,
              'N': -2, 'Q': 0, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -1,
              'V': 1, 'Y': 0, 'X': -1, 'Z': -1},
        'L': {'*': -5, 'A': -2, 'C': -2, 'B': -4, 'E': -3, 'D': -4,
              'G': -4, 'F': 1, 'I': 2, 'H': -3, 'K': -3, 'M': 3, 'L': 5,
              'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -1,
              'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3},
        'N': {'*': -5, 'A': -1, 'C': -2, 'B': 5, 'E': 0, 'D': 2, 'G': 0,
              'F': -4, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -4, 'N': 7,
              'Q': 0, 'P': -2, 'S': 1, 'R': -1, 'T': 0, 'W': -4, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 0},
        'Q': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2,
              'F': -4, 'I': -3, 'H': 1, 'K': 2, 'M': 0, 'L': -2, 'N': 0,
              'Q': 7, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -1, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 4},
        'P': {'*': -5, 'A': -1, 'C': -4, 'B': -2, 'E': -1, 'D': -1,
              'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -3,
              'L': -4, 'N': -2, 'Q': -1, 'P': 10, 'S': -1, 'R': -3,
              'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': -1},
        'S': {'*': -5, 'A': 1, 'C': -1, 'B': 0, 'E': -1, 'D': 0, 'G': 0,
              'F': -3, 'I': -3, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1,
              'Q': 0, 'P': -1, 'S': 5, 'R': -1, 'T': 2, 'W': -4, 'V': -2,
              'Y': -2, 'X': -1, 'Z': 0},
        'R': {'*': -5, 'A': -2, 'C': -4, 'B': -1, 'E': 0, 'D': -2, 'G': -3,
              'F': -3, 'I': -4, 'H': 0, 'K': 3, 'M': -2, 'L': -3, 'N': -1,
              'Q': 1, 'P': -3, 'S': -1, 'R': 7, 'T': -1, 'W': -3, 'V': -3,
              'Y': -1, 'X': -1, 'Z': 0},
        'T': {'*': -5, 'A': 0, 'C': -1, 'B': 0, 'E': -1, 'D': -1, 'G': -2,
              'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0,
              'Q': -1, 'P': -1, 'S': 2, 'R': -1, 'T': 5, 'W': -3, 'V': 0,
              'Y': -2, 'X': -1, 'Z': -1},
        'W': {'*': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -3, 'D': -5,
              'G': -3, 'F': 1, 'I': -3, 'H': -3, 'K': -3, 'M': -1, 'L': -2,
              'N': -4, 'Q': -1, 'P': -4, 'S': -4, 'R': -3, 'T': -3,
              'W': 15, 'V': -3, 'Y': 2, 'X': -1, 'Z': -2},
        'V': {'*': -5, 'A': 0, 'C': -1, 'B': -3, 'E': -3, 'D': -4, 'G': -4,
              'F': -1, 'I': 4, 'H': -4, 'K': -3, 'M': 1, 'L': 1, 'N': -3,
              'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 5,
              'Y': -1, 'X': -1, 'Z': -3},
        'Y': {'*': -5, 'A': -2, 'C': -3, 'B': -3, 'E': -2, 'D': -3,
              'G': -3, 'F': 4, 'I': -1, 'H': 2, 'K': -2, 'M': 0, 'L': -1,
              'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -1, 'T': -2, 'W': 2,
              'V': -1, 'Y': 8, 'X': -1, 'Z': -2},
        'X': {'*': -5, 'A': -1, 'C': -1, 'B': -1, 'E': -1, 'D': -1,
              'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1,
              'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1,
              'T': -1, 'W': -1, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1},
        'Z': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 5, 'D': 1, 'G': -2,
              'F': -4, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0,
              'Q': 4, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -2, 'V': -3,
              'Y': -2, 'X': -1, 'Z': 5}}


def local_pairwise_align_nucleotide(seq1, seq2, gap_open_penalty=5,
                                    gap_extend_penalty=2,
                                    match_score=2, mismatch_score=-3,
                                    substitution_matrix=None):
    """Locally align exactly two nucleotide seqs with Smith-Waterman

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    skbio.Alignment
        ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    skbio.core.alignment.local_pairwise_align_ssw
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
    # use the substitution matrix provided by the user, or compute from
    # match_score and mismatch_score if a substitution matrix was not provided
    if substitution_matrix is None:
        substitution_matrix = \
            _make_nt_substitution_matrix(match_score, mismatch_score)
    else:
        pass

    return local_pairwise_align(seq1, seq2, gap_open_penalty,
                                gap_extend_penalty, substitution_matrix)


def local_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                 gap_extend_penalty=1,
                                 substitution_matrix=None):
    """Locally align exactly two protein seqs with Smith-Waterman

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    skbio.Alignment
        ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_nucleotide
    skbio.core.alignment.local_pairwise_align_ssw
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
    if substitution_matrix is None:
        substitution_matrix = blosum50

    return local_pairwise_align(seq1, seq2, gap_open_penalty,
                                gap_extend_penalty, substitution_matrix)


def local_pairwise_align(seq1, seq2, gap_open_penalty,
                         gap_extend_penalty, substitution_matrix):
    """Locally align exactly two seqs with Smith-Waterman

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    skbio.Alignment
       ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.core.alignment.local_pairwise_align_ssw
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
    warn("You're using skbio's python implementation of Smith-Waterman "
         "alignment. This will be very slow (e.g., thousands of times slower) "
         "than skbio.core.alignment.local_pairwise_align_ssw.",
         EfficiencyWarning)

    if isinstance(seq1, unicode) or isinstance(seq1, str):
        seq1 = Alignment([BiologicalSequence(seq1)])
    elif isinstance(seq1, BiologicalSequence):
        seq1 = Alignment([seq1])
    else:
        pass

    if isinstance(seq2, unicode) or isinstance(seq2, str):
        seq2 = Alignment([BiologicalSequence(seq2)])
    elif isinstance(seq2, BiologicalSequence):
        seq2 = Alignment([seq2])
    else:
        pass

    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1, seq2, gap_open_penalty, gap_extend_penalty,
        substitution_matrix, new_alignment_score=0.0,
        init_matrices_f=_init_matrices_sw)

    end_row_position, end_col_position =\
        np.unravel_index(np.argmax(score_matrix), score_matrix.shape)

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = \
        _traceback(traceback_matrix, score_matrix, seq1, seq2,
                   end_row_position, end_col_position)
    start_end_positions = [(seq1_start_position, end_col_position-1),
                           (seq2_start_position, end_row_position-1)]

    return Alignment(aligned1 + aligned2, score=score,
                      start_end_positions=start_end_positions)

def global_pairwise_align_nucleotide(seq1, seq2, gap_open_penalty=5,
                                     gap_extend_penalty=2,
                                     match_score=1, mismatch_score=-2,
                                     substitution_matrix=None,
                                     penalize_terminal_gaps=False):
    """Globally align exactly two nucleotide seqs with Needleman-Wunsch

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    skbio.Alignment
        ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.core.alignment.local_pairwise_align_ssw
    global_pairwise_align
    global_pairwise_align_protein

    Notes
    -----
    Default ``match_score``, ``mismatch_score``, ``gap_open_penalty`` and
    ``gap_extend_penalty`` parameters are derived from the NCBI BLAST
    Server [1]_.

    References
    ----------
    .. [1] http://blast.ncbi.nlm.nih.gov/Blast.cgi

    """
    # use the substitution matrix provided by the user, or compute from
    # match_score and mismatch_score if a substitution matrix was not provided
    if substitution_matrix is None:
        substitution_matrix = \
            _make_nt_substitution_matrix(match_score, mismatch_score)
    else:
        pass

    return global_pairwise_align(seq1, seq2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix,
                                 penalize_terminal_gaps=penalize_terminal_gaps)


def global_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                  gap_extend_penalty=1,
                                  substitution_matrix=None,
                                  penalize_terminal_gaps=False):
    """Globally align exactly two protein seqs with Needleman-Wunsch

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    skbio.Alignment
        ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.core.alignment.local_pairwise_align_ssw
    global_pairwise_align
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
    if substitution_matrix is None:
        substitution_matrix = blosum50

    return global_pairwise_align(seq1, seq2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix,
                                 penalize_terminal_gaps=penalize_terminal_gaps)


def global_pairwise_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
                          substitution_matrix, penalize_terminal_gaps=False):
    """Globally align exactly two seqs with Needleman-Wunsch

    Parameters
    ----------
    seq1 : str or BiologicalSequence
        The first unaligned sequence.
    seq2 : str or BiologicalSequence
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
    penalize_terminal_gaps: bool, optional
        If True, will continue to penalize gaps even after one sequence has
        been aligned through its end. This behavior is true Needleman-Wunsch
        alignment, but results in (biologically irrelevant) artifacts when
        the sequences being aligned are of different length. This is ``False``
        by default, which is very likely to be the behavior you want in all or
        nearly all cases.

    Returns
    -------
    skbio.Alignment
        ``Alignment`` object containing the aligned sequences as well as
        details about the alignment.

    See Also
    --------
    local_pairwise_align
    local_pairwise_align_protein
    local_pairwise_align_nucleotide
    skbio.core.alignment.local_pairwise_align_ssw
    global_pairwise_align_protein
    global_pairwise_align_nucelotide

    Notes
    -----
    This algorithm (in a slightly more basic form) was originally described
    in [1]_. The scikit-bio implementation was validated against the
    EMBOSS needle web server [2]_.

    References
    ----------
    .. [1] A general method applicable to the search for similarities in
       the amino acid sequence of two proteins.
       Needleman SB, Wunsch CD.
       J Mol Biol. 1970 Mar;48(3):443-53.
    .. [2] http://www.ebi.ac.uk/Tools/psa/emboss_needle/

    """
    warn("You're using skbio's python implementation of Needleman-Wunsch "
         "alignment. This is known to be very slow (e.g., thousands of times "
         "slower than a native C implementation). We'll be adding a faster "
         "version soon (see https://github.com/biocore/scikit-bio/issues/254 "
         "to track progress on this).", EfficiencyWarning)

    if isinstance(seq1, unicode) or isinstance(seq1, str):
        seq1 = Alignment([BiologicalSequence(seq1)])
    elif isinstance(seq1, BiologicalSequence):
        seq1 = Alignment([seq1])
    else:
        pass

    if isinstance(seq2, unicode) or isinstance(seq2, str):
        seq2 = Alignment([BiologicalSequence(seq2)])
    elif isinstance(seq2, BiologicalSequence):
        seq2 = Alignment([seq2])
    else:
        pass

    if penalize_terminal_gaps:
        init_matrices_f = _init_matrices_nw
    else:
        init_matrices_f = _init_matrices_nw_no_terminal_gap_penalty

    score_matrix, traceback_matrix = \
        _compute_score_and_traceback_matrices(
            seq1, seq2, gap_open_penalty, gap_extend_penalty,
            substitution_matrix, new_alignment_score=-np.inf,
            init_matrices_f=init_matrices_f,
            penalize_terminal_gaps=penalize_terminal_gaps)

    end_row_position = traceback_matrix.shape[0] - 1
    end_col_position = traceback_matrix.shape[1] - 1

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = \
        _traceback(traceback_matrix, score_matrix, seq1, seq2,
                   end_row_position, end_col_position)
    start_end_positions = [(seq1_start_position, end_col_position-1),
                           (seq2_start_position, end_row_position-1)]

    return Alignment(aligned1 + aligned2, score=score,
                      start_end_positions=start_end_positions)

# Functions from here allow for generalized (global or local) alignment. I
# will likely want to put these in a single object to make the naming a little
# less clunky.

_traceback_encoding = {'match': 1, 'vertical-gap': 2, 'horizontal-gap': 3,
                       'uninitialized': -1, 'alignment-end': 0}


def _get_seq_ids(seq1, seq2):
    return _get_seq_id(seq1, "0"), _get_seq_id(seq2, "1")


def _get_seq_id(seq, default_id):
    try:
        result = seq.id
    except AttributeError:
        result = default_id
    else:
        if result is None or result.strip() == "":
            result = default_id
    return result


def _make_nt_substitution_matrix(match_score, mismatch_score, alphabet='ACGT'):
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


def _init_matrices_sw(seq1, seq2, gap_open_penalty, gap_extend_penalty):
    shape = (seq2.sequence_length()+1, seq1.sequence_length()+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, :] = _traceback_encoding['alignment-end']
    traceback_matrix[:, 0] = _traceback_encoding['alignment-end']
    return score_matrix, traceback_matrix


def _init_matrices_nw(seq1, seq2, gap_open_penalty, gap_extend_penalty):
    shape = (seq2.sequence_length()+1, seq1.sequence_length()+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, 0] = _traceback_encoding['alignment-end']

    # cache some values for quicker access
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    for i in range(1, shape[0]):
        score_matrix[i, 0] = -gap_open_penalty - ((i-1) * gap_extend_penalty)
        traceback_matrix[i, 0] = vgap

    for i in range(1, shape[1]):
        score_matrix[0, i] = -gap_open_penalty - ((i-1) * gap_extend_penalty)
        traceback_matrix[0, i] = hgap

    return score_matrix, traceback_matrix


def _init_matrices_nw_no_terminal_gap_penalty(
        seq1, seq2, gap_open_penalty, gap_extend_penalty):
    shape = (seq2.sequence_length()+1, seq1.sequence_length()+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, 0] = _traceback_encoding['alignment-end']

    # cache some values for quicker access
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    for i in range(1, shape[0]):
        traceback_matrix[i, 0] = vgap

    for i in range(1, shape[1]):
        traceback_matrix[0, i] = hgap

    return score_matrix, traceback_matrix


def _compute_score_and_traceback_matrices(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix,
        new_alignment_score=-np.inf, init_matrices_f=_init_matrices_nw,
        penalize_terminal_gaps=True):
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
    seq1_length = seq1.sequence_length()
    seq2_length = seq2.sequence_length()
    # cache some values for quicker/simpler access
    aend = _traceback_encoding['alignment-end']
    match = _traceback_encoding['match']
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    new_alignment_score = (new_alignment_score, aend)

    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    score_matrix, traceback_matrix = init_matrices_f(
        seq1, seq2, gap_open_penalty, gap_extend_penalty)

    # Iterate over the characters in seq2 (which corresponds to the horizontal
    # sequence in the matrix)
    for seq2_pos, seq2_chars in zip(range(1, seq2_length+1),
                                   seq2.iter_positions(str)):
        # Iterate over the characters in seq1 (which corresponds to the
        # horizontal sequence in the matrix)
        for seq1_pos, seq1_chars in zip(range(1, seq1_length+1),
                                       seq1.iter_positions(str)):
            # compute the score for a match/mismatch
            substitution_score = 0
            for seq2_char in seq2_chars:
                for seq1_char in seq1_chars:
                    if BiologicalSequence.is_gap(seq2_char) or\
                       BiologicalSequence.is_gap(seq1_char):
                           substitution_score += gap_substitution
                    else:
                        try:
                            substitution_score += \
                                substitution_matrix[seq1_char][seq2_char]
                        except KeyError:
                            offending_chars = \
                                [c for c in (seq1_char, seq2_char)
                                 if c not in substitution_matrix]
                            raise ValueError(
                                "One of the sequences contains a character that is not "
                                "contained in the substitution matrix. Are you using "
                                "an appropriate substitution matrix for your sequence "
                                "type (e.g., a nucleotide substitution matrix does not "
                                "make sense for aligning protein sequences)? Does your "
                                "sequence contain invalid characters? The offending "
                                "character(s) is: %s." % ', '.join(offending_chars))
            substitution_score /= (seq2.sequence_count() * seq1.sequence_count())

            diag_score = \
                (score_matrix[seq2_pos-1, seq1_pos-1] + substitution_score,
                 match)

            # compute the score for adding a gap in seq2 (vertical)
            if not penalize_terminal_gaps and (seq1_pos == seq1_length):
                # we've reached the end of seq1, so adding vertical gaps
                # (which become gaps in seq1) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                up_score = (score_matrix[seq2_pos-1, seq1_pos], vgap)
            elif traceback_matrix[seq2_pos-1, seq1_pos] == vgap:
                # gap extend, because the cell above was also a gap
                up_score = \
                    (score_matrix[seq2_pos-1, seq1_pos] - gap_extend_penalty,
                     vgap)
            else:
                # gap open, because the cell above was not a gap
                up_score = \
                    (score_matrix[seq2_pos-1, seq1_pos] - gap_open_penalty,
                     vgap)

            # compute the score for adding a gap in seq1 (horizontal)
            if not penalize_terminal_gaps and (seq2_pos == seq2_length):
                # we've reached the end of seq2, so adding horizontal gaps
                # (which become gaps in seq2) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                left_score = (score_matrix[seq2_pos, seq1_pos-1], hgap)
            elif traceback_matrix[seq2_pos, seq1_pos-1] == hgap:
                # gap extend, because the cell to the left was also a gap
                left_score = \
                    (score_matrix[seq2_pos, seq1_pos-1] - gap_extend_penalty,
                     hgap)
            else:
                # gap open, because the cell to the left was not a gap
                left_score = \
                    (score_matrix[seq2_pos, seq1_pos-1] - gap_open_penalty,
                     hgap)

            # identify the largest score, and use that information to populate
            # the score and traceback matrices
            best_score = _first_largest([new_alignment_score, left_score,
                                         diag_score, up_score])
            score_matrix[seq2_pos, seq1_pos] = best_score[0]
            traceback_matrix[seq2_pos, seq1_pos] = best_score[1]

    return score_matrix, traceback_matrix


def _traceback(traceback_matrix, score_matrix, aln1, aln2, start_row,
               start_col, gap_character='-'):
    # cache some values for simpler
    aend = _traceback_encoding['alignment-end']
    match = _traceback_encoding['match']
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    # initialize the result alignments
    aln1_sequence_count = aln1.sequence_count()
    aligned_seqs1 = []
    for e in range(aln1_sequence_count):
        aligned_seqs1.append([])

    aln2_sequence_count = aln2.sequence_count()
    aligned_seqs2 = []
    for e in range(aln2_sequence_count):
        aligned_seqs2.append([])

    current_row = start_row
    current_col = start_col

    best_score = score_matrix[current_row, current_col]
    current_value = None

    while current_value != aend:
        current_value = traceback_matrix[current_row, current_col]

        if current_value == match:
            for i in range(aln1_sequence_count):
                aligned_seqs1[i].append(aln1[i][current_col-1])
            for i in range(aln2_sequence_count):
                aligned_seqs2[i].append(aln2[i][current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == vgap:
            for i in range(aln1_sequence_count):
                aligned_seqs1[i].append('-')
            for i in range(aln2_sequence_count):
                aligned_seqs2[i].append(aln2[i][current_row-1])
            current_row -= 1
        elif current_value == hgap:
            for i in range(aln1_sequence_count):
                aligned_seqs1[i].append(aln1[i][current_col-1])
            for i in range(aln2_sequence_count):
                aligned_seqs2[i].append('-')
            current_col -= 1
        elif current_value == aend:
            continue
        else:
            raise ValueError(
                "Invalid value in traceback matrix: %s" % current_value)

    for i in range(aln1_sequence_count):
        aligned_seqs1[i] = BiologicalSequence(
            ''.join(map(str, aligned_seqs1[i][::-1])),
            id=_get_seq_id(aln1[0], str(i)))
    for i in range(aln2_sequence_count):
        aligned_seqs2[i] = BiologicalSequence(
            ''.join(map(str, aligned_seqs2[i][::-1])),
            id=_get_seq_id(aln2[0], str(aln1_sequence_count + i)))

    return (aligned_seqs1, aligned_seqs2, best_score,
            current_col, current_row)


def _first_largest(scores):
    """ Similar to max, but returns the first element achieving the high score

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
