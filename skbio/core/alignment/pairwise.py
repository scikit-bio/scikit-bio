#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
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
                                    match_score=2, mismatch_score=-3):
    """Locally align exactly two nucleotide seqs with Smith-Waterman

       Parameters
       ----------
       seq1 : str or BiologicalSequence
           The first unaligned sequence.
       seq2 : str or BiologicalSequence
           The second unaligned sequence.
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
           Penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive).
       match_score : int, float
           The score to add for a match between a pair of bases (this is added
           to the previous best alignment score, so is typically positive).
       mismatch_score : int, float
           The score to add for a mismatch between a pair of bases (this is
           added to the previous best alignment score, so is typically
           negative).

       Returns
       -------
       skbio.Alignment
           ``Alignment`` object containing the aligned sequences as well as
           details about the alignment.

       Notes
       -----
       Default ``match_score``, ``mismatch_score``, ``gap_open_penalty`` and
       ``gap_extend_penalty`` parameters are derived from the NCBI BLAST
       Server.

    """
    substitution_matrix = \
        _make_nt_substitution_matrix(match_score, mismatch_score)

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
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
           Penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive).
       substitution_matrix: 2D dict (or similar)
           Lookup for substitution scores (these values are added to the
           previous best alignment score); default is BLOSUM 50.

       Returns
       -------
       skbio.Alignment
           ``Alignment`` object containing the aligned sequences as well as
           details about the alignment.

       Notes
       -----
       Default ``gap_open_penalty`` and ``gap_extend_penalty`` parameters are
       derived from the NCBI BLAST Server.

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
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
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
    result = Alignment(
        [BiologicalSequence(aligned1, id=0),
         BiologicalSequence(aligned2, id=1)],
        score=score, start_end_positions=start_end_positions)
    return result


def global_pairwise_align_nucleotide(seq1, seq2, gap_open_penalty=5,
                                     gap_extend_penalty=2,
                                     match_score=1, mismatch_score=-2):
    """Globally align exactly two nucleotide seqs with Needleman-Wunsch

       Parameters
       ----------
       seq1 : str or BiologicalSequence
           The first unaligned sequence.
       seq2 : str or BiologicalSequence
           The second unaligned sequence.
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
           Penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive).
       match_score : int, float
           The score to add for a match between a pair of bases (this is added
           to the previous best alignment score, so is typically positive).
       mismatch_score : int, float
           The score to add for a mismatch between a pair of bases (this is
           added to the previous best alignment score, so is typically
           negative).

       Returns
       -------
       skbio.Alignment
           ``Alignment`` object containing the aligned sequences as well as
           details about the alignment.

       Notes
       -----
       Default ``match_score``, ``mismatch_score``, ``gap_open_penalty`` and
       ``gap_extend_penalty`` parameters are derived from the NCBI BLAST
       Server.

    """
    substitution_matrix = \
        _make_nt_substitution_matrix(match_score, mismatch_score)

    return global_pairwise_align(seq1, seq2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix)


def global_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                  gap_extend_penalty=1,
                                  substitution_matrix=None):
    """Globally align exactly two protein seqs with Needleman-Wunsch

       Parameters
       ----------
       seq1 : str or BiologicalSequence
           The first unaligned sequence.
       seq2 : str or BiologicalSequence
           The second unaligned sequence.
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
           Penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive).
       substitution_matrix: 2D dict (or similar)
           Lookup for substitution scores (these values are added to the
           previous best alignment score); default is BLOSUM 50.

       Returns
       -------
       skbio.Alignment
           ``Alignment`` object containing the aligned sequences as well as
           details about the alignment.

       Notes
       -----
       Default ``gap_open_penalty`` and ``gap_extend_penalty`` parameters are
       derived from the NCBI BLAST Server.

    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    return global_pairwise_align(seq1, seq2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix)


def global_pairwise_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
                          substitution_matrix):
    """Globally align exactly two seqs with Needleman-Wunsch

       Parameters
       ----------
       seq1 : str or BiologicalSequence
           The first unaligned sequence.
       seq2 : str or BiologicalSequence
           The second unaligned sequence.
       gap_open_penalty : int, float
           Penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive).
       gap_extend_penalty : int, float
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

    score_matrix, traceback_matrix = \
        _compute_score_and_traceback_matrices(
            seq1, seq2, gap_open_penalty, gap_extend_penalty,
            substitution_matrix, new_alignment_score=-np.inf,
            init_matrices_f=_init_matrices_nw)

    end_row_position = traceback_matrix.shape[0] - 1
    end_col_position = traceback_matrix.shape[1] - 1

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = \
        _traceback(traceback_matrix, score_matrix, seq1, seq2,
                   end_row_position, end_col_position)
    start_end_positions = [(seq1_start_position, end_col_position-1),
                           (seq2_start_position, end_row_position-1)]
    result = Alignment(
        [BiologicalSequence(aligned1, id=0),
         BiologicalSequence(aligned2, id=1)],
        score=score, start_end_positions=start_end_positions)
    return result

# Functions from here allow for generalized (global or local) alignment. I
# will likely want to put these in a single object to make the naming a little
# less clunky.

_traceback_encoding = {'match': 1, 'vertical-gap': 2, 'horizontal-gap': 3,
                       'uninitialized': -1, 'alignment-end': 0}


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
    shape = (len(seq2)+1, len(seq1)+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int16)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, 0] = _traceback_encoding['alignment-end']
    return score_matrix, traceback_matrix


def _init_matrices_nw(seq1, seq2, gap_open_penalty, gap_extend_penalty):
    shape = (len(seq2)+1, len(seq1)+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int16)
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


def _compute_score_and_traceback_matrices(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix,
        new_alignment_score=-np.inf, init_matrices_f=_init_matrices_nw):
    """Return dynamic programming (score) and traceback matrices
    """
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
    # Iterate over the characters in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in 'Biological Sequence
    # Analysis'
    for i, c2 in zip(range(1, len(seq2)+1), seq2):
        # Iterate over the characters in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in 'Biological Sequence
        # Analysis'
        for j, c1 in zip(range(1, len(seq1)+1), seq1):
            substitution_score = substitution_matrix[c1][c2]
            diag_score = (score_matrix[i-1][j-1] + substitution_score, match)
            if traceback_matrix[i-1][j] == vgap:
                # gap extend, because the cell above was also a gap
                up_score = (score_matrix[i-1][j] - gap_extend_penalty, vgap)
            else:
                # gap open, because the cell above was not a gap
                up_score = (score_matrix[i-1][j] - gap_open_penalty, vgap)
            if traceback_matrix[i][j-1] == hgap:
                # gap extend, because the cell to the left was also a gap
                left_score = (score_matrix[i][j-1] - gap_extend_penalty, hgap)
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (score_matrix[i][j-1] - gap_open_penalty, hgap)
            best_score = _first_largest([new_alignment_score, left_score,
                                         diag_score, up_score])
            score_matrix[i][j] = best_score[0]
            traceback_matrix[i][j] = best_score[1]
    return score_matrix, traceback_matrix


def _traceback(traceback_matrix, score_matrix, seq1, seq2, start_row,
               start_col, gap_character='-'):
    # cache some values for simpler
    aend = _traceback_encoding['alignment-end']
    match = _traceback_encoding['match']
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    aligned_seq1 = []
    aligned_seq2 = []

    current_row = start_row
    current_col = start_col

    best_score = score_matrix[current_row][current_col]

    while True:
        current_value = traceback_matrix[current_row][current_col]

        if current_value == match:
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == vgap:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
        elif current_value == hgap:
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append('-')
            current_col -= 1
        elif current_value == aend:
            break
        else:
            raise ValueError(
                "Invalid value in traceback matrix: %s" % current_value)

    return (''.join(map(str, aligned_seq1[::-1])),
            ''.join(map(str, aligned_seq2[::-1])),
            best_score, current_col, current_row)


def _first_largest(scores):
    """ Similar to max, but returns the first element achieving the high score

        If max receives a tuple, it will break a score for the highest value
        of entry[i] with entry[i+1]. We don't want that here - to better match
        with the results of other tools, it's better to be able to set the
        order of the choices.
    """
    result = scores[0]
    for score, direction in scores[1:]:
        if score > result[0]:
            result = (score, direction)
    return result
