#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

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
              'Y': -2, 'X': -1, 'Z': 5}
}
nt_substitution_matrix = {'A': {'A':  1, 'C': -2, 'G': -2, 'T': -2},
                          'C': {'A': -2, 'C':  1, 'G': -2, 'T': -2},
                          'G': {'A': -2, 'C': -2, 'G':  1, 'T': -2},
                          'T': {'A': -2, 'C': -2, 'G': -2, 'T':  1}}

def make_nt_substitution_matrix(match, mismatch, alphabet='ACGT'):
    result = {}
    for c1 in alphabet:
        row = {}
        for c2 in alphabet:
            if c1 == c2:
                row[c2] = match
            else:
                row[c2] = mismatch
        result[c1] = row
    return result

def local_pairwise_align_nucleotide(sequence1, sequence2, gap_open_penalty=5,
                                    gap_extend_penalty=2,
                                    substitution_matrix=None):
    """Locally align nucleotide seqs using Smith-Waterman w/ affine gap scoring

       Parameters
       ----------
       sequence1 : str or BiologicalSequence
           The first unaligned sequence
       sequence2 : str or BiologicalSequence
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the
           previous best alignment score); default is nt_substitution_matrix

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2

       Examples
       --------
       >>> from skbio import local_pairwise_align_nucleotide
       >>> s1 = "GCGTGCCTAAGGTATGCAAG"
       >>> s2 = "ACGTGCCTAGGTACGCAAG"
       >>> r = local_pairwise_align_nucleotide(s1, s2)
       >>> print(r[0])
       CGTGCCTAAGGTATGCAAG
       >>> print(r[1])
       CGTGCCTA-GGTACGCAAG

    """
    if substitution_matrix is None:
        substitution_matrix = nt_substitution_matrix

    return local_pairwise_align(sequence1, sequence2, gap_open_penalty,
                                gap_extend_penalty, substitution_matrix)


def local_pairwise_align_protein(sequence1, sequence2, gap_open_penalty=11,
                                 gap_extend_penalty=1,
                                 substitution_matrix=None):
    """Locally align two protein seqs (Smith-Waterman w affine gap scoring)

       Parameters
       ----------
       sequence1 : string
           The first unaligned sequence
       sequence2 : string
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the
           previous best alignment score)

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2

       Examples
       --------
       >>> from skbio import local_pairwise_align_protein
       >>> s1 = "HEAGAWGHEE"
       >>> s2 = "PAWHEAE"
       >>> r = local_pairwise_align_protein(s1, s2, 8, 8)
       >>> print(r[0])
       AWGHE
       >>> print(r[1])
       AW-HE

    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    return local_pairwise_align(sequence1, sequence2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix)


def local_pairwise_align(sequence1, sequence2, gap_open_penalty,
                         gap_extend_penalty, substitution_matrix):
    """Locally align sequences using Smith-Waterman w/ affine gap scoring

       Parameters
       ----------
       sequence1 : str or BiologicalSequence
           The first unaligned sequence
       sequence2 : str or BiologicalSequence
           The second unaligned sequence
       gap_open_penalty : int, float
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar)
           lookup for substitution scores (these values are added to the
           previous best alignment score); default is nt_substitution_matrix

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
       int
          The start position of the alignment in sequence 1
       int
          The start position of the alignment in sequence 2
    """
    sw_matrix, traceback_matrix = \
        _local_dynamic_programming_and_traceback(sequence1,
                                                 sequence2,
                                                 gap_open_penalty,
                                                 gap_extend_penalty,
                                                 substitution_matrix)

    return _local_traceback(traceback_matrix,sw_matrix,sequence1,sequence2)


def global_pairwise_align_nucleotide(seq1, seq2, gap_open_penalty=5,
                                    gap_extend_penalty=2,
                                    substitution_matrix=None):
    """Globally align nucleotide seqs with Needleman-Wunsch

       Parameters
       ----------
       sequence1 : str or BiologicalSequence
           The first unaligned sequence
       sequence2 : str or BiologicalSequence
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the
           previous best alignment score); default is nt_substitution_matrix

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment

    """
    if substitution_matrix is None:
        substitution_matrix = nt_substitution_matrix

    return nw_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
                     substitution_matrix)

def global_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                  gap_extend_penalty=1,
                                  substitution_matrix=None):
    """Globally align protein seqs with Needleman-Wunsch

       Parameters
       ----------
       sequence1 : str or BiologicalSequence
           The first unaligned sequence
       sequence2 : str or BiologicalSequence
           The second unaligned sequence
       gap_open_penalty : int, float, optional
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float, optional
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar), optional
           lookup for substitution scores (these values are added to the
           previous best alignment score); default is nt_substitution_matrix

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment

       Examples
       --------
       >>> from skbio import global_pairwise_align_protein
       >>> s1 = "HEAGAWGHEE"
       >>> s2 = "PAWHEAE"
       >>> r = global_pairwise_align_protein(s1, s2, 8, 8)
       >>> print(r[0])
       HEAGAWGHE-E
       >>> print(r[1])
       -PA--W-HEAE

    """
    if substitution_matrix is None:
        substitution_matrix = blosum50

    return nw_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
                     substitution_matrix)

def nw_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
             substitution_matrix):
    """Globally align seqs with Needleman-Wunsch

       Parameters
       ----------
       sequence1 : str or BiologicalSequence
           The first unaligned sequence
       sequence2 : str or BiologicalSequence
           The second unaligned sequence
       gap_open_penalty : int, float
           penalty for opening a gap (this is substracted from previous best
           alignment score, so is typically positive)
       gap_extend_penalty : int, float
           penalty for extending a gap (this is substracted from previous best
           alignment score, so is typically positive)
       substitution_matrix: 2D dict (or similar)
           lookup for substitution scores (these values are added to the
           previous best alignment score); default is nt_substitution_matrix

       Returns
       -------
       string
          The first aligned sequence
       string
          The second aligned sequence
       float
          The score of the alignment
    """
    nw_matrix, traceback_matrix = \
        _global_dynamic_programming_and_traceback(\
            seq1, seq2, gap_open_penalty, gap_extend_penalty,
            substitution_matrix)
    aligned_seq1, aligned_seq2, score = _global_traceback(\
        traceback_matrix, nw_matrix, seq1, seq2)
    return aligned_seq1, aligned_seq2, score


## Everything under here will have a lot of shared code with different
## methods, suggesting that an Aligner class is probably the way to go.

def _local_dynamic_programming_and_traceback(seq1, seq2, gap_open_penalty,
                                             gap_extend_penalty,
                                             substitution_matrix):
    """Return dynamic programming and traceback matrices

    These are generated for local alignment using Smith-Waterman.
    """
    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    sw_matrix = [[0 for i in range(0,len(seq1)+1)]]
    traceback_matrix = [[None] + [None for i in range(0,len(seq1))]]
    # Iterate over the amino acids in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in 'Biological Sequence
    # Analysis'
    for i,aa2 in zip(range(1,len(seq2)+1),seq2):
        # Initialize the current row of the matrix
        current_row = [0]
        current_traceback_matrix_row = [None]
        # Iterate over the amino acids in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in 'Biological Sequence
        # Analysis'
        new_alignment_score = (0,None)
        for j,aa1 in zip(range(1,len(seq1)+1),seq1):
            substitution_score = substitution_matrix[aa1][aa2]
            diag_score = (sw_matrix[i-1][j-1] + substitution_score,'\\')
            if traceback_matrix[i-1][j] == '|':
                # gap extend, because the cell above was also a gap
                up_score = (sw_matrix[i-1][j] - gap_extend_penalty,'|')
            else:
                # gap open, because the cell above was not a gap
                up_score = (sw_matrix[i-1][j] - gap_open_penalty,'|')
            if current_traceback_matrix_row[-1] == '-':
                # gap extend, because the cell to the left was also a gap
                left_score = (current_row[-1] - gap_extend_penalty,'-')
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (current_row[-1] - gap_open_penalty,'-')
            best_score = _first_largest([new_alignment_score, left_score,
                                         diag_score, up_score])
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        sw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return sw_matrix, traceback_matrix

def _local_traceback(traceback_matrix,sw_matrix,seq1,seq2,gap_character='-'):

    aligned_seq1 = []
    aligned_seq2 = []

    current_row = None
    current_col = None
    best_score = 0.0
    for i in range(len(sw_matrix[0])):
        for j in range(len(sw_matrix)):
            current_score = sw_matrix[j][i]
            if current_score > best_score:
                best_score = current_score
                current_row = j
                current_col = i

    while True:
        current_value = traceback_matrix[current_row][current_col]

        if current_value == '\\':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == '|':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
        elif current_value == '-':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append('-')
            current_col -= 1
        elif current_value == None:
            break
        else:
            raise ValueError(
                "Invalid value in traceback matrix: %s" % current_value)

    return (''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]),
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

def _global_dynamic_programming_and_traceback(seq1, seq2, gap_open_penalty,
    gap_extend_penalty, substitution_matrix):
    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    nw_row1 = [0.]
    for i in range(1, len(seq1)+1):
        nw_row1.append(-gap_open_penalty - ((i-1) * gap_extend_penalty))
    nw_matrix = [nw_row1]
    traceback_matrix = [[None] + ['-' for i in range(0,len(seq1))]]
    # Iterate over the amino acids in sequence two (which will correspond
    # to the vertical sequence in the matrix)
    # Note that i corresponds to column numbers, as in the 'Biological Sequence
    # Analysis' example
    for i,aa2 in zip(range(1,len(seq2)+1),seq2):
        # Initialize the current row of the matrix
        current_row = [(-gap_open_penalty - ((i-1) * gap_extend_penalty))]
        current_traceback_matrix_row = ['|']
        # Iterate over the amino acids in sequence one (which will
        # correspond to the horizontal sequence in the matrix)
        # Note that j corresponds to row numbers, as in the 'Biological Sequence
        # Analysis' example
        for j,aa1 in zip(range(1,len(seq1)+1),seq1):
            substitution_score = substitution_matrix[aa1][aa2]
            diag_score = (nw_matrix[i-1][j-1] + substitution_score,'\\')
            if traceback_matrix[i-1][j] == '|':
                # gap extend, because the cell above was also a gap
                up_score = (nw_matrix[i-1][j] - gap_extend_penalty,'|')
            else:
                # gap open, because the cell above was not a gap
                up_score = (nw_matrix[i-1][j] - gap_open_penalty,'|')
            if current_traceback_matrix_row[-1] == '-':
                # gap extend, because the cell to the left was also a gap
                left_score = (current_row[-1] - gap_extend_penalty,'-')
            else:
                # gap open, because the cell to the left was not a gap
                left_score = (current_row[-1] - gap_open_penalty,'-')
            best_score = _first_largest([left_score, diag_score, up_score])
            current_row.append(best_score[0])
            current_traceback_matrix_row.append(best_score[1])
        # append the current row to the matrix
        nw_matrix.append(current_row)
        traceback_matrix.append(current_traceback_matrix_row)
    return nw_matrix, traceback_matrix

def _global_traceback(traceback_matrix,nw_matrix,seq1,seq2,gap_character='-'):

    aligned_seq1 = []
    aligned_seq2 = []

    current_row = len(traceback_matrix) - 1
    current_col = len(traceback_matrix[0]) - 1

    best_score = nw_matrix[current_row][current_col]

    while True:
        current_value = traceback_matrix[current_row][current_col]

        if current_value == '\\':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
            current_col -= 1
        elif current_value == '|':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[current_row-1])
            current_row -= 1
        elif current_value == '-':
            aligned_seq1.append(seq1[current_col-1])
            aligned_seq2.append('-')
            current_col -= 1
        elif current_value == None:
            break
        else:
            raise ValueError, "Invalid value in traceback matrix: %s" % current_value

    return ''.join(aligned_seq1[::-1]), ''.join(aligned_seq2[::-1]), best_score
