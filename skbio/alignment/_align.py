# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.alignment import TabularMSA, AlignPath
from skbio.sequence import DNA
from ._c_align import _align_score


def pair_score(
    seq1=None,
    seq2=None,
    path=None,
    sub_score=(2, -3),
    gap_score=(-5, -2),
    term_score=False,
):
    """Calculate score of a pair of aligned sequences.

    Parameters
    ----------
    seq1 : Sequence or str, optional
        The first aligned sequence.

    seq2 : Sequence or str, optional
        The second aligned sequence.

    path : PairAlignPath, optional
        Paiwise alignment path. If provided, `seq1` and `seq2` are considered as
        unaligned. Can be provided alone without `seq1` and `seq2` if `sub_score` is a
        single number.

    sub_score : int, float, array_like of (2,), SubstitutionMatrix, or str, optional
        Score of a match, mismatch or substitution. It can be one of the following:
        - Single number: A uniform score for match and mismatch.
        - Tuple of two numbers: Match score and mismatch score.
        - SubstitutionMatrix: A matrix of substitution scores.
        - String: Name of the substitution matrix that can be recognized by
          `SubstitutionMatrix.by_name`.
        Default is (2, -3).

    gap_score : int, float, or array_like of (2,), optional
        Score of a gap. The value is usually negative, indicating a penalty to the
        alignment score. It can be one of the following:
        - One number: Linear gap penalty. Each gap character is penalized by this value
          (x). A continuous gap of length n has a total penalty of x * n.
        - Two numbers: Affine gap penalty. The two numbers (a, b) represent gap opening
          penalty and gap extension penalty. A continuous gap of length n has a total
          penalty of a + b * n.
        Default is (-5, -2).

    term_score : bool, optional
        Whether to score terminal gaps.

    Returns
    -------
    float
        Pairwise alignment score.

    """
    pass


def align_score(seq1, seq2, subMatrix, gap_open, gap_extend):
    """
    Computes the alignment score between two sequences using the given substitution
    matrix and gap penalties.

    Parameters:
    ----------
    seq1 : string or any scikit-bio DNA compatible class
        The first sequence, represented as a string, DNA, or any other DNA-compatible
        class.
    seq2 : string or any scikit-bio DNA compatible class
        The second sequence, represented as a string, DNA, or any other DNA-compatible
        class.
    subMatrix : skbio.sequence.SubstitutionMatrix
        Substitution matrix, taken from the scikit-bio sequence module.
    gap_open : float
        The penalty for opening a gap.
    gap_extend : float
        The penalty for extending an existing gap.

    Returns:
    -------
    float
        The computed alignment score.
    """
    # Has sequences through substitution matrix
    seq1_idx = subMatrix._char_hash[DNA(seq1)._bytes]
    seq2_idx = subMatrix._char_hash[DNA(seq2)._bytes]
    # Run though helper function to get ouput
    return _align_score(
        seq1_idx, seq2_idx, subMatrix._data.astype(np.float32), gap_open, gap_extend
    )
