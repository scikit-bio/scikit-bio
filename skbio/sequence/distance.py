"""
Sequence distance metrics (:mod:`skbio.sequence.distance`)
==========================================================

.. currentmodule:: skbio.sequence.distance

This module contains functions for computing distances between scikit-bio
``Sequence`` objects. These functions can be used directly or supplied to other
parts of the scikit-bio API that accept a sequence distance metric as input,
such as :meth:`skbio.sequence.Sequence.distance` and
:meth:`skbio.stats.distance.DistanceMatrix.from_iterable`.

Functions
---------

.. autosummary::
   :toctree: generated/

   hamming

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np
import scipy.spatial.distance

import skbio
from skbio.util._decorator import experimental


@experimental(as_of='0.4.2')
def hamming(seq1, seq2):
    """Compute Hamming distance between two sequences.

    The Hamming distance between two equal-length sequences is the proportion
    of differing characters.

    Parameters
    ----------
    seq1, seq2 : Sequence
        Sequences to compute Hamming distance between.

    Returns
    -------
    float
        Hamming distance between `seq1` and `seq2`.

    Raises
    ------
    TypeError
        If `seq1` and `seq2` are not ``Sequence`` instances.
    TypeError
        If `seq1` and `seq2` are not the same type.
    ValueError
        If `seq1` and `seq2` are not the same length.

    See Also
    --------
    scipy.spatial.distance.hamming

    Notes
    -----
    ``np.nan`` will be returned if the sequences do not contain any characters.

    This function does not make assumptions about the sequence alphabet in use.
    Each sequence object's underlying sequence of characters are used to
    compute Hamming distance. Characters that may be considered equivalent in
    certain contexts (e.g., `-` and `.` as gap characters) are treated as
    distinct characters when computing Hamming distance.

    Examples
    --------
    >>> from skbio import Sequence
    >>> from skbio.sequence.distance import hamming
    >>> seq1 = Sequence('AGGGTA')
    >>> seq2 = Sequence('CGTTTA')
    >>> hamming(seq1, seq2)
    0.5

    """
    _check_seqs(seq1, seq2)

    # Hamming requires equal length sequences. We are checking this here
    # because the error you would get otherwise is cryptic.
    if len(seq1) != len(seq2):
        raise ValueError(
            "Hamming distance can only be computed between sequences of equal "
            "length (%d != %d)" % (len(seq1), len(seq2)))

    # scipy throws a RuntimeWarning when computing Hamming distance on length 0
    # input.
    if not seq1:
        distance = np.nan
    else:
        distance = scipy.spatial.distance.hamming(seq1.values, seq2.values)

    return float(distance)


@experimental(as_of='0.4.2')
def kmer_distance(seq1, seq2, k=3, overlap=True):
    """Compute the kmer distance between a pair of sequences
    Parameters
    ----------
    seq1 : skbio.Sequence
    seq2 : skbio.Sequence
    k : int, optional
        The word length.
    overlap : bool, optional
        Defines whether the k-words should be overlapping or not
        overlapping.
    Returns
    -------
    float
        Fraction of the set of k-mers from both seq1 and
        seq2 that are unique to either seq1 or
        seq2.
    Raises
    ------
    ValueError
        If k < 1.
    Notes
    -----
    k-mer counts are not incorporated in this distance metric.
    """
    _check_seqs(seq1, seq2)
    seq1_kmers = set(map(str, seq1.iter_kmers(k, overlap)))
    seq2_kmers = set(map(str, seq2.iter_kmers(k, overlap)))
    all_kmers = seq1_kmers | seq2_kmers
    shared_kmers = seq1_kmers & seq2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique


def _check_seqs(seq1, seq2):
    # Asserts both sequences are skbio.sequence objects
    for seq in seq1, seq2:
        if not isinstance(seq, skbio.Sequence):
            raise TypeError(
                "`seq1` and `seq2` must be Sequence instances, not %r"
                % type(seq).__name__)

    # Asserts sequences have the same type
    if type(seq1) is not type(seq2):
        raise TypeError(
            "Sequences must have matching type. Type %r does not match type %r"
            % (type(seq1).__name__, type(seq2).__name__))
