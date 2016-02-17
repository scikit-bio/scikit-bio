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
    for seq in seq1, seq2:
        if not isinstance(seq, skbio.Sequence):
            raise TypeError(
                "`seq1` and `seq2` must be Sequence instances, not %r"
                % type(seq).__name__)

    if type(seq1) is not type(seq2):
        raise TypeError(
            "Sequences must have matching type. Type %r does not match type %r"
            % (type(seq1).__name__, type(seq2).__name__))

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
