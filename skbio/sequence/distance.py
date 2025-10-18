"""Sequence distance metrics (:mod:`skbio.sequence.distance`)
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
   :toctree:

   hamming
   kmer_distance
   jc69_correct

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

import skbio


def hamming(
    seq1,
    seq2,
    proportion=True,
    gap_policy="unique",
    degenerate_policy="unique",
):
    r"""Compute the Hamming distance between two sequences.

    The Hamming distance [1]_ between two equal-length sequences is the number of
    differing characters. It is often normalized to a proportion of the sequence
    length. This proportion is usually referred to as *p*-distance in bioinformatics.

    Parameters
    ----------
    seq1, seq2 : Sequence
        Sequences to compute the Hamming distance between.
    proportion : bool, optional
        If True (default), normalize to a proportion of the sequence length.
    gap_policy : {'unique', 'ignore'}, optional
        How to handle gaps in the sequences. "unique" (default) treats gaps as unique
        characters in the calculation. "ignore" excludes positions with either or both
        gaps from the calculation.
    degenerate_policy : {'unique', 'ignore'}, optional
        How to handle degenerate characters in the sequences. "unique" (default) treats
        them as unique characters. "ignore" excludes positions with degenerate
        characters from the calculation.

    Returns
    -------
    float
        Hamming distance between ``seq1`` and ``seq2``.

    Raises
    ------
    TypeError
        If the sequences are not ``Sequence`` instances.
    TypeError
        If the sequences are not the same type.
    ValueError
        If the sequences are not the same length.
    AttributeError
        If gap or degenerate characters are not defined for the sequences.

    See Also
    --------
    jc69_correct
    scipy.spatial.distance.hamming

    Notes
    -----
    ``np.nan`` will be returned if the sequences do not contain any characters.

    This function does not make assumptions about the sequence alphabet in use. Each
    sequence object's underlying sequence of characters are used to compute Hamming
    distance. Characters that may be considered equivalent in certain contexts (e.g.,
    `-` and `.` as gap characters) are treated as distinct characters when computing
    Hamming distance.

    References
    ----------
    .. [1] Hamming, R. W. (1950). Error detecting and error correcting codes. The Bell
       system technical journal, 29(2), 147-160.

    Examples
    --------
    >>> from skbio.sequence import Sequence
    >>> from skbio.sequence.distance import hamming
    >>> seq1 = Sequence('AGGGTA')
    >>> seq2 = Sequence('CGTTTA')
    >>> hamming(seq1, seq2)
    0.5

    """
    _check_seqs(seq1, seq2)

    L = len(seq1)

    # Hamming requires equal length sequences. We are checking this here
    # because the error you would get otherwise is cryptic.
    if len(seq2) != L:
        raise ValueError(
            "Hamming distance can only be computed between sequences of equal "
            "length (%d != %d)" % (len(seq1), len(seq2))
        )

    if L == 0:
        return np.nan

    # Create a Boolean mask of gap and/or degenerate characters.
    if gap_policy == "ignore":
        mask = seq1.gaps() | seq2.gaps()
        if degenerate_policy == "ignore":
            mask |= seq1.degenerates() | seq2.degenerates()
    elif degenerate_policy == "ignore":
        mask = seq1.degenerates() | seq2.degenerates()
    else:
        mask = None

    # Reduce sequences to valid positions, if necessary.
    if mask is not None:
        valid = ~mask
        L_valid = np.sum(valid)
        if L_valid == 0:
            return np.nan
        elif L_valid < L:
            seq1 = seq1[valid]
            seq2 = seq2[valid]
            L = L_valid

    # Count different positions
    dist = np.sum(seq1.values != seq2.values)
    if proportion:
        dist /= L

    return float(dist)


def jc69_correct(dist, chars=4):
    r"""Perform Jukes-Cantor (JC69) correction of a raw Hamming distance.

    The JC69 model [1]_ estimates the true sequence distance (number of substitutions
    per site) by correcting the observed sequence distance (*p*-distance) to account
    for repeated substitutions at the same site (a.k.a., saturation). For nucleotide
    sequences (4 characters), the JC69-corrected distance :math:`D` is calculated as:

    .. math::
        D = -\frac{3}{4} ln(1 - \frac{4}{3} p)

    Parameters
    ----------
    dist : float
        Uncorrected normalized Hamming distance (a.k.a., *p*-distance) between two
        sequences.
    chars : int, optional
        Number of definite characters in the alphabet. Default is 4, which is for
        nucleotide sequences (``A``, ``C``, ``G`` and ``T``/``U``).

    Returns
    -------
    float
        JC69-corrected distance.

    Raises
    ------
    ValueError
        If ``chars`` is less than 2.

    See Also
    --------
    hamming

    Notes
    -----
    JC69 is the most basic evolutionary model, assuming equal character frequencies and
    equal substitution rates between characters.

    The functions returns ``nan`` if ``dist >= (chars - 1) / chars``. This happens when
    the two sequences are too divergent and subsitutions are over-saturated for the
    meaningful estimation of the true evolutionary distance.

    References
    ----------
    .. [1] Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules.
       Mammalian protein metabolism, 3(21), 132.

    Examples
    --------
    >>> from skbio.sequence import DNA
    >>> from skbio.sequence.distance import hamming, jc69_correct
    >>> seq1 = DNA('AGGGTA')
    >>> seq2 = DNA('CGTTTA')
    >>> p_dist = hamming(seq1, seq2)
    >>> p_dist
    0.5

    >>> jc69_dist = jc69_correct(p_dist)
    >>> jc69_dist  # doctest: +ELLIPSIS
    0.823959...

    """
    if chars < 2:
        raise ValueError("`chars` must be at least 2.")
    frac = (chars - 1) / chars
    if dist >= frac:
        return np.nan
    elif dist == 0:
        return 0.0
    else:
        return -frac * np.log(1.0 - dist / frac)


def kmer_distance(seq1, seq2, k, overlap=True):
    r"""Compute the *k*-mer distance between a pair of sequences.

    The *k*-mer distance between two sequences is the fraction of *k*-mer that are
    unique to either sequence.

    Parameters
    ----------
    seq1, seq2 : Sequence
        Sequences to compute *k*-mer distance between.
    k : int
        The *k*-mer length.
    overlap : bool, optional
        Defines whether the *k*-mers should be overlapping or not.

    Returns
    -------
    float
        *k*-mer distance between ``seq1`` and ``seq2``.

    Raises
    ------
    ValueError
        If ``k`` is less than 1.
    TypeError
        If the sequences are not ``Sequence`` instances.
    TypeError
        If the sequences are not the same type.

    Notes
    -----
    *k*-mers counts are not incorporated in this distance metric.

    ``np.nan`` will be returned if there are no kmers defined for the
    sequences.

    Examples
    --------
    >>> from skbio.sequence import Sequence
    >>> seq1 = Sequence('ATCGGCGAT')
    >>> seq2 = Sequence('GCAGATGTG')
    >>> kmer_distance(seq1, seq2, 3)  # doctest: +ELLIPSIS
    0.9230769230...

    """
    _check_seqs(seq1, seq2)
    seq1_kmers = set(map(str, seq1.iter_kmers(k, overlap=overlap)))
    seq2_kmers = set(map(str, seq2.iter_kmers(k, overlap=overlap)))
    all_kmers = seq1_kmers | seq2_kmers
    if not all_kmers:
        return np.nan
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
                % type(seq).__name__
            )

    # Asserts sequences have the same type
    if type(seq1) is not type(seq2):
        raise TypeError(
            "Sequences must have matching type. Type %r does not match type %r"
            % (type(seq1).__name__, type(seq2).__name__)
        )
