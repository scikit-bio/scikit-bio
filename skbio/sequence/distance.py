"""Sequence distance metrics (:mod:`skbio.sequence.distance`)
==========================================================

.. currentmodule:: skbio.sequence.distance

This module contains functions for computing distances between scikit-bio
``Sequence`` objects. These functions can be used directly or supplied to other
parts of the scikit-bio API that accept a sequence distance metric as input,
such as :meth:`skbio.sequence.Sequence.distance` and
:meth:`skbio.stats.distance.DistanceMatrix.from_iterable`.

Generic distance metrics
------------------------

.. autosummary::
   :toctree:

   hamming
   kmer_distance


Distance metrics for nucleotide sequences
-----------------------------------------

   jc69


Utility functions
-----------------

   jc69_correct

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, Any, TYPE_CHECKING

import functools

import numpy as np
from scipy.spatial.distance import pdist

from skbio.sequence import Sequence, GrammaredSequence, DNA, RNA, Protein
from skbio.sequence._alphabet import _encode_alphabet

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import ArrayLike

# ----------------------------------------------------------------------------
# Functions of this module are organized in the following structure: Each
# distance metric "xyz" has three functions:
#
# - `xyz` is a public-facing function that takes a pair of sequences as input
#   and outputs a single number.
# - `_xyz_full` is a private function that takes a 2-D array representing the
#   ASCII codes of multiple aligned sequences as input, and generates a 1-D
#   array representing a condensed distance matrix between the sequences.
# - `_xyz_pair` resembles `_xyz_full` but additionally consumes a Boolean mask
#   of the ASCII code array representing valid sites.
#
# The two private functions are called by `skbio.alignment.align_dists` and
# this is more efficient than calling the public function between each pair of
# sequences.
# ----------------------------------------------------------------------------


def _metric_specs(
    seqtype: Optional[type] = None,
    equal: bool = False,
    alphabet: Optional[Union[str, "ArrayLike"]] = None,
):
    r"""Specifications of a sequence distance metric.

    Parameters
    ----------
    func : callable
        Function that calculates a sequence distance metric.

    seqtype : type or tuple of types, optional
        Valid sequence types. Can be a single type (such as ``Protein``) or a tuple of
        multiple types (such as ``(DNA, RNA)``). If None (default), any ``Sequence``
        objects, grammared or not, are valid.

    equal : bool, optional
        If True, the two sequences must have the same length. Default is False.

    alphabet : str, 1D array_like, {'nongap', 'definite', 'canonical'}, optional
        An alphabet of valid characters to be considered by the metric. Can be a string
        or array-like of characters or their ASCII codes. Three special keywords are
        recognized: "nongap" (excluding gap characters such as "-" and "*"), "definite"
        (definite characters, i.e., non-degenerate and non-gap characters), and
        "canonical" (canonical characters, such as the four nucleobases and the 20
        basic amino acids). Keywords can be used only when sequences are grammared.
        If None (default), all characters are valid.

        Invalid characters will be removed from the sequences prior to calculation. If
        ``equal=True`` is set, positions with either or both invalid characters in the
        two sequences will be removed.

    Returns
    -------
    callable
        Decorated function.

    Notes
    -----
    This function serves as a decorator for individual functions that calculate
    sequence distance metrics. The two positional arguments of a decorated
    function must be two ``Sequence`` objects. Additional arguments may follow.

    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(seq1, seq2, *args, **kwargs):
            # check if sequences are skbio.sequence objects
            for seq in seq1, seq2:
                if not isinstance(seq, Sequence):
                    raise TypeError(
                        "Sequences must be skbio.sequence.Sequence instances, not "
                        f"{type(seq).__name__!r}."
                    )

            # check if sequences have the same type
            if type(seq1) is not type(seq2):
                raise TypeError(
                    f"Sequences must have matching type. {type(seq1).__name__!r} "
                    f"does not match {type(seq2).__name__!r}."
                )

            # check if sequences have the expected type
            if seqtype is not None:
                for seq in seq1, seq2:
                    if not isinstance(seq, seqtype):
                        if isinstance(seqtype, tuple):
                            names_ = tuple(x.__name__ for x in seqtype)
                        else:
                            names_ = seqtype.__name__
                        raise TypeError(
                            f"Sequences must be {names_!r} instances, not "
                            f"{type(seq).__name__!r}."
                        )

            # check if sequences have the same length
            if equal and len(seq1) != len(seq2):
                raise ValueError(
                    f"{func.__name__!r} can only be calculated between equal-length "
                    f"sequences. {len(seq1)} != {len(seq2)}."
                )

            # filter sequences by a given alphabet
            if alphabet is not None:
                if alphabet in ("nongap", "definite", "canonical"):
                    valid = getattr(seq1, f"_{alphabet}_hash")
                else:
                    encoded = _encode_alphabet(alphabet)
                    valid = np.zeros((Sequence._num_ascii_codes,), dtype=bool)
                    valid[encoded] = True

                if equal:
                    pos = valid[seq1._bytes] & valid[seq2._bytes]
                    seq1, seq2 = seq1[pos], seq2[pos]
                else:
                    seq1, seq2 = seq1[valid[seq1._bytes]], seq2[valid[seq2._bytes]]

            # call function to calculate sequence distance
            return func(seq1, seq2, *args, **kwargs)

        # transfer parameters to function attributes
        wrapper._is_metric = True
        wrapper._seqtype = seqtype
        wrapper._equal = equal
        wrapper._alphabet = alphabet

        return wrapper

    return decorator


@_metric_specs(equal=True)
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
    if (L := len(seq1)) == 0:
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


def _get_valid(seq1, seq2):
    L = len(seq1)
    if L == 0:
        return 0, None, None
    valid = seq1.definites() | seq2.definites()
    L_valid = np.sum(valid)
    if L_valid == 0:
        return 0, None, None
    elif L_valid < L:
        return L_valid, seq1[valid], seq2[valid]
    else:
        return L, seq1, seq2


@_metric_specs(equal=True, seqtype=GrammaredSequence)
def p_dist(seq1, seq2):
    """Calculate p-distance between two aligned sequences."""
    L, seq1, seq2 = _get_valid(seq1, seq2)
    return np.count_nonzero(seq1 != seq2) / L


def _p_dist(seqs):
    """Compute pairwise p-distances between multiple sequences."""
    return pdist(seqs, metric="hamming")


def _p_dist_pair(seqs, mask):
    n = seqs.shape[0]
    n_1 = n - 1
    dm = np.empty((n * n_1 // 2,))
    start = 0
    for i in range(n_1):
        end = start + n_1 - i
        target = dm[start:end]

        mask_ = mask[i] & mask[i + 1 :]  ##
        diff = (seqs[i] != seqs[i + 1 :]) & mask_  ##, ##
        p = np.count_nonzero(diff, axis=1)  ## faster than np.sum, no out=.

        L = np.sum(mask_, axis=1)
        np.divide(p, L, out=target, where=L > 0)

        start = end
    return dm


@_metric_specs(equal=True, seqtype=(DNA, RNA))
def jc69(seq1, seq2):
    """Compute the JC69 distance between two sequences."""
    return jc69_correct(p_dist(seq1, seq2))


def _jc69(seqs, chars=4):
    """Compute pairwise JC69 distances between sequences."""
    return jc69_correct(_p_dist(seqs), chars)


def _jc69_pair(seqs, mask, chars=4):
    """Compute pairwise JC69 distances between sequences."""
    return jc69_correct(_p_dist_pair(seqs, mask), chars)


def jc69_correct(dists, chars=4):
    r"""Perform Jukes-Cantor (JC69) correction of raw Hamming distances.

    The JC69 model [1]_ estimates the true evolutionary distance (number of
    substitutions per site) between two sequences by correcting the observed sequence
    distance (*p*-distance) to account for repeated substitutions at the same site
    (i.e., saturation). For nucleotide sequences (4 characters), the JC69-corrected
    distance :math:`D` is calculated as:

    .. math::
        D = -\frac{3}{4} ln(1 - \frac{4}{3} p)

    Parameters
    ----------
    dists : float or array_like
        Uncorrected normalized Hamming distances (a.k.a., *p*-distance) between pairs
        of sequences.
    chars : int, optional
        Number of definite characters in the alphabet. Default is 4, which is for
        nucleotide sequences ("A", "C", "G" and "T"/"U").

    Returns
    -------
    float or ndarray
        JC69-corrected distances.

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
    estimation of the true evolutionary distance.

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
    is_scalar = np.isscalar(dists)
    arr = np.asarray(dists)

    # copy of array?
    mask = arr < frac
    vals = arr[mask]
    vals /= -frac
    vals += 1.0
    np.log(vals, out=vals)
    vals *= -frac
    res = np.full(arr.shape, np.nan)
    res[mask] = vals

    if is_scalar:
        return res.item()
    return res


@_metric_specs(equal=True)
def transitions(seq1, seq2, proportion=True):
    L, seq1, seq2 = _get_valid(seq1, seq2)


@_metric_specs()
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
    seq1_kmers = set(map(str, seq1.iter_kmers(k, overlap=overlap)))
    seq2_kmers = set(map(str, seq2.iter_kmers(k, overlap=overlap)))
    all_kmers = seq1_kmers | seq2_kmers
    if not all_kmers:
        return np.nan
    shared_kmers = seq1_kmers & seq2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique
