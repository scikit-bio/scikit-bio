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
   pdist
   logdet
   kmer_distance


Nucleotide distance metrics
---------------------------

.. autosummary::
   :toctree:

   jc69
   f81
   k2p
   f84
   tn93


Utility functions
-----------------

.. autosummary::
   :toctree:

   jc69_correct

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, TYPE_CHECKING
from inspect import signature
import functools

import numpy as np

from skbio.sequence import Sequence, GrammaredSequence, DNA, RNA, Protein
from skbio.sequence._alphabet import _encode_alphabet

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import ArrayLike

# ----------------------------------------------------------------------------
# Functions of this module are organized in the following structure: Each
# distance metric "xyz" has two functions:
#
# - `xyz` is a public-facing function that takes a pair of sequences as input
#   and outputs a single number. It has several attributes specifying its
#   behavior. See `_metric_specs` for details.
#
# - `_xyz` is a private function that consumes a 2-D array representing the
#   ASCII codes of multiple aligned sequences, and generates a condensed
#   distance matrix between them. It is typically called by `skbio.alignment.
#   align_dists`. Optionally, it consumes a Boolean mask of the input array,
#   representing valid sites.
#
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
        name = func.__name__
        has_freqs = "freqs" in signature(func).parameters

        @functools.wraps(func)
        def wrapper(seq1, seq2, *args, **kwargs):
            # Check if sequences are skbio.sequence objects.
            type1 = type(seq1)
            if not issubclass(type1, Sequence):
                raise TypeError(
                    "Sequences must be skbio.sequence.Sequence instances, not "
                    f"{type1.__name__!r}."
                )

            # Check if sequences have the same type.
            if type1 is not type(seq2):
                raise TypeError(
                    f"Sequences must have matching type. {type1.__name__!r} "
                    f"does not match {type(seq2).__name__!r}."
                )

            # Check if sequences have the expected type.
            _check_seqtype(name, type1, seqtype)

            # Check if sequences have the same length. This is for metrics that take
            # aligned sequences.
            if equal:
                if (L1 := len(seq1)) != len(seq2):
                    raise ValueError(
                        f"{name!r} can only be calculated between equal-length "
                        f"sequences. {L1} != {len(seq2)}."
                    )

                # Return a NaN if sequences are empty.
                if L1 == 0:
                    return float("nan")

            # Get valid characters.
            valid = _char_hash(alphabet, type1)

            # Calculate character frequencies if needed.
            if has_freqs:
                if (freqs := kwargs.get("freqs")) is None:
                    kwargs["freqs"] = _char_freqs((seq1._bytes, seq2._bytes), valid)
                else:
                    kwargs["freqs"] = _check_freqs(freqs)

            # Filter sequences to valid characters.
            if valid is not None:
                if equal:
                    pos = valid[seq1._bytes] & valid[seq2._bytes]
                    seq1, seq2 = seq1[pos], seq2[pos]
                    if len(seq1) == 0:
                        return float("nan")
                else:
                    seq1, seq2 = seq1[valid[seq1._bytes]], seq2[valid[seq2._bytes]]

            # Call function to calculate sequence distance.
            return func(seq1, seq2, *args, **kwargs)

        # Transfer parameters to function attributes.
        wrapper._is_metric = True
        wrapper._seqtype = seqtype
        wrapper._equal = equal
        wrapper._alphabet = alphabet
        wrapper._has_freqs = has_freqs

        return wrapper

    return decorator


def _check_seqtype(name, this, valid=None):
    """Check if the input sequences have the expected sequence type.

    Parameters
    ----------
    name : str
        Name of the metric.
    this : type
        Type of the current sequence.
    valid : type or tuple of type, optional
        Valid sequence type(s).

    Raises
    ------
    TypeError
        If sequence type is incompatible with the given metric.

    """
    if valid is None or issubclass(this, valid):
        return
    if isinstance(valid, tuple):
        types = tuple(x.__name__ for x in valid)
    else:
        types = valid.__name__
    raise TypeError(
        f"{name!r} is compatible with {types!r} sequences, not {this.__name__!r}."
    )


def _char_hash(alphabet, seqtype):
    """Get a hash table of valid characters for sequence filtering.

    Parameters
    ----------
    alphabet : str, 1D array_like, {'nongap', 'definite', 'canonical'}, or None
        An alphabet of valid characters to be considered by the metric.
    seqtype : type
        Sequence type, which stores grammar information.

    Returns
    -------
    ndarray of bool of shape (128,)
        Boolean mask of ASCII codes of valid characters (True).

    """
    if alphabet is None:
        return
    elif alphabet in ("nongap", "definite", "canonical"):
        return getattr(seqtype, f"_{alphabet}_hash")
    else:
        encoded = _encode_alphabet(alphabet)
        valid = np.zeros((seqtype._num_ascii_codes,), dtype=bool)
        valid[encoded] = True
        return valid


# Mappings of canonical characters to indices
_char_indices = {}


def _char_index(seqtype):
    """Get a mapping of characters to indices.

    It generates a vector in which indices are ASCII codes and values are indices
    within the alphabet. Only canonical characters are included in the indices.
    Other characters are masked with k (k is the number of canonical characters).

    For example, DNA will have: A: 0, C: 1, G: 2, T: 3, others: 4.

    Parameters
    ----------
    seqtype : type
        Sequence type. Must be a subclass of `GrammaredSequence`, which defines
        canonical
        grammar information.

    Returns
    -------
    nchar : int
        Number of canonical characters.
    index : ndarray of int of shape (128,)
        Mapping of characters (ASCII codes) to indices.

    """
    if seqtype not in _char_indices:
        chars = seqtype._canonical_codes
        nchar = chars.size
        index = np.full(seqtype._num_ascii_codes, nchar, dtype=int)
        index[chars] = np.arange(nchar)
        _char_indices[seqtype] = (nchar, index)
    return _char_indices[seqtype]


def _char_freqs(seqs, valid=None):
    """Calculate relative frequencies of characters in sequences.

    Parameters
    ----------
    seqs : array_like of int
        Input sequences as ASCII codes (0-127).
    valid : ndarray of bool of shape (128,), optional
        Boolean mask of valid ASCII codes (True). If omitted, all ASCII codes (n = 128)
        will be included in the result.

    Returns
    -------
    ndarray of float of shape (n_alphabet,)
        Relative frequencies of characters.

    """
    # This is an efficient way to count character frequencies. `ravel` makes a memory
    # view of the data without copying, unless needed. `bincount` counts frequencies of
    # all ASCII codes. Then only frequencies of valid characters will be returned.
    freqs = np.bincount(np.ravel(seqs), minlength=128)
    if valid is not None:
        chars = np.flatnonzero(valid)
        freqs = freqs[chars]
    total = np.sum(freqs)
    if total > 0:
        return freqs / total
    else:
        return np.full(freqs.size, np.nan)


def _check_freqs(freqs, nonzero=False):
    """Validate character frequencies."""
    freqs = np.asarray(freqs)

    if nonzero:
        if not np.all(freqs > 0.0):
            raise ValueError("Character frequencies must all be positive.")
    else:
        if not np.all(freqs >= 0.0):
            raise ValueError("Character frequencies must all be non-negative.")

    if not np.isclose(np.sum(freqs), 1.0):
        raise ValueError("Character frequencies must sum to 1.0.")

    return freqs


def _build_dm(func, seqs, mask):
    n = seqs.shape[0]
    n_1 = n - 1
    dm = np.empty((n * n_1 // 2,))
    # if L == 0:
    #     dm.fill(np.nan)
    #     return dm
    # masked = mask is not None
    start = 0
    for i in range(n_1):
        end = start + n_1 - i
        # target = dm[start:end]

        func(seqs, mask, i, dm[start:end])

        # diff = (seqs[i] != seqs[i + 1 :])

        # if masked:
        #     sites = mask[i] & mask[i + 1 :]
        #     L = np.count_nonzero(sites, axis=1)

        #     diff &= sites

        #     # non-empty sequences
        #     filled = L > 0
        #     all_valid = np.all(filled):
        #     if not all_valid:
        #         target[~filled] = np.nan

        # else:
        #     all_valid = True

        # # identical sequences
        # p = np.count_nonzero(diff, axis=1)
        # ident = p == 0

        # if all_valid:
        #     target[:] = func(p, L)
        # else:
        #     target[valid] = func(p, L)

        start = end
    return dm


@_metric_specs(equal=True)
def hamming(seq1, seq2, proportion=True):
    r"""Compute the Hamming distance between two sequences.

    The Hamming distance [1]_ between two equal-length sequences is the number of
    differing characters. It is often normalized to a proportion of the sequence
    length.

    Parameters
    ----------
    seq1, seq2 : Sequence
        Sequences to compute the Hamming distance between.
    proportion : bool, optional
        If True (default), normalize to a proportion of the sequence length.

        .. versionadded:: 0.7.2

    Returns
    -------
    float
        Hamming distance between the two sequences.

    Raises
    ------
    TypeError
        If the sequences are not ``Sequence`` instances.
    TypeError
        If the sequences are not the same type.
    ValueError
        If the sequences are not the same length.

    See Also
    --------
    pdist
    scipy.spatial.distance.hamming

    Notes
    -----
    This function does not make assumptions about the sequence alphabet in use. All
    characters of each sequence, including gaps and ambiguous codes, are used to
    compute Hamming distance. Characters that may be considered equivalent in certain
    contexts (e.g., "-" and "." as gap characters) are treated as distinct characters
    when computing Hamming distance.

    If this behavior is not desired, consider using :func:`pdist` instead.

    NaN will be returned if the sequences do not contain any characters.

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
    return _hamming(
        np.vstack((seq1._bytes, seq2._bytes)), None, None, proportion=proportion
    ).item()


def _hamming(seqs, mask, seqtype, proportion=True):
    """Compute pairwise Hamming distances between sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Sequences to compute the Hamming distance between.
    mask : ndarray of bool of shape (n_sequences, n_positions)
        A placeholder. Hamming distance always consider all characters.
    seqtype : type
        A placeholder. Hamming distance is irrelevant to sequence type.
    proportion : bool, optional
        If True (default), normalize to a proportion of the sequence length.

    Returns
    -------
    ndarray of float of shape (C(n_sequences, 2),)
        Hamming distance matrix in condensed form.

    """
    # When `proportion=True`, the result is equivalent to SciPy's `pdist(seqs,
    # metric="hamming")`. But it seems that the following code is faster.
    npos = seqs.shape[1]

    def func(seqs, mask, i, out):
        subs = seqs[i] != seqs[i + 1 :]
        p = np.count_nonzero(subs, axis=1)

        # normalize by sequence length
        if proportion:
            np.divide(p, npos, out=out)
        else:
            out[:] = p

    return _build_dm(func, seqs, mask)


@_metric_specs(equal=True, seqtype=GrammaredSequence, alphabet="definite")
def pdist(seq1, seq2):
    r"""Calculate the *p*-distance between two aligned sequences.

    .. versionadded:: 0.7.2

    *p*-distance is the proportion of differing sites between two aligned sequences. It
    is equivalent to the normalized Hamming distance, but only considers definite
    characters (i.e., leaving out gaps and degenerate characters).

    .. math::
        p = \frac{\text{No. of differing sites}}{\text{Total no. of sites}}

    Parameters
    ----------
    seq1, seq2 : GrammaredSequence
        Sequences to compute the *p*-distance between.

    Returns
    -------
    float
        *p*-distance between the two sequences.

    Raises
    ------
    See ``hamming``.

    See Also
    --------
    hamming
    jc69

    Notes
    -----
    *p*-distance is the simplest measurement of the evolutionary distance (number of
    substitutions per site) between two sequences. It is also referred to as the *raw
    distance*.

    *p*-distance effectively estimates the evolutionary distance between two closely
    related sequences, where the number of observed substitutions is small. However,
    it may underestimate the true evolutionary distance when the two sequences are
    divergent and substitutions became saturated. This limitation may be overcome by
    adopting metrics that correct for multiple putative substitutions per site (such
    as JC69) and other biases.

    This function should not be confused with the :func:`~scipy.spatial.distance.pdist`
    function of SciPy.

    Examples
    --------
    >>> from skbio.sequence import DNA
    >>> from skbio.sequence.distance import pdist
    >>> seq1 = DNA('AGGGTA')
    >>> seq2 = DNA('CGTTTA')
    >>> pdist(seq1, seq2)
    0.5

    """
    return _pdist(np.vstack((seq1._bytes, seq2._bytes)), None, None).item()


def _pdist(seqs, mask, seqtype):
    """Compute pairwise p-distances between sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Sequences to compute the Hamming distance between.
    mask : ndarray of bool of shape (n_sequences, n_positions)
        Boolean mask of valid sites (True). If provided, each sequence pair will be
        filtered to positions that are valid in both of them.
    seqtype : type
        A placeholder. p-distance is irrelevant to sequence type.

    Returns
    -------
    ndarray of float of shape (C(n_sequences, 2),)
        p-distance matrix in condensed form.

    """
    npos = seqs.shape[1]

    def func(seqs, mask, i, out):
        subs = seqs[i] != seqs[i + 1 :]
        if mask is None:
            L = npos
            filled = True
        else:
            sites = mask[i] & mask[i + 1 :]
            subs &= sites
            L = np.count_nonzero(sites, axis=1)
            filled = L > 0
            if np.all(filled):
                filled = True

        # count substitutions
        p = np.count_nonzero(subs, axis=1)

        # normalize by sequence length
        np.divide(p, L, out=out, where=filled)

        if filled is not True:
            out[~filled] = np.nan

    return _build_dm(func, seqs, mask)


@_metric_specs(equal=True, seqtype=(DNA, RNA), alphabet="definite")
def jc69(seq1, seq2):
    r"""Calculate the JC69 distance between two aligned nucleotide sequences.

    .. versionadded:: 0.7.2

    The Jukes-Cantor 1969 (JC69) model estimates the evolutionary distance (number of
    substitutions per site) between two nucleotide sequences by correcting the observed
    proportion of differing sites (i.e., *p*-distance) to account for multiple putative
    substitutions at the same site (i.e., saturation). It is calculated as:

    .. math::
        D = -\frac{3}{4} ln(1 - \frac{4}{3} p)

    Parameters
    ----------
    seq1, seq2 : {DNA, RNA}
        Sequences to compute the JC69 distance between.

    Returns
    -------
    float
        JC69 distance between the two sequences.

    See Also
    --------
    jc69_correct
    pdist
    f81

    Notes
    -----
    The Jukes-Cantor 1969 (JC69) model was originally described in [1]_.

    JC69 is a basic evolutionary model for nucleotide sequences. It assumes equal base
    frequencies and equal substitution rates between bases. It models sequence
    evolution as a continuous-time Markov chain, and corrects the observed distance
    (*p*-distance) for repeated substitutions to estimate the true distance.

    This function returns NaN if :math:`p \geq 0.75`. This happens when the two
    sequences are too divergent and substitutions are over-saturated for reliable
    estimation of the evolutionary distance.

    References
    ----------
    .. [1] Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules.
       Mammalian Protein Metabolism, 3(21), 132.

    """
    return _jc69(np.vstack((seq1._bytes, seq2._bytes)), None, None).item()


def _jc69(seqs, mask, seqtype, chars=4):
    """Compute pairwise JC69 distances between sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Sequences to compute the Hamming distance between.
    mask : ndarray of bool of shape (n_sequences, n_positions)
        Boolean mask of valid sites (True). If provided, each sequence pair will be
        filtered to positions that are valid in both of them.
    seqtype : type
        A placeholder. p-distance is irrelevant to sequence type. Calculation of JC69
        distance does not distinguish between DNA and RNA sequences.
    chars : int, optional
        Number of definite characters in the alphabet. Default is 4 (for nucleotides).

    Returns
    -------
    ndarray of float of shape (C(n_sequences, 2),)
        JC69 distance matrix in condensed form.

    """
    dm = _pdist(seqs, mask, None)
    _p_correct(dm, 0.75)
    return dm


def jc69_correct(dists, chars=4, inplace=False):
    r"""Perform Jukes-Cantor (JC69) correction of a raw distance.

    .. versionadded:: 0.7.2

    The Jukes-Cantor 1969 (JC69) model estimates the evolutionary distance (number of
    substitutions per site) between two sequences by correcting the observed proportion
    of differing sites (*p*-distance) to account for multiple putative substitutions at
    the same site (i.e., saturation). The distance is calculated as:

    .. math::
        D = -\alpha ln(1 - \frac{p}{\alpha})

    Where :math:`\alpha = \frac{n - 1}{n}`, in which :math:`n` is the number of
    characters in the alphabet.

    Parameters
    ----------
    dists : float or array_like
        Uncorrected normalized Hamming distances (a.k.a., *p*-distance) between pairs
        of sequences.
    chars : int, optional
        Number of definite characters in the alphabet. Default is 4, which is for
        nucleotide sequences ("A", "C", "G" and "T/U").
    inplace : bool, optional
        If True and ``dists`` is an array, perform correction in place. Otherwise,
        create a copy of the data. Default is False.

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
    pdist

    Notes
    -----
    JC69 is the most basic evolutionary model, assuming equal character frequencies and
    equal substitution rates between characters.

    The JC69 model was developed for nucleotide sequences (4 characters, therefore
    :math:`\alpha = 0.75`). Technically, it is possible to apply JC69 correction on
    protein sequences, with ``chars=20``. However, this estimation may be inaccurate
    due to the profound variation of amino acid frequencies and substitution rates.

    This function returns NaN if :math:`p \geq \alpha`. This happens when the two
    sequences are too divergent and substitutions are over-saturated for reliable
    estimation of the true evolutionary distance.

    This function assumes that ``dists`` are valid *p*-distances within the range of
    [0, 1]. No validation will be performed. Input values outside this range will
    result in unexpected outputs.

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

    if inplace and isinstance(dists, np.ndarray):
        arr = dists
    else:
        arr = np.array(dists, copy=True)

    # perform correction
    _p_correct(arr, frac)

    if is_scalar:
        return arr.item()
    return arr


def _p_correct(dists, frac):
    """Correct p-distances in place using a JC69-like equation."""
    # There are three categories of values:
    # 1. Values exceeding frac become NaN, as these sites are saturated and the
    #    equation cannot estimate the true substitution number.
    # 2. Values == 0 remain 0.
    # 3. Values that suffice 0 < x < frac are subject to the equation.
    nonzero = dists > 0
    saturated = dists >= frac

    if np.any(saturated):
        dists[saturated] = np.nan
        all_valid = False
        if np.all(nonzero):
            valid = ~saturated
        else:
            valid = nonzero & ~saturated
    else:
        if np.all(nonzero):
            all_valid = True
            valid = None
        else:
            all_valid = False
            valid = nonzero

    if all_valid:
        vals = dists
    else:
        vals = dists[valid]

    # apply correction in place
    vals /= -frac
    vals += 1.0
    np.log(vals, out=vals)
    vals *= -frac

    if not all_valid:
        dists[valid] = vals


@_metric_specs(equal=True, seqtype=(DNA, RNA), alphabet="definite")
def f81(seq1, seq2, freqs=None):
    r"""Calculate the F81 distance between two aligned nucleotide sequences.

    .. versionadded:: 0.7.2

    The Felsenstein 1981 (F81) model assumes equal substitution rates and allows
    differential base frequencies (:math:`\pi`). The distance is calculated as:

    .. math::
        D = -\alpha ln(1 - \frac{p}{\alpha})

    Where :math:`p` is the proportion of differing sites (i.e., *p*-distance). Factor
    :math:`\alpha` is calculated as:

    .. math::
        \alpha = 1 - \pi_A^2 - \pi_C^2 - \pi_G^2 - \pi_T^2

    Parameters
    ----------
    seq1, seq2 : {DNA, RNA}
        Sequences to compute the F81 distance between.
    freqs : array_like of float of shape (4,), optional
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1. If not provided, the observed frequencies from the two input sequences
        combined will be used.

    Returns
    -------
    float
        F81 distance between the two sequences.

    See Also
    --------
    jc69
    f81

    Notes
    -----
    The Felsenstein 1981 (F81) model was described in [1]_ in the context of maximum
    likelihood estimation. The above equation for F81 distance calculation was adopted
    from [2]_, which is consistent with the implementation in ``ape::dist.dna``. The
    same equation was also described in [3]_ under the equal-input model.

    F81 is an extension of the JC69 model by allowing varying base frequencies. When
    the observed or user-provided based frequencies are equal (e.g., by specifying
    ``freqs=(.25, .25, .25, .25)``), the result will be identical to that of JC69.

    This function returns NaN if :math:`p \geq \alpha`.

    References
    ----------
    .. [1] Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
       likelihood approach. Journal of Molecular Evolution, 17(6), 368-376.

    .. [2] McGuire, G., Prentice, M. J., & Wright, F. (1999). Improved error bounds for
       genetic distances from DNA sequences. Biometrics, 55(4), 1064-1070.

    .. [3] Tajima, F., & Nei, M. (1984). Estimation of evolutionary distance between
       nucleotide sequences. Molecular Biology and Evolution, 1(3), 269-285.

    """
    return _f81(
        np.vstack((seq1._bytes, seq2._bytes)), None, type(seq1), freqs=freqs
    ).item()


def _f81(seqs, mask, seqtype, freqs):
    """Compute pairwise F81 distances between sequences."""
    frac = 1.0 - np.sum(np.asarray(freqs) ** 2)
    arr = _pdist(seqs, mask, None)
    _p_correct(arr, frac)
    return arr


@_metric_specs(equal=True, seqtype=(DNA, RNA), alphabet="definite")
def k2p(seq1, seq2):
    r"""Calculate the K2P distance between two aligned nucleotide sequences.

    .. versionadded:: 0.7.2

    The Kimura 2-parameter (K2P, a.k.a. K80) model allows differential rates of
    transitions (substitutions between two purines or between two pyrimidines) versus
    transversions (substitutions between a purine and a pyrimidine), while assuming
    equal base frequencies. The distance is calculated as:

    .. math::
        D = -\frac{1}{2} ln((1 - 2P - Q) \sqrt{1 - 2Q})

    Where :math:`P` and :math:`Q` are the proportions of transitions and transversions,
    respectively.

    Parameters
    ----------
    seq1, seq2 : {DNA, RNA}
        Sequences to compute the K2P distance between.

    Returns
    -------
    float
        K2P distance between the two sequences.

    See Also
    --------
    jc69
    f84

    Notes
    -----
    The Kimura 2-parameter model (K2P or K80) was originally described in [1]_.

    K2P is an extension of the JC69 model by modeling differential transition and
    transversion rates. Meanwhile, K2P can be considered as a special case of the F84
    model by assuming equal base frequencies.

    This function returns NaN if either :math:`1 - 2P - Q` or :math:`1 - 2Q` is zero or
    negative, which implicates over-saturation of substitutions.

    References
    ----------
    .. [1] Kimura, M. (1980). A simple method for estimating evolutionary rates of base
       substitutions through comparative studies of nucleotide sequences. Journal of
       Molecular Evolution, 16(2), 111-120.

    """
    return _k2p(np.vstack((seq1._bytes, seq2._bytes)), None, type(seq1)).item()


def _k2p(seqs, mask, seqtype):
    """Compute pairwise K2P distances between sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Sequences to compute the K2P distance between.
    mask : ndarray of bool of shape (n_sequences, n_positions)
        Boolean mask of valid sites (True). If provided, each sequence pair will be
        filtered to positions that are valid in both of them.
    seqtype : type
        Sequence type (DNA or RNA). Default is DNA.

    Returns
    -------
    ndarray of float of shape (C(n_sequences, 2),)
        K2P distance matrix in condensed form.

    """
    # This function uses the absolute difference between two ASCII codes to determine
    # their substitution type. Match: 0. Transition: A<->G: 6, C<->T: 17, C<->U: 18.
    # Transversion: all others. C = 17 (DNA) or 18 (RNA).
    C = 17 + issubclass(seqtype, RNA)

    npos = seqs.shape[1]
    logL_75 = 0.75 * np.log(npos)

    def func(seqs, mask, i, out):
        # Identify substitutions (difference != 0).
        # Casting to int16 because the subtraction of two uint8's can underflow.
        subs = np.abs(np.subtract(seqs[i], seqs[i + 1 :], dtype=np.int16))

        # Pairwise deletion of masked sites.
        if mask is not None:
            sites = mask[i] & mask[i + 1 :]
            subs[~sites] = 0
            L = np.count_nonzero(sites, axis=1)

        else:
            L = npos

        # substitutions
        p = np.count_nonzero(subs, axis=1)
        ident = p == 0

        # transitions (ts) (P = ts / L)
        ts = np.count_nonzero(subs == 6, axis=1) + np.count_nonzero(subs == C, axis=1)

        # transversions (tv) (Q = tv / L)
        tv = p - ts

        # convert into L * (1 - 2P - Q)
        ts *= -2
        ts -= tv
        ts += L

        # convert into L * (1 - 2Q)
        tv *= -2
        tv += L

        # Four groups of sequence pairs are to be treated differentially:
        #   1. Empty sequences: K = NaN
        #   2. Identical sequences: K = 0
        #   3. Overly divergent sequences: K = NaN
        #   4. Others (normal case): use the equation

        # Identify normal pairs (4).
        valid = (ts > 0) & (tv > 0) & ~ident

        has_inval = not np.all(valid)
        if has_inval:
            ts, tv = ts[valid], tv[valid]
            if mask is not None:
                ident &= L != 0
                L = L[valid]

        # Simplified calculation from the original formula:
        # K = -0.5 ln((1 - 2P - Q) * sqrt(1 - 2Q))
        #   = -0.5 ln((L - 2ts - tv) / L * ((L - 2tv) / L)**0.5)
        #   = -0.5 (ln(L - 2ts - tv) + 0.5 ln(L - 2tv) - 1.5 ln(L))
        #   = 0.75 ln(L) - 0.5 ln(L - 2ts - tv) - 0.25 ln(L - 2tv)
        c0 = logL_75 if mask is None else 0.75 * np.log(L)
        res = c0 - 0.5 * np.log(ts) - 0.25 * np.log(tv)

        if has_inval:
            out[valid] = res
            out[~valid] = np.nan
            out[ident] = 0.0
        else:
            out[:] = res

    return _build_dm(func, seqs, mask)


@_metric_specs(equal=True, seqtype=(DNA, RNA), alphabet="definite")
def f84(seq1, seq2, freqs=None):
    r"""Calculate the F84 distance between two aligned nucleotide sequences.

    .. versionadded:: 0.7.2

    The Felsenstein 1984 (F84) model allows differential rates of transitions and
    transversions, and differential base frequencies (:math:`\pi`). The distance is
    calculated as:

    .. math::
        D = -2Aln(1 - \frac{1}{2A}P - \frac{A-B}{2AC}Q) + 2(A-B-C)ln(1-\frac{1}{2C}Q)

    Where :math:`P` and :math:`Q` are the proportions of transitions and transversions,
    respectively. And:

    .. math::
        \begin{align}
        &A = \frac{\pi_A\pi_G}{\pi_A+\pi_G} + \frac{\pi_C\pi_T}{\pi_C+\pi_T} \\
        &B = \pi_A\pi_G + \pi_C\pi_T \\
        &C = (\pi_A+\pi_G)(\pi_C+\pi_T)
        \end{align}

    Parameters
    ----------
    seq1, seq2 : {DNA, RNA}
        Sequences to compute the F84 distance between.
    freqs : array_like of float of shape (4,), optional
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1. If not provided, the observed frequencies from the two input sequences
        combined will be used.

    Returns
    -------
    float
        F84 distance between the two sequences.

    See Also
    --------
    k2p
    tn93

    Notes
    -----
    The Felsenstein 1984 (F84) model for calculating sequence distance was initially
    implemented in the Phylip package [1]_. The model was then described in [2]_ and
    [3]_. The above equation was adopted from [4]_, which is consistent with the
    implementation in ``ape::dist.dna``.

    F84 is an extension of the K2P model that allows unequal base frequencies. When
    the observed or user-provided based frequencies are equal (e.g., by specifying
    ``freqs=(.25, .25, .25, .25)``), the result will be identical to that of K2P.

    F84 may be considered as a special case of the TN93 model where purine and
    pyrimidine transition rates are equal.

    This function returns NaN if either of the two logarithm arguments is zero or
    negative.

    References
    ----------
    .. [1] https://phylipweb.github.io/phylip/doc/dnadist.html

    .. [2] Kishino, H., & Hasegawa, M. (1989). Evaluation of the maximum likelihood
       estimate of the evolutionary tree topologies from DNA sequence data, and the
       branching order in Hominoidea. Journal of molecular evolution, 29(2), 170-179.

    .. [3] Felsenstein, J., & Churchill, G. A. (1996). A Hidden Markov Model approach
       to variation among sites in rate of evolution. Molecular Biology and Evolution,
       13(1), 93-104.

    .. [4] McGuire, G., Prentice, M. J., & Wright, F. (1999). Improved error bounds for
       genetic distances from DNA sequences. Biometrics, 55(4), 1064-1070.

    """
    return _f84(
        np.vstack((seq1._bytes, seq2._bytes)), None, type(seq1), freqs=freqs
    ).item()


def _f84(seqs, mask, seqtype, freqs):
    """Compute pairwise F84 distances between sequences."""
    npos = seqs.shape[1]

    piA, piC, piG, piT = freqs

    piR = piA + piG
    piY = piC + piT
    piAxG = piA * piG
    piCxT = piC * piT

    # intermediates (Eq. 2 of McGuire et al., 1999)
    A = piAxG / piR + piCxT / piY
    B = piAxG + piCxT
    C = piR * piY

    A_inv = 1.0 / A
    C_inv = 1.0 / C
    BpC = B + C
    AmB_AxC = (A - B) * A_inv * C_inv
    BCx2 = 2.0 * BpC

    nposx2 = 2.0 * npos
    log2 = np.log(2.0)

    # coefficients
    part0 = BCx2 * (np.log(npos) + log2)
    c1 = 2.0 * A
    c2 = 2.0 * BpC - c1

    X = 17 + issubclass(seqtype, RNA)
    masked = mask is not None

    def func(seqs, mask, i, out):
        subs = np.abs(np.subtract(seqs[i], seqs[i + 1 :], dtype=np.int16))

        if masked:
            sites = mask[i] & mask[i + 1 :]
            subs[~sites] = 0
            L = np.count_nonzero(sites, axis=1)
            Lx2 = 2.0 * L
        else:
            L = npos
            Lx2 = nposx2

        # substitutions
        p = np.count_nonzero(subs, axis=1)
        ident = p == 0

        # transitions
        P = np.count_nonzero(subs == 6, axis=1) + np.count_nonzero(subs == X, axis=1)

        # transversions
        Q = p - P

        # arguments (must be positive) * 2L
        a1 = Lx2 - A_inv * P - AmB_AxC * Q
        a2 = Lx2 - C_inv * Q

        # identify normal cases
        # There is no need to test L > 0. If L = 0, a1 and a2 are guaranteed to be 0,
        # per IEEE-754 (0 times any finite number is exactly 0).
        valid = (a1 > 0) & (a2 > 0) & ~ident

        has_inval = not np.all(valid)
        if has_inval:
            a1, a2 = a1[valid], a2[valid]
            if masked:
                ident &= L != 0
                L = L[valid]

        c0 = BCx2 * (np.log(L) + log2) if masked else part0

        # The formula (Eq. 3 of McGuire et al., 1999)
        res = c0 - c1 * np.log(a1) - c2 * np.log(a2)

        if has_inval:
            out[valid] = res
            out[~valid] = np.nan
            out[ident] = 0.0
        else:
            out[:] = res

    return _build_dm(func, seqs, mask)


@_metric_specs(equal=True, seqtype=(DNA, RNA), alphabet="definite")
def tn93(seq1, seq2, freqs=None):
    r"""Calculate the TN93 distance between two aligned nucleotide sequences.

    .. versionadded:: 0.7.2

    The Tamura and Nei (1993) (TN93) model assumes differential rates of the two
    types of transitions: between purines (R) (i.e., A <-> G) and between pyrimidines
    (Y) (i.e., C <-> T/U), and transversions (i.e., between a purine and a pyrimidine).
    It also allows varying base frequencies (:math:`\pi`). The distance is calculated
    as:

    .. math::
        \begin{align}
        D = &-2\frac{\pi_A\pi_G}{\pi_R}
            ln(1-\frac{\pi_R}{2\pi_A\pi_G}P_1-\frac{1}{2\pi_R}Q) \\
            &-2\frac{\pi_C\pi_T}{\pi_Y}
            ln(1-\frac{\pi_Y}{2\pi_C\pi_T}P_2-\frac{1}{2\pi_Y}Q) \\
            &-2(\pi_R\pi_Y-\frac{\pi_A\pi_G\pi_Y}{\pi_R}-\frac{\pi_C\pi_T\pi_R}{\pi_Y})
            ln(1-\frac{1}{2\pi_R\pi_Y}Q)
        \end{align}

    Where :math:`P_1` and :math:`P_2` are the proportions of purine and pyrimidine
    transitions, respectively. :math:`Q` is the proportion of transversions.

    Parameters
    ----------
    seq1, seq2 : {DNA, RNA}
        Sequences to compute the TN93 distance between.
    freqs : array_like of float of shape (4,), optional
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1. If not provided, the observed frequencies from the two input sequences
        combined will be used.

    Returns
    -------
    float
        TN93 distance between the two sequences.

    See Also
    --------
    k2p
    f84

    Notes
    -----
    The Tamura and Nei 1993 (TN93) model was originally described in [1]_.

    This function returns NaN if any of the three logarithm arguments is zero or
    negative, which implicates over-saturation of substitutions.

    References
    ----------
    .. [1] Tamura, K., & Nei, M. (1993). Estimation of the number of nucleotide
       substitutions in the control region of mitochondrial DNA in humans and
       chimpanzees. Molecular Biology and Evolution, 10(3), 512-526.

    """
    return _tn93(
        np.vstack((seq1._bytes, seq2._bytes)), None, type(seq1), freqs=freqs
    ).item()


def _tn93(seqs, mask, seqtype, freqs):
    """Compute pairwise TN93 distances between sequences.

    Parameters
    ----------
    seqs : ndarray of uint8 of shape (n_sequences, n_positions)
        Sequences to compute the TN93 distance between.
    mask : ndarray of bool of shape (n_sequences, n_positions)
        Boolean mask of valid sites (True). If provided, each sequence pair will be
        filtered to positions that are valid in both of them.
    seqtype : type
        Sequence type (DNA or RNA). Default is DNA.
    freqs : array_like of float of shape (4,)
        Relative frequencies of nucleobases A, C, G, and T/U.

    Returns
    -------
    ndarray of float of shape (C(n_sequences, 2),)
        TN93 distance matrix in condensed form.

    """
    npos = seqs.shape[1]

    masked = mask is not None

    piA, piC, piG, piT = freqs

    # frequencies of purines (R) and pyrimidines (Y)
    piR = piA + piG
    piY = piC + piT

    piAxG = piA * piG
    piCxT = piC * piT
    piRxY = piR * piY
    piR_AxG = piR / piAxG
    piY_CxT = piY / piCxT

    # coefficients
    c1 = 2.0 * piAxG / piR
    c2 = 2.0 * piCxT / piY
    c3 = 2.0 * (piRxY - piY / piR_AxG - piR / piY_CxT)

    c123 = c1 + c2 + c3
    nposx2 = 2.0 * npos
    log2 = np.log(2.0)
    part0 = c123 * (np.log(npos) + log2)

    X = 17 + issubclass(seqtype, RNA)

    def func(seqs, mask, i, out):
        subs = np.abs(np.subtract(seqs[i], seqs[i + 1 :], dtype=np.int16))

        if masked:
            sites = mask[i] & mask[i + 1 :]
            subs[~sites] = 0
            L = np.count_nonzero(sites, axis=1)
            Lx2 = 2.0 * L
        else:
            L = npos
            Lx2 = nposx2

        # substitutions
        p = np.count_nonzero(subs, axis=1)
        ident = p == 0

        # purine transitions
        P1 = np.count_nonzero(subs == 6, axis=1)

        # pyrimidine transitions
        P2 = np.count_nonzero(subs == X, axis=1)

        # transversions
        Q = p - P1 - P2

        # arguments (must be positive) * 2L
        Lx2 = 2.0 * L
        a1 = Lx2 - P1 * piR_AxG - Q / piR
        a2 = Lx2 - P2 * piY_CxT - Q / piY
        a3 = Lx2 - Q / piRxY

        # identify normal cases
        valid = (a1 > 0) & (a2 > 0) & (a3 > 0) & ~ident

        has_inval = not np.all(valid)
        if has_inval:
            a1, a2, a3 = a1[valid], a2[valid], a3[valid]
            if masked:
                ident &= L != 0
                L = L[valid]

        c0 = c123 * (np.log(L) + log2) if masked else part0

        # The formula
        res = c0 - c1 * np.log(a1) - c2 * np.log(a2) - c3 * np.log(a3)

        if has_inval:
            out[valid] = res
            out[~valid] = np.nan
            out[ident] = 0.0
        else:
            out[:] = res

    return _build_dm(func, seqs, mask)


@_metric_specs(equal=True, seqtype=GrammaredSequence, alphabet="canonical")
def logdet(seq1, seq2, pseudocount=None):
    r"""Calculate the LogDet distance between two aligned sequences.

    .. versionadded:: 0.7.2

    The LogDet estimator of evolutionary distance is robust to compositional biases and
    nonstationary (i.e., changing over time) character frequencies. The distance is
    calculated as:

    .. math::
        D = -\frac{1}{k}ln det \mathbf{F} - ln k

    Where :math:`\mathbf{F}` is a :math:`k \times k` matrix of proportions of character
    pairs between the two sequences. :math:`k` is the number of canonical characters in
    the alphabet of the specific sequence type (e.g., 4 for nucleotide and 20 for
    protein).

    Parameters
    ----------
    seq1, seq2 : GrammaredSequence
        Sequences to compute the LogDet distance between.
    pseudocount : float, optional
        A small positive value added to the count of each character pair to avoid
        logarithm of zero. Can prevent a NaN result when some character pairs are
        missing from the sequences. Default is None.

    Returns
    -------
    float
        LogDet distance between the two sequences.

    See Also
    --------
    paralin

    Notes
    -----
    The LogDet (meaning "logarithm of determinant") transformation was originally
    described in [1]_. The above equation of LogDet distance was adopted from [2]_,
    which is consistent with the implementation in ``ape::dist.dna``.

    The function returns NaN when the determinant is 0 or negative.

    Although the LogDet distance is mostly used for analyzing nucleotide sequences (_k_
    = 4), this function is applicable to any grammared sequences with an arbitrarily
    sized alphabet, in accordance with the original paper [1]_.

    However, a large alphabet (e.g., protein sequences, with _k_ = 20) often results
    in a sparse _F_ matrix, due to unobserved character pairs in the sequences.
    Consequently, the LogDet distance will be NaN. This can be mitigated by specifying
    a pseudocount (e.g., 0.5) to regularize the matrix.

    Additionally, the LogDet distance tends to over-estimate the evolutionary distance
    when the character frequencies are highly unequal. See [3]_ for a discussion and
    a modified LogDet distance to account for unequal character frequencies.

    The LogDet distance between two identical sequences is 0 only when the character
    frequencies are equal, which is rarely the case in real data. However, the
    ``DistanceMatrix`` data structure forces the diagonal (distance from each sequence
    to itself) to be 0. Be cautious of this when self-distance is involved in the
    subsequent analysis.

    References
    ----------
    .. [1] Lockhart, P. J., Steel, M. A., Hendy, M. D., & Penny, D. (1994). Recovering
       evolutionary trees under a more realistic model of sequence evolution. Molecular
       Biology and Evolution, 11(4), 605-612.

    .. [2] Gu, X., & Li, W. H. (1996). Bias-corrected paralinear and LogDet distances
       and tests of molecular clocks and phylogenies under nonstationary nucleotide
       frequencies. Molecular Biology and Evolution, 13(10), 1375-1383.

    .. [3] Tamura, K., & Kumar, S. (2002). Evolutionary distance estimation under
       heterogeneous substitution pattern among lineages. Molecular Biology and
       Evolution, 19(10), 1727-1736.

    """
    return _logdet(
        np.vstack((seq1._bytes, seq2._bytes)), None, type(seq1), pseudocount
    ).item()


def _logdet(seqs, mask, seqtype, pseudocount=None):
    """Compute pairwise LogDet distances between sequences."""
    nchar, index = _char_index(seqtype)
    seqs = index[seqs]

    masked = mask is not None
    n = seqs.shape[0]
    n_1 = n - 1

    # number of valid characters + others (see below)
    k = nchar + masked

    c1 = -1.0 / nchar
    part0 = np.log(nchar)

    # an offset added to sequence 2's to facilitate bincount (see below)
    offvec = np.arange(n).reshape(-1, 1) * k**2

    def func(seqs, mask, i, out):
        # number of sequence pairs
        m = n_1 - i

        # The following code computes the count matrix between all sequence pairs.
        # `np.bincount` is highly efficient for counting, but it only operates on a 1D
        # array. To overcome this, two tricks are applied:
        # 1) For each pair of sequences (i vs j), compute i * k + j, do bincount, then
        #    reshape to (k, k), which is the count matrix.
        # 2) For all pairs of sequences (i vs J), offset J by k^2 each, flatten, do
        #    bincount, reshape to (-1, ...), which returns the counts of each sequence
        #    pair.
        # The outcome is a 3D array of shape (m, k, k).
        pairs = seqs[i] * k + seqs[i + 1 :] + offvec[:m]
        mat = np.bincount(pairs.ravel(), minlength=offvec[m, 0]).reshape(-1, k, k)

        # In pairwise deletion mode, all characters other than the canonical ones are
        # given the same and maximum code, which will be trimmed off from the count
        # matrix.
        # Complete deletion mode does not need this as all characters are canonical,
        # which saves some compute (25 -> 16 for nucleotide sequences).
        if masked:
            mat = mat[:, :nchar, :nchar]

        mat = mat.astype(float)

        if pseudocount:
            mat += pseudocount

        total = mat.sum(axis=(1, 2), keepdims=True)

        # drop empty sequences
        # Note that if a pseudocount is added, empty sequences become non-empty.
        has_empty = False
        if masked and not pseudocount:
            filled = (total > 0).ravel()
            if has_empty := not np.all(filled):
                mat, total = mat[filled], total[filled]

        # convert counts to proportions
        mat /= total

        # `np.logabsdet` is more stable than `det` -> `log`.
        sign, logabsdet = np.linalg.slogdet(mat)

        # The determinant should be positive (sign = 1). Otherwise result is NaN.
        valid = sign == 1.0
        if has_inval := not np.all(valid):
            logabsdet = logabsdet[valid]

        # The formula (Eq. 2 of Gu & Li, 1996)
        res = c1 * logabsdet - part0

        if has_empty:
            if has_inval:
                filled[filled] = valid  # combine Boolean masks
            valid = filled
            has_inval = True
        elif not has_inval:
            valid = slice(None)

        out[valid] = res
        if has_inval:
            out[~valid] = np.nan

    return _build_dm(func, seqs, None)


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
        *k*-mer distance between the two sequences.

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
