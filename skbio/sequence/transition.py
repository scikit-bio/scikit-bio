"""Transition probability models (:mod:`skbio.sequence.transition`)
================================================================

.. currentmodule:: skbio.sequence.transition

This module provides functions for calculating transition probability matrices (TPMs,
a.k.a. substitution probability matrices) under different substitution models for a
specified evolutionary distance (expected number of substitutions per site, which is
often represented as branch length in a phylogenetic tree).

A TPM gives the probability that each ancestral state (rows) is observed as each
descendant state (columns) after the specified evolutionary distance. For
continuous-time Markov models, the TPM is obtained from the instantaneous rate matrix,
:math:`Q`, as :math:`P(t) = e^{Qt}`. Each row of the matrix sums to one.

Models differ in the assumptions they make about instantaneous substitution rates and
equilibrium state frequencies, resulting in different sets of model parameters.

Transition probability matrices
-------------------------------

.. autosummary::
   :toctree:

   jc69
   k2p
   f81
   hky85
   tn93

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from skbio.sequence.distance import _check_freqs
from skbio.sequence import SubstitutionMatrix
import skbio.sequence as sk_seqtype

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Callable
    from numpy.typing import ArrayLike, NDArray


def _tpm_wrapper(
    func: Callable,
    d: float,
    seqtype: str,
    allowed_seqtypes: str | tuple,
    *args: tuple,
    **kwargs: dict,
):
    """Wrapper for transition probability model functions.

    This helper function transforms distance `d` so it is always is a vector for the
    internal interface and transition probability matrix output so it is always a valid
    `SubstitutionMatrix` instance.

    Parameters
    ----------
    func : callable
        Transition probability model function.
    d : float
        Evolutionary distance between sequences.
    seqtype : type or tuple of types, optional
        Name of the sequence type. Should be one of values in `allowed_seqtypes`.
    allowed_seqtypes : string or tuple of strings, optional
        Valid sequence types for the model. Can be a single type (such as `"Protein"`)
        or a tuple of multiple types (such as `("DNA", "RNA")`).
    kwargs : dict, optional
        Transition model-specific keyword parameters. Refer to the documentation of
        the chosen model.

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    """
    # Validate sequence type
    if not isinstance(allowed_seqtypes, tuple):
        allowed_seqtypes = (allowed_seqtypes,)

    if seqtype not in allowed_seqtypes:
        types = ", ".join(allowed_seqtypes)
        raise TypeError(f"Sequence should be: {types}. Got {seqtype} instead.")

    seqtype_cls = getattr(sk_seqtype, seqtype)

    # Get the states for the substitution matrix
    states = sorted(seqtype_cls.definite_chars)

    # Vectorize the distance (input: scalar, output: 1D array)
    d = np.atleast_1d(_check_d(d))

    # Call the transition probability model function to get the probabilities
    probs = func(d, *args, **kwargs)

    return SubstitutionMatrix(states, probs[0])


def _check_d(d: float):
    """Validate that the evolutionary distance `d` is non-negative."""
    if d < 0:
        raise ValueError("Distance must be non-negative.")
    return d


def _check_kappa(kappa: float):
    """Validate that the ts/tv ratio `kappa` is between 0 and 1."""
    if kappa <= 0.0 or kappa > 1.0:
        raise ValueError(
            "Parameter 'kappa' must be greater than 0 and less or equal to 1."
        )
    return kappa


def jc69(d: float, seqtype: str = "DNA") -> SubstitutionMatrix:
    r"""Calculate the JC69 transition probability matrix for a given distance.

    .. versionadded:: 0.7.4

    The Jukes-Cantor 1969 (JC69) model assumes equal nucleotide frequencies
    and equal substitution rates between all nucleotides.

    The transition probability between nucleotide states for sequences
    separated by evolutionary distance :math:`d` (expected substitutions
    per site) is given by:

    .. math::

        P_{ij} = \begin{cases}
            \frac{1}{4} + \frac{3}{4} e^{-\frac{4d}{3}} & i = j \\
            \frac{1}{4} - \frac{1}{4} e^{-\frac{4d}{3}} & i \neq j
        \end{cases}

    where :math:`i, j \in \{A, C, G, T/U\}`, and :math:`i` and :math:`j`
    denote ancestral and descendant states, respectively.

    Parameters
    ----------
    d : float
        Evolutionary distance between sequences.
    seqtype : {'DNA', 'RNA'}, optional
        Sequence type. Used to label matrix states as nucleotides (A, C, G, T/U).
        Default is "DNA".

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    See Also
    --------
    skbio.sequence.distance.jc69

    Notes
    -----
    The Jukes-Cantor 1969 (JC69) model was originally described in [1]_.

    It is a continuous-time Markov chain model that assumes equal base frequencies and
    equal substitution rates between all nucleotides.

    References
    ----------
    .. [1] Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules.
       Mammalian Protein Metabolism, 3(21), 132.

    """
    return _tpm_wrapper(_jc69, d, seqtype=seqtype, allowed_seqtypes=("DNA", "RNA"))


def _jc69(d: NDArray) -> NDArray:
    """Calculate the JC69 transition probability matrices for given distances.

    Parameters
    ----------
    d : ndarray of shape (n,)
        Evolutionary distances between sequences.

    Returns
    -------
    ndarray of shape (n, 4, 4)
        Transition probability matrices.

    """
    e = np.exp(-4.0 * d / 3.0)

    same = 0.25 * (1.0 + 3.0 * e)
    diff = 0.25 * (1.0 - e)

    P = np.empty((len(d), 4, 4), dtype=np.float64)
    P[:] = diff[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = same[:, None]

    return P


def k2p(d: float, kappa: float, seqtype: str = "DNA") -> SubstitutionMatrix:
    r"""Calculate the K2P transition probability matrix for a given distance.

    .. versionadded:: 0.7.4

    The Kimura 2-parameter (K2P, a.k.a. K80) model assumes separate rates for
    transitions :math:`(A \leftrightarrow G, C \leftrightarrow T/U)` and
    transversions (all other changes of substitutions), with a fixed
    transition/transversion rate ratio :math:`\kappa`.

    The transition probabilities between nucleotide states for sequences
    separated by evolutionary distance :math:`d` are given by:

    .. math::
        \begin{aligned}
        &P_{no\, change} = \frac{1}{4} + \frac{1}{4} e^{\frac{-4d}{3}}
            + \frac{1}{2} e^{\frac{-2(\kappa + 1)d}{3}} \\
        &P_{transition} = \frac{1}{4} + \frac{1}{4} e^{\frac{-4d}{3}}
            - \frac{1}{2} e^{\frac{-2(\kappa + 1)d}{3}} \\
        &P_{transversion} = \frac{1}{4} - \frac{1}{4} e^{\frac{-4d}{3}}
        \end{aligned}

    Where :math:`\kappa` is ratio between transition and transversion rates.

    Parameters
    ----------
    d : float
        Evolutionary distance between sequences.
    kappa : float
        Ratio between transition and transversion rates. Should be between 0 and 1.
    seqtype : {'DNA', 'RNA'}, optional
        Sequence type. Used to label matrix states as nucleotides (A, C, G, T/U).
        Default is "DNA".

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    See Also
    --------
    skbio.sequence.distance.k2p
    jc69

    Notes
    -----
    The Kimura 2-parameter model (K2P or K80) was originally described in [1]_.

    K2P is an extension of the JC69 model by modeling differential transition and
    transversion rates. Meanwhile, K2P can be considered as a special case of the F84
    model by assuming equal base frequencies.

    References
    ----------
    .. [1] Kimura, M. (1980). A simple method for estimating evolutionary rates of base
       substitutions through comparative studies of nucleotide sequences. Journal of
       Molecular Evolution, 16(2), 111-120.

    """

    return _tpm_wrapper(
        _k2p,
        d,
        seqtype=seqtype,
        allowed_seqtypes=("DNA", "RNA"),
        kappa=_check_kappa(kappa),
    )


def _k2p(d: NDArray, kappa: float) -> NDArray:
    """Calculate the K2P transition probability matrices for given distances.

    Parameters
    ----------
    d : ndarray of shape (n,)
        Evolutionary distances between sequences.
    kappa : float
        Ratio between transition and transversion rates.

    Returns
    -------
    ndarray of shape (n, 4, 4)
        Transition probability matrices.

    """
    e1 = np.exp(-4.0 * d / 3.0)
    e2 = np.exp(-2.0 * (kappa + 1) * d / 3.0)

    same = 0.25 + 0.25 * e1 + 0.5 * e2
    transition = 0.25 + 0.25 * e1 - 0.5 * e2
    transversion = 0.25 - 0.25 * e1

    P = np.empty((len(d), 4, 4), dtype=np.float64)

    P[:] = transversion[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = same[:, None]

    # transitions: A<->G and C<->T/U
    P[:, 0, 2] = transition
    P[:, 2, 0] = transition
    P[:, 1, 3] = transition
    P[:, 3, 1] = transition

    return P


def f81(d: float, freqs: ArrayLike, seqtype: str = "DNA") -> SubstitutionMatrix:
    r"""
    Calculate the F81 transition probability matrix for a given distance.

    .. versionadded:: 0.7.4

    The Felsenstein 1981 (F81) model assumes equal substitution rates among
    nucleotide pairs but allows unequal equilibrium base frequencies :math:`\pi`.

    The transition probability between nucleotide states for sequences
    separated by evolutionary distance :math:`d` is given by:

    .. math::
        P_{ij} = \begin{cases}
            \pi_i + (1 - \pi_i) e^{-\frac{d}{b}} & i = j \\
            \pi_j (1 - e^{-\frac{d}{b}}) & i \neq j
        \end{cases}

    Where :math:`i,j \in \{A, C, G, T/U\}`, :math:`i` is ancestral nucleotide,
    :math:`j` is descendant nucleotide. The normalization constant :math:`b`
    is defined as:

    .. math::
        b = 1 - \pi_A^2 - \pi_C^2 - \pi_G^2 - \pi_T^2

    Parameters
    ----------
    d : float
        Evolutionary distance between sequences.
    freqs : array_like of float of shape (4,)
        Relative frequencies of nucleobases A, C, G, and T/U, respectively.
        Should sum to 1.
    seqtype : {'DNA', 'RNA'}, optional
        Sequence type. Used to label matrix states as nucleotides (A, C, G, T/U).
        Default is "DNA".

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    See Also
    --------
    skbio.sequence.distance.f81
    jc69

    Notes
    -----
    The Felsenstein 1981 (F81) model was originally introduced in [1]_.

    It extends JC69 by allowing unequal equilibrium base frequencies while
    retaining equal substitution rates between all nucleotide pairs.

    F81 reduces to JC69 when all base frequencies are equal
    (:math:`\pi_i = 0.25`).


    References
    ----------
    .. [1] Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
       likelihood approach. Journal of Molecular Evolution, 17(6), 368-376.

    """

    return _tpm_wrapper(
        _f81,
        d,
        seqtype=seqtype,
        allowed_seqtypes=("DNA", "RNA"),
        freqs=_check_freqs(freqs),
    )


def _f81(d: NDArray, freqs: NDArray) -> NDArray:
    """Calculate the F81 transition probability matrices for given distances.

    Parameters
    ----------
    d : ndarray of shape (n,)
        Evolutionary distances between sequences.
    freqs : ndarray of shape (4,)
        Stationary frequencies of the four nucleotides.

    Returns
    -------
    ndarray of shape (n, 4, 4)
        Transition probability matrices.

    """
    freqs = np.asarray(freqs, dtype=np.float64)

    scale_factor = 1.0 / (1.0 - np.sum(freqs**2))

    e = np.exp(-scale_factor * d)

    n = len(d)
    P = np.empty((n, 4, 4), dtype=np.float64)

    P[:] = freqs[None, None, :] * (1.0 - e)[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = freqs[None, :] + (1.0 - freqs)[None, :] * e[:, None]

    return P


def hky85(
    d: float, freqs: ArrayLike, kappa: float, seqtype: str = "DNA"
) -> SubstitutionMatrix:
    r"""Calculate the HKY85 transition probability matrix for a given distance.

    .. versionadded:: 0.7.4

    The Hasegawa, Kishino, and Yano 1985 (HKY85) allows differential base
    frequencies (:math:`\pi`) and assumes fixed ratio :math:`\kappa` between
    rates of transitions (substitutions between two purines or between two
    pyrimidines) and transversions (substitutions between a purine and a
    pyrimidine). Transition probability for a nucleotide for sequences with
    a distance :math:`d` between them can be calculated like:

    .. math::
        P_{ij} = e^{\frac{-(\kappa + 1)d}{2b}}\delta_{ij} +
        e^{-\frac{d}{b}}\left (1 - e^{\frac{-(\kappa - 1)d}{2b}}\right )
        \frac{\pi_j\epsilon_{ij}}{\sum_k \epsilon_{jk}\pi_k} +
        \left (1 - e^{-\frac{d}{b}}\right )\pi_j

    Where :math:`i,j \in \{A, C, G, T/U\}`, :math:`i` is ancestral nucleotide,
    :math:`j` is descendant nucleotide, :math:`\delta_{ij}` is Kronecker
    delta function and is 1 if :math:`i=j` and 0 otherwise, \epsilon_{ij}
    is purine/pyrimidine indicator function and is 1 if :math:`i` and :math:`j`
    in the same nucloetide class and is 0 otherwise. The normalization constant
    :math:`b` is defined as:

    .. math::
        b = 1 - \pi_A^2 - \pi_C^2 - \pi_G^2 - \pi_T^2

    Parameters
    ----------
    d : float
        Evolutionary distance between sequences.
    freqs : array_like of float of shape (4,)
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1.
    kappa : float
        Ratio between transition and transversion rates. Should be between 0 and 1.
    seqtype : {'DNA', 'RNA'}, optional
        Sequence type. Used to label matrix states as nucleotides (A, C, G, T/U).
        Default is "DNA".

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    See Also
    --------
    skbio.sequence.distance.f84
    f81
    tn93

    Notes
    -----
    The HKY85 model was described in [1]_. The transition probability formulas
    were described in [2]_.

    It is a generalization of F81 that incorporates transition/transversion
    bias via :math:`\kappa`. When :math:`\kappa = 1`, the model reduces to
    a model without transition/transversion bias and becomes equivalent to
    F81 under standard normalization.

    References
    ----------
    .. [1] Hasegawa, M., Kishino, H., & Yano, T. A. (1985). Dating of the human-ape
       splitting by a molecular clock of mitochondrial DNA. Journal of Molecular
       Evolution, 22(2), 160-174.
    .. [2] Felsenstein, J. (2004). Inferring Phylogenies. 2003. Sinauer Associates,
       Sunderland, Massachusetts.

    """

    return _tpm_wrapper(
        _hky85,
        d,
        seqtype=seqtype,
        allowed_seqtypes=("DNA", "RNA"),
        freqs=_check_freqs(freqs),
        kappa=_check_kappa(kappa),
    )


def _hky85(d: NDArray, freqs: NDArray, kappa: float) -> NDArray:
    r"""Calculate theHKY85 transition probability matrix for given distances.

    Parameters
    ----------
    d : ndarray of shape (n,)
        Evolutionary distances between sequences.
    freqs : ndarray of shape (4,)
        Stationary frequencies of the four nucleotides.
    kappa : float
        Ratio between transition and transversion rates.

    Returns
    -------
    ndarray of shape (n, 4, 4)
        Transition probability matrices.

    """
    return _tn93(d, freqs=freqs, kappa_r=kappa, kappa_y=kappa)


def tn93(
    d: float, freqs: ArrayLike, kappa_r: float, kappa_y: float, seqtype: str = "DNA"
) -> SubstitutionMatrix:
    r"""Calculate the TN93 transition probability matrix for a given distance.

    .. versionadded:: 0.7.4

    The Tamura-Nei 1993 (TN93) allows differential base frequencies (:math:`\pi`) and
    assumes different ratios :math:`\kappa` between rates of transitions (substitutions
    between two purines or between two pyrimidines) and transversions (substitutions
    between a purine and a pyrimidine) for each class of nucleotides. Transition
    probability for a nucleotide for sequences with a distance :math:`d` between them
    can be calculated as:

    .. math::
        P_{ij} = e^{-\frac{\kappa_i d}{2b}}\delta_{ij} +
            e^{\frac{-d}{b}}\left (1 - e^{-\frac{(\kappa_i - 1)d}{2b}}\right )
            \frac{\pi_j\epsilon_{ij}}{\sum_k \epsilon_{jk}\pi_k} +
            \left (1 - e^{\frac{-d}{b}}\right )\pi_j

    Where :math:`i,j \in \{A, C, G, T/U\}`, :math:`i` is ancestral nucleotide,
    :math:`j` is descendant nucleotide, :math:`\delta_{ij}` is Kronecker delta function
    and is 1 if :math:`i=j` and 0 otherwise, \epsilon_{ij} is purine/pyrimidine
    indicator function which is 1 if :math:`i` and :math:`j` are in the same nucloetide
    class and is 0 otherwise, :math:`\kappa_i` is a transition/transversion rate ratio
    for purines or pyrimidines, depending on a substitution type. The normalization
    constant :math:`b` is defined as:

    .. math::
        b = 1 - \pi_A^2 - \pi_C^2 - \pi_G^2 - \pi_T^2

    Parameters
    ----------
    d : float
        Evolutionary distance between sequences.
    freqs : array_like of float of shape (4,)
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1.
    kappa_r : float
        Transition/transversion rate ratio of purines. Should be between 0 and 1.
    kappa_y : float
        Transition/transversion rate ratio of pyrimidines. Should be between 0 and 1.
    seqtype : {'DNA', 'RNA'}, optional
        Sequence type. Used to label matrix states as nucleotides (A, C, G, T/U).
        Default is "DNA".

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides, columns are
        descendant nucleotides.

    See Also
    --------
    skbio.sequence.distance.tn93
    hky85

    Notes
    -----
    The Tamura-Nei 1993 (TN93) model was described in [1]_. The transition probability
    formulas were taken from [2]_.

    It generalizes HKY85 by allowing different transition rates for purine and
    pyrimidine substitutions while maintaining a single transversion rate class. HKY85
    is recovered when :math:`\kappa_R = \kappa_Y`.

    References
    ----------
    .. [1] Tamura, K., & Nei, M. (1993). Estimation of the number of nucleotide
       substitutions in the control region of mitochondrial DNA in humans and
       chimpanzees. Molecular Biology and Evolution, 10(3), 512-526.
    .. [2] Felsenstein, J. (2004). Inferring Phylogenies. 2003. Sinauer
       Associates, Sunderland, Massachusetts.

    """

    return _tpm_wrapper(
        _tn93,
        d,
        seqtype=seqtype,
        allowed_seqtypes=("DNA", "RNA"),
        freqs=_check_freqs(freqs),
        kappa_r=_check_kappa(kappa_r),
        kappa_y=_check_kappa(kappa_y),
    )


def _tn93(d: NDArray, freqs: NDArray, kappa_r: float, kappa_y: float) -> NDArray:
    r"""Tamura-Nei 1993 (TN93) transition probability matrix for branch length d.

    Parameters
    ----------
    d : ndarray of shape (n,)
        Evolutionary distances between sequences..
    freqs : ndarray of shape (4,)
        Stationary frequencies of the four nucleotides.
    kappa_r : float
        Transition/transversion rate ratio of purines.
    kappa_y : float
        Transition/transversion rate ratio of pyrimidines.

    """
    pi_R = freqs[0] + freqs[2]
    pi_Y = freqs[1] + freqs[3]
    adjusted_freqs = np.array(
        [freqs[0] / pi_R, freqs[1] / pi_Y, freqs[2] / pi_R, freqs[3] / pi_Y]
    )[None, None, :]
    freqs = freqs[None, None, :]

    # Note that kappa is vector here
    kappa = np.array([kappa_r, kappa_y] * 2)[None, :, None]

    same_class = np.fromfunction(lambda i, j: (i % 2) == (j % 2), (4, 4), dtype=int)[
        None, :, :
    ]

    identity = np.eye(4)[None, :, :]

    scale_factor = 1.0 / (1.0 - np.sum(freqs**2))

    e1 = np.exp(-scale_factor * (kappa + 1.0) * d / 2.0)
    e2 = np.exp(-scale_factor * d)[:, None, None]
    e3 = np.exp(-scale_factor * (kappa - 1.0) * d / 2.0)

    P = np.empty((len(d), 4, 4), dtype=np.float64)

    P = (
        e1 * identity
        + e2 * (1.0 - e3) * adjusted_freqs * same_class
        + (1.0 - e2) * freqs
    )

    return P
