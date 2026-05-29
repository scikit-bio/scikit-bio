"""Transition probability models (:mod:`skbio.tree.tpm`)
==========================================================

.. currentmodule:: skbio.tree.tpm

Lorem ipsum.


Nucleotide transition probability models
----------------------------------------

.. autosummary::
   :toctree:

   jc69
   k2p
   f81

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from typing import TYPE_CHECKING
from inspect import signature
import functools

import numpy as np
import numpy.typing as npt

from skbio.sequence.distance import _char_hash, _check_freqs
from collections.abc import Callable
from skbio.sequence import (
    Sequence,
    GrammaredSequence,
    DNA,
    RNA,
    Protein,
    SubstitutionMatrix,
)
import skbio.sequence as sk_seqtype

if TYPE_CHECKING:
    from numpy.typing import ArrayLike


def _wrap_vector_tpm_function(
    function: Callable,
    d: float,
    seqtype: str,
    allowed_seqtypes: str | tuple,
    *args: tuple,
    **kwargs: dict,
):
    """
    Transforms distance d so it is always is a vector for the internal
    interface and transition probability matrix output so it is always a
    valid PairwiseMatrix.
    Parameters
    ----------
    function : function
        Transition probability model function.
    d : float
        distance between sequences.
    kwargs : tuple, optional
        Transition model-specific positional parameters. Refer to
        the documentation of the chosen model.
    seqtype : type or tuple of types, optional
        Valid sequence types. Can be a single type (such as ``Protein``) or a tuple of
        multiple types (such as ``(DNA, RNA)``). If None (default), any ``Sequence``
        objects, grammared or not, are valid.
    kwargs : dict, optional
        Transition model-specific keyword parameters. Refer to
        the documentation of the chosen model.
    """

    d = np.atleast_1d(d)

    probabilities = function(d, *args, **kwargs)

    if not isinstance(allowed_seqtypes, tuple):
        allowed_seqtypes = (allowed_seqtypes,)

    if seqtype not in allowed_seqtypes:
        types = ", ".join(allowed_seqtypes)
        raise TypeError(f"Sequence should be: {types}. Got {seqtype} instead")

    seqtype_cls = getattr(sk_seqtype, seqtype)
    states = sorted(seqtype_cls.definite_chars)

    return SubstitutionMatrix(states, probabilities[0])


def jc69(
    d: float,
    seqtype: str = "DNA",
) -> SubstitutionMatrix:
    r"""
    Calculate the JC69 transition probability matrix for a given distance.

    .. versionadded:: 0.7.3

    The Jukes-Cantor 1969 (JC69) model assumes equal substitution
    rates between nucleotides. Transition probability for a
    nucleotide for sequences with distance :math:`d` between them can be
    calculated like:

    .. math::

        P_{ij} = \begin{cases}
            \frac{1}{4} + \frac{3}{4} e^{-\frac{4d}{3}} & i = j \\
            \frac{1}{4} - \frac{1}{4} e^{-\frac{4d}{3}} & i \neq j
        \end{cases}

    Where :math:`i,j \in \{A, C, G, T\}`, :math:`i` is ancestral nucleotide,
    :math:`j` is descendant nucleotide.

    Parameters
    ----------
    d : float
        distance or distances between sequences.
    seqtype : str
        String that holds type of the sequence. "DNA" (default) and "RNA"
        are valid options. Needed to assign proper nucleotide letters to
        PairwiseMatrix ids.

    Returns
    -------
    PairwiseMatrix
        Transition probability matrix. Rows are ancestral nucleotides,
        columns are descendant nucleotides.


    Notes
    -----
    The Jukes-Cantor 1969 (JC69) model was originally described in [1]_.

    JC69 is a basic evolutionary model for nucleotide sequences. It assumes equal base
    frequencies and equal substitution rates between bases. It models sequence
    evolution as a continuous-time Markov chain, and corrects the observed distance
    (*p*-distance) for repeated substitutions to estimate the true distance.

    References
    ----------
    .. [1] Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules.
       Mammalian Protein Metabolism, 3(21), 132.

    """

    return _wrap_vector_tpm_function(
        _jc69, d, allowed_seqtypes=("DNA", "RNA"), seqtype=seqtype
    )


def _jc69(
    d: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    e = np.exp(-4.0 * d / 3.0)

    same = 0.25 * (1.0 + 3.0 * e)
    diff = 0.25 * (1.0 - e)

    P = np.empty((len(d), 4, 4), dtype=np.float64)
    P[:] = diff[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = same[:, None]

    return P


def k2p(
    d: float,
    alpha: float,
    beta: float,
    seqtype: str = "DNA",
) -> SubstitutionMatrix:
    r"""
    Calculate the K2P transition probability matrix for a given distance.

    .. versionadded:: 0.7.3

    The Kimura 2-parameter (K2P, a.k.a. K80) model assumes separate rates for
    transitions :math:`(A \leftrightarrow G, C \leftrightarrow T)` and transversions
    (all other changes of nucleotides). Transition probability for a
    nucleotide for sequences with distance :math:`d` between them can be
    calculated like:

    .. math::
        \begin{aligned}
        &P_{no\, change} = \frac{1}{4} + \frac{1}{4} e^{-4 \beta d}
            + \frac{1}{2} e^{-2 (\alpha + \beta) d} \\
        &P_{transition} = \frac{1}{4} + \frac{1}{4} e^{-4 \beta d}
            - \frac{1}{2} e^{-2 (\alpha + \beta) d} \\
        &P_{transversion} = \frac{1}{4} - \frac{1}{4} e^{-4 \beta d}
        \end{aligned}

    Where :math:`\alpha` is transition rate, :math:`\beta` is transversion rate.

    Parameters
    ----------
    d : float
        distance or distances between sequences.
    alpha : float
        transition rate
    beta : float
        transversion rate
    seqtype : str
        String that holds type of the sequence. "DNA" (default) and "RNA"
        are valid options. Needed to assign proper nucleotide letters to
        PairwiseMatrix ids.

    Returns
    -------
    PairwiseMatrix
        Transition probability matrix. Rows are ancestral nucleotides,
        columns are descendant nucleotides.


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

    return _wrap_vector_tpm_function(
        _k2p, d, alpha, beta, allowed_seqtypes=("DNA", "RNA"), seqtype=seqtype
    )


def _k2p(
    d: npt.NDArray[np.float64],
    alpha: float,
    beta: float,
) -> npt.NDArray[np.float64]:
    e1 = np.exp(-4.0 * beta * d)
    e2 = np.exp(-2.0 * (alpha + beta) * d)

    same = 0.25 + 0.25 * e1 + 0.5 * e2
    transition = 0.25 + 0.25 * e1 - 0.5 * e2
    transversion = 0.25 - 0.25 * e1

    P = np.empty((len(d), 4, 4), dtype=np.float64)

    P[:] = transversion[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = same[:, None]

    # transitions: A<->G and C<->T
    P[:, 0, 1] = transition
    P[:, 1, 0] = transition
    P[:, 2, 3] = transition
    P[:, 3, 2] = transition

    return P


def f81(
    d: float,
    freqs: ArrayLike,
    seqtype: str = "DNA",
) -> SubstitutionMatrix:
    r"""
    Calculate the F81 transition probability matrix for a given distance.

    .. versionadded:: 0.7.3

    The Felsenstein 1981 (F81) model assumes equal substitution rates and allows
    differential base frequencies (:math:`\pi`). Transition probability for a
    nucleotide for sequences with distance :math:`d` between them can be
    calculated like:

    .. math::

        P_{ij} = \begin{cases}
            \pi_i + (1 - \pi_i) e^{-\frac{d}{b}} & i = j \\
            \pi_j (1 - e^{-\frac{d}{b}}) & i \neq j
        \end{cases}

    Where :math:`i,j \in \{A, C, G, T\}`, :math:`i` is ancestral nucleotide,
    :math:`j` is descendant nucleotide. Factor :math:`b` can be calculated as:

    .. math::
        b = 1 - \pi_A^2 - \pi_C^2 - \pi_G^2 - \pi_T^2

    Parameters
    ----------
    d : float
        distance or distances between sequences.
    freqs : array_like of float of shape (4,)
        Relative frequencies of nucleobases A, C, G, and T/U, respectively. Should sum
        to 1.
    seqtype : str
        String that holds type of the sequence. "DNA" (default) and "RNA"
        are valid options. Needed to assign proper nucleotide letters to
        SubstitutionMatrix dimensions.

    Returns
    -------
    SubstitutionMatrix
        Transition probability matrix. Rows are ancestral nucleotides,
        columns are descendant nucleotides.


    Notes
    -----
    The Felsenstein 1981 (F81) model was described in [1]_ in the context of maximum
    likelihood estimation.

    F81 is an extension of the JC69 model by allowing varying base frequencies. When
    the observed or user-provided based frequencies are equal (e.g., by specifying
    ``freqs=(.25, .25, .25, .25)``), the result will be identical to that of JC69.

    References
    ----------
    .. [1] Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
       likelihood approach. Journal of Molecular Evolution, 17(6), 368-376.

    """

    return _wrap_vector_tpm_function(
        _k2p,
        d,
        freqs=_check_freqs(freqs),
        allowed_seqtypes=("DNA", "RNA"),
        seqtype=seqtype,
    )


def _f81(
    d: float | npt.NDArray[np.float64],
    freqs: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """
    F81 transition probability matrix for branch length d.

    Parameters
    ----------
    d
        Branch length(s).
    freqs
        Stationary frequencies of states.
        Shape: (4,)
    """

    freqs = np.asarray(freqs, dtype=np.float64)

    b = 1.0 / (1.0 - np.sum(freqs**2))

    e = np.exp(-b * d)

    n = len(d)
    P = np.empty((n, 4, 4), dtype=np.float64)

    P[:] = freqs[None, None, :] * (1.0 - e)[:, None, None]

    idx = np.arange(4)
    P[:, idx, idx] = freqs[None, :] + (1.0 - freqs)[None, :] * e[:, None]

    return P
