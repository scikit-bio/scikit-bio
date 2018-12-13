r"""
Composition Statistics (:mod:`skbio.stats.composition`)
=======================================================

.. currentmodule:: skbio.stats.composition

This module provides functions for compositional data analysis.

Many 'omics datasets are inherently compositional - meaning that they
are best interpreted as proportions or percentages rather than
absolute counts.

Formally, :math:`x` is a composition if :math:`\sum_{i=0}^D x_{i} = c`
and :math:`x_{i} > 0`, :math:`1 \leq i \leq D` and :math:`c` is a real
valued constant and there are :math:`D` components for each
composition. In this module :math:`c=1`. Compositional data can be
analyzed using Aitchison geometry. [1]_

However, in this framework, standard real Euclidean operations such as
addition and multiplication no longer apply. Only operations such as
perturbation and power can be used to manipulate this data.

This module allows two styles of manipulation of compositional data.
Compositional data can be analyzed using perturbation and power
operations, which can be useful for simulation studies. The
alternative strategy is to transform compositional data into the real
space.  Right now, the centre log ratio transform (clr) and
the isometric log ratio transform (ilr) [2]_ can be used to accomplish
this. This transform can be useful for performing standard statistical
tools such as parametric hypothesis testing, regressions and more.

The major caveat of using this framework is dealing with zeros.  In
the Aitchison geometry, only compositions with nonzero components can
be considered. The multiplicative replacement technique [3]_ can be
used to substitute these zeros with small pseudocounts without
introducing major distortions to the data.

Functions
---------

.. autosummary::
   :toctree:

   closure
   multiplicative_replacement
   perturb
   perturb_inv
   power
   inner
   clr
   clr_inv
   ilr
   ilr_inv
   alr
   alr_inv
   centralize
   ancom
   sbp_basis

References
----------
.. [1] V. Pawlowsky-Glahn, J. J. Egozcue, R. Tolosana-Delgado (2015),
   Modeling and Analysis of Compositional Data, Wiley, Chichester, UK

.. [2] J. J. Egozcue.,  "Isometric Logratio Transformations for
   Compositional Data Analysis" Mathematical Geology, 35.3 (2003)

.. [3] J. A. Martin-Fernandez,  "Dealing With Zeros and Missing Values in
   Compositional Data Sets Using Nonparametric Imputation",
   Mathematical Geology, 35.3 (2003)


Examples
--------

>>> import numpy as np

Consider a very simple environment with only 3 species. The species
in the environment are equally distributed and their proportions are
equivalent:

>>> otus = np.array([1./3, 1./3., 1./3])

Suppose that an antibiotic kills off half of the population for the
first two species, but doesn't harm the third species. Then the
perturbation vector would be as follows

>>> antibiotic = np.array([1./2, 1./2, 1])

And the resulting perturbation would be

>>> perturb(otus, antibiotic)
array([ 0.25,  0.25,  0.5 ])

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import scipy.stats
import skbio.util
from skbio.util._decorator import experimental


@experimental(as_of="0.4.0")
def closure(mat):
    """
    Performs closure to ensure that all elements add up to 1.

    Parameters
    ----------
    mat : array_like
       a matrix of proportions where
       rows = compositions
       columns = components

    Returns
    -------
    array_like, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Raises
    ------
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import closure
    >>> X = np.array([[2, 2, 6], [4, 4, 2]])
    >>> closure(X)
    array([[ 0.2,  0.2,  0.6],
           [ 0.4,  0.4,  0.2]])

    """
    mat = np.atleast_2d(mat)
    if np.any(mat < 0):
        raise ValueError("Cannot have negative proportions")
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    if np.all(mat == 0, axis=1).sum() > 0:
        raise ValueError("Input matrix cannot have rows with all zeros")
    mat = mat / mat.sum(axis=1, keepdims=True)
    return mat.squeeze()


@experimental(as_of="0.4.0")
def multiplicative_replacement(mat, delta=None):
    r"""Replace all zeros with small non-zero values

    It uses the multiplicative replacement strategy [1]_ ,
    replacing zeros with a small positive :math:`\delta`
    and ensuring that the compositions still add up to 1.


    Parameters
    ----------
    mat: array_like
       a matrix of proportions where
       rows = compositions and
       columns = components
    delta: float, optional
       a small number to be used to replace zeros
       If delta is not specified, then the default delta is
       :math:`\delta = \frac{1}{N^2}` where :math:`N`
       is the number of components

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Raises
    ------
    ValueError
       Raises an error if negative proportions are created due to a large
       `delta`.

    Notes
    -----
    This method will result in negative proportions if a large delta is chosen.

    References
    ----------
    .. [1] J. A. Martin-Fernandez. "Dealing With Zeros and Missing Values in
           Compositional Data Sets Using Nonparametric Imputation"


    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import multiplicative_replacement
    >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
    >>> multiplicative_replacement(X)
    array([[ 0.1875,  0.375 ,  0.375 ,  0.0625],
           [ 0.0625,  0.4375,  0.4375,  0.0625]])

    """
    mat = closure(mat)
    z_mat = (mat == 0)

    num_feats = mat.shape[-1]
    tot = z_mat.sum(axis=-1, keepdims=True)

    if delta is None:
        delta = (1. / num_feats)**2

    zcnts = 1 - tot * delta
    if np.any(zcnts) < 0:
        raise ValueError('The multiplicative replacement created negative '
                         'proportions. Consider using a smaller `delta`.')
    mat = np.where(z_mat, delta, zcnts * mat)
    return mat.squeeze()


@experimental(as_of="0.4.0")
def perturb(x, y):
    r"""
    Performs the perturbation operation.

    This operation is defined as

    .. math::
        x \oplus y = C[x_1 y_1, \ldots, x_D y_D]

    :math:`C[x]` is the closure operation defined as

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : array_like, float
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : array_like, float
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb
    >>> x = np.array([.1,.3,.4, .2])
    >>> y = np.array([1./6,1./6,1./3,1./3])
    >>> perturb(x,y)
    array([ 0.0625,  0.1875,  0.5   ,  0.25  ])

    """
    x, y = closure(x), closure(y)
    return closure(x * y)


@experimental(as_of="0.4.0")
def perturb_inv(x, y):
    r"""
    Performs the inverse perturbation operation.

    This operation is defined as

    .. math::
        x \ominus y = C[x_1 y_1^{-1}, \ldots, x_D y_D^{-1}]

    :math:`C[x]` is the closure operation defined as

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]


    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : array_like
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : array_like
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb_inv
    >>> x = np.array([.1,.3,.4, .2])
    >>> y = np.array([1./6,1./6,1./3,1./3])
    >>> perturb_inv(x,y)
    array([ 0.14285714,  0.42857143,  0.28571429,  0.14285714])
    """
    x, y = closure(x), closure(y)
    return closure(x / y)


@experimental(as_of="0.4.0")
def power(x, a):
    r"""
    Performs the power operation.

    This operation is defined as follows

    .. math::
        `x \odot a = C[x_1^a, \ldots, x_D^a]

    :math:`C[x]` is the closure operation defined as

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : array_like, float
        a matrix of proportions where
        rows = compositions and
        columns = components
    a : float
        a scalar float

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import power
    >>> x = np.array([.1,.3,.4, .2])
    >>> power(x, .1)
    array([ 0.23059566,  0.25737316,  0.26488486,  0.24714631])

    """
    x = closure(x)
    return closure(x**a).squeeze()


@experimental(as_of="0.4.0")
def inner(x, y):
    r"""
    Calculates the Aitchson inner product.

    This inner product is defined as follows

    .. math::
        \langle x, y \rangle_a =
        \frac{1}{2D} \sum\limits_{i=1}^{D} \sum\limits_{j=1}^{D}
        \ln\left(\frac{x_i}{x_j}\right) \ln\left(\frac{y_i}{y_j}\right)

    Parameters
    ----------
    x : array_like
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : array_like
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    numpy.ndarray
         inner product result

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import inner
    >>> x = np.array([.1, .3, .4, .2])
    >>> y = np.array([.2, .4, .2, .2])
    >>> inner(x, y)  # doctest: +ELLIPSIS
    0.2107852473...
    """
    x = closure(x)
    y = closure(y)
    a, b = clr(x), clr(y)
    return a.dot(b.T)


@experimental(as_of="0.4.0")
def clr(mat):
    r"""
    Performs centre log ratio transformation.

    This function transforms compositions from Aitchison geometry to
    the real space. The :math:`clr` transform is both an isometry and an
    isomorphism defined on the following spaces

    :math:`clr: S^D \rightarrow U`

    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`

    It is defined for a composition :math:`x` as follows:

    .. math::
        clr(x) = \ln\left[\frac{x_1}{g_m(x)}, \ldots, \frac{x_D}{g_m(x)}\right]

    where :math:`g_m(x) = (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.

    Parameters
    ----------
    mat : array_like, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    -------
    numpy.ndarray
         clr transformed matrix

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])

    """
    mat = closure(mat)
    lmat = np.log(mat)
    gm = lmat.mean(axis=-1, keepdims=True)
    return (lmat - gm).squeeze()


@experimental(as_of="0.4.0")
def clr_inv(mat):
    r"""
    Performs inverse centre log ratio transformation.

    This function transforms compositions from the real space to
    Aitchison geometry. The :math:`clr^{-1}` transform is both an isometry,
    and an isomorphism defined on the following spaces

    :math:`clr^{-1}: U \rightarrow S^D`

    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`

    This transformation is defined as follows

    .. math::
        clr^{-1}(x) = C[\exp( x_1, \ldots, x_D)]

    Parameters
    ----------
    mat : array_like, float
       a matrix of real values where
       rows = transformed compositions and
       columns = components

    Returns
    -------
    numpy.ndarray
         inverse clr transformed matrix

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr_inv
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr_inv(x)
    array([ 0.21383822,  0.26118259,  0.28865141,  0.23632778])

    """
    return closure(np.exp(mat))


@experimental(as_of="0.4.0")
def ilr(mat, basis=None, check=True):
    r"""
    Performs isometric log ratio transformation.

    This function transforms compositions from Aitchison simplex to
    the real space. The :math: ilr` transform is both an isometry,
    and an isomorphism defined on the following spaces

    :math:`ilr: S^D \rightarrow \mathbb{R}^{D-1}`

    The ilr transformation is defined as follows

    .. math::
        ilr(x) =
        [\langle x, e_1 \rangle_a, \ldots, \langle x, e_{D-1} \rangle_a]

    where :math:`[e_1,\ldots,e_{D-1}]` is an orthonormal basis in the simplex.

    If an orthornormal basis isn't specified, the J. J. Egozcue orthonormal
    basis derived from Gram-Schmidt orthogonalization will be used by
    default.

    Parameters
    ----------
    mat: numpy.ndarray
       a matrix of proportions where
       rows = compositions and
       columns = components

    basis: numpy.ndarray, float, optional
        orthonormal basis for Aitchison simplex
        defaults to J.J.Egozcue orthonormal basis.

    check: bool
        Specifies if the basis is orthonormal.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import ilr
    >>> x = np.array([.1, .3, .4, .2])
    >>> ilr(x)
    array([-0.7768362 , -0.68339802,  0.11704769])

    Notes
    -----
    If the `basis` parameter is specified, it is expected to be a basis in the
    Aitchison simplex.  If there are `D-1` elements specified in `mat`, then
    the dimensions of the basis needs be `D-1 x D`, where rows represent
    basis vectors, and the columns represent proportions.
    """
    mat = closure(mat)
    if basis is None:
        basis = clr_inv(_gram_schmidt_basis(mat.shape[-1]))
    else:
        if len(basis.shape) != 2:
            raise ValueError("Basis needs to be a 2D matrix, "
                             "not a %dD matrix." %
                             (len(basis.shape)))
        if check:
            _check_orthogonality(basis)

    return inner(mat, basis)


@experimental(as_of="0.4.0")
def ilr_inv(mat, basis=None, check=True):
    r"""
    Performs inverse isometric log ratio transform.

    This function transforms compositions from the real space to
    Aitchison geometry. The :math:`ilr^{-1}` transform is both an isometry,
    and an isomorphism defined on the following spaces

    :math:`ilr^{-1}: \mathbb{R}^{D-1} \rightarrow S^D`

    The inverse ilr transformation is defined as follows

    .. math::
        ilr^{-1}(x) = \bigoplus\limits_{i=1}^{D-1} x \odot e_i

    where :math:`[e_1,\ldots, e_{D-1}]` is an orthonormal basis in the simplex.

    If an orthonormal basis isn't specified, the J. J. Egozcue orthonormal
    basis derived from Gram-Schmidt orthogonalization will be used by
    default.


    Parameters
    ----------
    mat: numpy.ndarray, float
       a matrix of transformed proportions where
       rows = compositions and
       columns = components

    basis: numpy.ndarray, float, optional
        orthonormal basis for Aitchison simplex
        defaults to J.J.Egozcue orthonormal basis

    check: bool
        Specifies if the basis is orthonormal.


    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import ilr
    >>> x = np.array([.1, .3, .6,])
    >>> ilr_inv(x)
    array([ 0.34180297,  0.29672718,  0.22054469,  0.14092516])

    Notes
    -----
    If the `basis` parameter is specified, it is expected to be a basis in the
    Aitchison simplex.  If there are `D-1` elements specified in `mat`, then
    the dimensions of the basis needs be `D-1 x D`, where rows represent
    basis vectors, and the columns represent proportions.
    """

    if basis is None:
        basis = _gram_schmidt_basis(mat.shape[-1] + 1)
    else:
        if len(basis.shape) != 2:
            raise ValueError("Basis needs to be a 2D matrix, "
                             "not a %dD matrix." %
                             (len(basis.shape)))
        if check:
            _check_orthogonality(basis)
        # this is necessary, since the clr function
        # performs np.squeeze()
        basis = np.atleast_2d(clr(basis))

    return clr_inv(np.dot(mat, basis))


@experimental(as_of="0.5.5")
def alr(mat, denominator_idx=0):
    r"""
    Performs additive log ratio transformation.

    This function transforms compositions from a D-part Aitchison simplex to
    a non-isometric real space of D-1 dimensions. The argument
    `denominator_col` defines the index of the column used as the common
    denominator. The :math: `alr` transformed data are amenable to multivariate
    analysis as long as statistics don't involve distances.

    :math:`alr: S^D \rightarrow \mathbb{R}^{D-1}`

    The alr transformation is defined as follows

    .. math::
        alr(x) = \left[ \ln \frac{x_1}{x_D}, \ldots,
        \ln \frac{x_{D-1}}{x_D} \right]

    where :math:`D` is the index of the part used as common denominator.

    Parameters
    ----------
    mat: numpy.ndarray
       a matrix of proportions where
       rows = compositions and
       columns = components

    denominator_idx: int
       the index of the column (2D-matrix) or position (vector) of
       `mat` which should be used as the reference composition. By default
       `denominator_idx=0` to specify the first column or position.

    Returns
    -------
    numpy.ndarray
         alr-transformed data projected in a non-isometric real space
         of D-1 dimensions for a D-parts composition

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import alr
    >>> x = np.array([.1, .3, .4, .2])
    >>> alr(x)
    array([ 1.09861229,  1.38629436,  0.69314718])
    """
    mat = closure(mat)
    if mat.ndim == 2:
        mat_t = mat.T
        numerator_idx = list(range(0, mat_t.shape[0]))
        del numerator_idx[denominator_idx]
        lr = np.log(mat_t[numerator_idx, :]/mat_t[denominator_idx, :]).T
    elif mat.ndim == 1:
        numerator_idx = list(range(0, mat.shape[0]))
        del numerator_idx[denominator_idx]
        lr = np.log(mat[numerator_idx]/mat[denominator_idx])
    else:
        raise ValueError("mat must be either 1D or 2D")
    return lr


@experimental(as_of="0.5.5")
def alr_inv(mat, denominator_idx=0):
    r"""
    Performs inverse additive log ratio transform.

    This function transforms compositions from the non-isometric real space of
    alrs to Aitchison geometry.

    :math:`alr^{-1}: \mathbb{R}^{D-1} \rightarrow S^D`

    The inverse alr transformation is defined as follows

    .. math::
         alr^{-1}(x) = C[exp([y_1, y_2, ..., y_{D-1}, 0])]

    where :math:`C[x]` is the closure operation defined as

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    mat: numpy.ndarray
       a matrix of alr-transformed data
    denominator_idx: int
       the index of the column (2D-composition) or position (1D-composition) of
       the output where the common denominator should be placed. By default
       `denominator_idx=0` to specify the first column or position.

    Returns
    -------
    numpy.ndarray
         Inverse alr transformed matrix or vector where rows sum to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import alr, alr_inv
    >>> x = np.array([.1, .3, .4, .2])
    >>> alr_inv(alr(x))
    array([ 0.1,  0.3,  0.4,  0.2])
    """
    mat = np.array(mat)
    if mat.ndim == 2:
        mat_idx = np.insert(mat, denominator_idx,
                            np.repeat(0, mat.shape[0]), axis=1)
        comp = np.zeros(mat_idx.shape)
        comp[:, denominator_idx] = 1 / (np.exp(mat).sum(axis=1) + 1)
        numerator_idx = list(range(0, comp.shape[1]))
        del numerator_idx[denominator_idx]
        for i in numerator_idx:
            comp[:, i] = comp[:, denominator_idx] * np.exp(mat_idx[:, i])
    elif mat.ndim == 1:
        mat_idx = np.insert(mat, denominator_idx, 0, axis=0)
        comp = np.zeros(mat_idx.shape)
        comp[denominator_idx] = 1 / (np.exp(mat).sum(axis=0) + 1)
        numerator_idx = list(range(0, comp.shape[0]))
        del numerator_idx[denominator_idx]
        for i in numerator_idx:
            comp[i] = comp[denominator_idx] * np.exp(mat_idx[i])
    else:
        raise ValueError("mat must be either 1D or 2D")
    return comp


@experimental(as_of="0.4.0")
def centralize(mat):
    r"""Center data around its geometric average.

    Parameters
    ----------
    mat : array_like, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    -------
    numpy.ndarray
         centered composition matrix

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import centralize
    >>> X = np.array([[.1,.3,.4, .2],[.2,.2,.2,.4]])
    >>> centralize(X)
    array([[ 0.17445763,  0.30216948,  0.34891526,  0.17445763],
           [ 0.32495488,  0.18761279,  0.16247744,  0.32495488]])

    """
    mat = closure(mat)
    cen = scipy.stats.gmean(mat, axis=0)
    return perturb_inv(mat, cen)


@experimental(as_of="0.4.1")
def ancom(table, grouping,
          alpha=0.05,
          tau=0.02,
          theta=0.1,
          multiple_comparisons_correction='holm-bonferroni',
          significance_test=None,
          percentiles=(0.0, 25.0, 50.0, 75.0, 100.0)):
    r""" Performs a differential abundance test using ANCOM.

    This is done by calculating pairwise log ratios between all features
    and performing a significance test to determine if there is a significant
    difference in feature ratios with respect to the variable of interest.

    In an experiment with only two treatments, this tests the following
    hypothesis for feature :math:`i`

    .. math::

        H_{0i}: \mathbb{E}[\ln(u_i^{(1)})] = \mathbb{E}[\ln(u_i^{(2)})]

    where :math:`u_i^{(1)}` is the mean abundance for feature :math:`i` in the
    first group and :math:`u_i^{(2)}` is the mean abundance for feature
    :math:`i` in the second group.

    Parameters
    ----------
    table : pd.DataFrame
        A 2D matrix of strictly positive values (i.e. counts or proportions)
        where the rows correspond to samples and the columns correspond to
        features.
    grouping : pd.Series
        Vector indicating the assignment of samples to groups.  For example,
        these could be strings or integers denoting which group a sample
        belongs to.  It must be the same length as the samples in `table`.
        The index must be the same on `table` and `grouping` but need not be
        in the same order.
    alpha : float, optional
        Significance level for each of the statistical tests.
        This can can be anywhere between 0 and 1 exclusive.
    tau : float, optional
        A constant used to determine an appropriate cutoff.
        A value close to zero indicates a conservative cutoff.
        This can can be anywhere between 0 and 1 exclusive.
    theta : float, optional
        Lower bound for the proportion for the W-statistic.
        If all W-statistics are lower than theta, then no features
        will be detected to be differentially significant.
        This can can be anywhere between 0 and 1 exclusive.
    multiple_comparisons_correction : {None, 'holm-bonferroni'}, optional
        The multiple comparison correction procedure to run.  If None,
        then no multiple comparison correction procedure will be run.
        If 'holm-boniferroni' is specified, then the Holm-Boniferroni
        procedure [1]_ will be run.
    significance_test : function, optional
        A statistical significance function to test for significance between
        classes.  This function must be able to accept at least two 1D
        array_like arguments of floats and returns a test statistic and a
        p-value. By default ``scipy.stats.f_oneway`` is used.
    percentiles : iterable of floats, optional
        Percentile abundances to return for each feature in each group. By
        default, will return the minimum, 25th percentile, median, 75th
        percentile, and maximum abundances for each feature in each group.

    Returns
    -------
    pd.DataFrame
        A table of features, their W-statistics and whether the null hypothesis
        is rejected.

        `"W"` is the W-statistic, or number of features that a single feature
        is tested to be significantly different against.

        `"Reject null hypothesis"` indicates if feature is differentially
        abundant across groups (`True`) or not (`False`).

    pd.DataFrame
        A table of features and their percentile abundances in each group. If
        ``percentiles`` is empty, this will be an empty ``pd.DataFrame``. The
        rows in this object will be features, and the columns will be a
        multi-index where the first index is the percentile, and the second
        index is the group.

    See Also
    --------
    multiplicative_replacement
    scipy.stats.ttest_ind
    scipy.stats.f_oneway
    scipy.stats.wilcoxon
    scipy.stats.kruskal

    Notes
    -----
    The developers of this method recommend the following significance tests
    ([2]_, Supplementary File 1, top of page 11): if there are 2 groups, use
    the standard parametric t-test (``scipy.stats.ttest_ind``) or
    non-parametric Wilcoxon rank sum test (``scipy.stats.wilcoxon``).
    If there are more than 2 groups, use parametric one-way ANOVA
    (``scipy.stats.f_oneway``) or nonparametric Kruskal-Wallis
    (``scipy.stats.kruskal``). Because one-way ANOVA is equivalent
    to the standard t-test when the number of groups is two, we default to
    ``scipy.stats.f_oneway`` here, which can be used when there are two or
    more groups.  Users should refer to the documentation of these tests in
    SciPy to understand the assumptions made by each test.

    This method cannot handle any zero counts as input, since the logarithm
    of zero cannot be computed.  While this is an unsolved problem, many
    studies, including [2]_, have shown promising results by adding
    pseudocounts to all values in the matrix. In [2]_, a pseudocount of 0.001
    was used, though the authors note that a pseudocount of 1.0 may also be
    useful. Zero counts can also be addressed using the
    ``multiplicative_replacement`` method.

    References
    ----------
    .. [1] Holm, S. "A simple sequentially rejective multiple test procedure".
       Scandinavian Journal of Statistics (1979), 6.
    .. [2] Mandal et al. "Analysis of composition of microbiomes: a novel
       method for studying microbial composition", Microbial Ecology in Health
       & Disease, (2015), 26.

    Examples
    --------
    First import all of the necessary modules:

    >>> from skbio.stats.composition import ancom
    >>> import pandas as pd

    Now let's load in a DataFrame with 6 samples and 7 features (e.g.,
    these may be bacterial OTUs):

    >>> table = pd.DataFrame([[12, 11, 10, 10, 10, 10, 10],
    ...                       [9,  11, 12, 10, 10, 10, 10],
    ...                       [1,  11, 10, 11, 10, 5,  9],
    ...                       [22, 21, 9,  10, 10, 10, 10],
    ...                       [20, 22, 10, 10, 13, 10, 10],
    ...                       [23, 21, 14, 10, 10, 10, 10]],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'],
    ...                      columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6',
    ...                               'b7'])

    Then create a grouping vector. In this example, there is a treatment group
    and a placebo group.

    >>> grouping = pd.Series(['treatment', 'treatment', 'treatment',
    ...                       'placebo', 'placebo', 'placebo'],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'])

    Now run ``ancom`` to determine if there are any features that are
    significantly different in abundance between the treatment and the placebo
    groups. The first DataFrame that is returned contains the ANCOM test
    results, and the second contains the percentile abundance data for each
    feature in each group.

    >>> ancom_df, percentile_df = ancom(table, grouping)
    >>> ancom_df['W']
    b1    0
    b2    4
    b3    0
    b4    1
    b5    1
    b6    0
    b7    1
    Name: W, dtype: int64

    The W-statistic is the number of features that a single feature is tested
    to be significantly different against.  In this scenario, `b2` was detected
    to have significantly different abundances compared to four of the other
    features. To summarize the results from the W-statistic, let's take a look
    at the results from the hypothesis test. The `Reject null hypothesis`
    column in the table indicates whether the null hypothesis was rejected,
    and that a feature was therefore observed to be differentially abundant
    across the groups.

    >>> ancom_df['Reject null hypothesis']
    b1    False
    b2     True
    b3    False
    b4    False
    b5    False
    b6    False
    b7    False
    Name: Reject null hypothesis, dtype: bool

    From this we can conclude that only `b2` was significantly different in
    abundance between the treatment and the placebo. We still don't know, for
    example, in which group `b2` was more abundant. We therefore may next be
    interested in comparing the abundance of `b2` across the two groups.
    We can do that using the second DataFrame that was returned. Here we
    compare the median (50th percentile) abundance of `b2` in the treatment and
    placebo groups:

    >>> percentile_df[50.0].loc['b2']
    Group
    placebo      21.0
    treatment    11.0
    Name: b2, dtype: float64

    We can also look at a full five-number summary for ``b2`` in the treatment
    and placebo groups:

    >>> percentile_df.loc['b2'] # doctest: +NORMALIZE_WHITESPACE
    Percentile  Group
    0.0         placebo      21.0
    25.0        placebo      21.0
    50.0        placebo      21.0
    75.0        placebo      21.5
    100.0       placebo      22.0
    0.0         treatment    11.0
    25.0        treatment    11.0
    50.0        treatment    11.0
    75.0        treatment    11.0
    100.0       treatment    11.0
    Name: b2, dtype: float64

    Taken together, these data tell us that `b2` is present in significantly
    higher abundance in the placebo group samples than in the treatment group
    samples.

    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError('`table` must be a `pd.DataFrame`, '
                        'not %r.' % type(table).__name__)
    if not isinstance(grouping, pd.Series):
        raise TypeError('`grouping` must be a `pd.Series`,'
                        ' not %r.' % type(grouping).__name__)

    if np.any(table <= 0):
        raise ValueError('Cannot handle zeros or negative values in `table`. '
                         'Use pseudocounts or ``multiplicative_replacement``.'
                         )

    if not 0 < alpha < 1:
        raise ValueError('`alpha`=%f is not within 0 and 1.' % alpha)

    if not 0 < tau < 1:
        raise ValueError('`tau`=%f is not within 0 and 1.' % tau)

    if not 0 < theta < 1:
        raise ValueError('`theta`=%f is not within 0 and 1.' % theta)

    if multiple_comparisons_correction is not None:
        if multiple_comparisons_correction != 'holm-bonferroni':
            raise ValueError('%r is not an available option for '
                             '`multiple_comparisons_correction`.'
                             % multiple_comparisons_correction)

    if (grouping.isnull()).any():
        raise ValueError('Cannot handle missing values in `grouping`.')

    if (table.isnull()).any().any():
        raise ValueError('Cannot handle missing values in `table`.')

    percentiles = list(percentiles)
    for percentile in percentiles:
        if not 0.0 <= percentile <= 100.0:
            raise ValueError('Percentiles must be in the range [0, 100], %r '
                             'was provided.' % percentile)

    duplicates = skbio.util.find_duplicates(percentiles)
    if duplicates:
        formatted_duplicates = ', '.join(repr(e) for e in duplicates)
        raise ValueError('Percentile values must be unique. The following'
                         ' value(s) were duplicated: %s.' %
                         formatted_duplicates)

    groups = np.unique(grouping)
    num_groups = len(groups)

    if num_groups == len(grouping):
        raise ValueError(
            "All values in `grouping` are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' variance because each group of samples "
            "contains only a single sample).")

    if num_groups == 1:
        raise ValueError(
            "All values the `grouping` are the same. This method cannot "
            "operate on a grouping vector with only a single group of samples"
            "(e.g., there are no 'between' variance because there is only a "
            "single group).")

    if significance_test is None:
        significance_test = scipy.stats.f_oneway

    table_index_len = len(table.index)
    grouping_index_len = len(grouping.index)
    mat, cats = table.align(grouping, axis=0, join='inner')
    if (len(mat) != table_index_len or len(cats) != grouping_index_len):
        raise ValueError('`table` index and `grouping` '
                         'index must be consistent.')

    n_feat = mat.shape[1]

    _logratio_mat = _log_compare(mat.values, cats.values, significance_test)
    logratio_mat = _logratio_mat + _logratio_mat.T

    # Multiple comparisons
    if multiple_comparisons_correction == 'holm-bonferroni':
        logratio_mat = np.apply_along_axis(_holm_bonferroni,
                                           1, logratio_mat)
    np.fill_diagonal(logratio_mat, 1)
    W = (logratio_mat < alpha).sum(axis=1)
    c_start = W.max() / n_feat
    if c_start < theta:
        reject = np.zeros_like(W, dtype=bool)
    else:
        # Select appropriate cutoff
        cutoff = c_start - np.linspace(0.05, 0.25, 5)
        prop_cut = np.array([(W > n_feat*cut).mean() for cut in cutoff])
        dels = np.abs(prop_cut - np.roll(prop_cut, -1))
        dels[-1] = 0

        if (dels[0] < tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[1]
        elif (dels[0] >= tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[2]
        elif (dels[1] >= tau) and (dels[2] < tau) and (dels[3] < tau):
            nu = cutoff[3]
        else:
            nu = cutoff[4]
        reject = (W >= nu*n_feat)

    feat_ids = mat.columns
    ancom_df = pd.DataFrame(
        {'W': pd.Series(W, index=feat_ids),
         'Reject null hypothesis': pd.Series(reject, index=feat_ids)})

    if len(percentiles) == 0:
        return ancom_df, pd.DataFrame()
    else:
        data = []
        columns = []
        for group in groups:
            feat_dists = mat[cats == group]
            for percentile in percentiles:
                columns.append((percentile, group))
                data.append(np.percentile(feat_dists, percentile, axis=0))
        columns = pd.MultiIndex.from_tuples(columns,
                                            names=['Percentile', 'Group'])
        percentile_df = pd.DataFrame(
            np.asarray(data).T, columns=columns, index=feat_ids)
        return ancom_df, percentile_df


def _holm_bonferroni(p):
    """ Performs Holm-Bonferroni correction for pvalues
    to account for multiple comparisons

    Parameters
    ---------
    p: numpy.array
        array of pvalues

    Returns
    -------
    numpy.array
        corrected pvalues
    """
    K = len(p)
    sort_index = -np.ones(K, dtype=np.int64)
    sorted_p = np.sort(p)
    sorted_p_adj = sorted_p*(K-np.arange(K))
    for j in range(K):
        idx = (p == sorted_p[j]) & (sort_index < 0)
        num_ties = len(sort_index[idx])
        sort_index[idx] = np.arange(j, (j+num_ties), dtype=np.int64)

    sorted_holm_p = [min([max(sorted_p_adj[:k]), 1])
                     for k in range(1, K+1)]
    holm_p = [sorted_holm_p[sort_index[k]] for k in range(K)]
    return holm_p


def _log_compare(mat, cats,
                 significance_test=scipy.stats.ttest_ind):
    """ Calculates pairwise log ratios between all features and performs a
    significiance test (i.e. t-test) to determine if there is a significant
    difference in feature ratios with respect to the variable of interest.

    Parameters
    ----------
    mat: np.array
       rows correspond to samples and columns correspond to
       features (i.e. OTUs)
    cats: np.array, float
       Vector of categories
    significance_test: function
        statistical test to run

    Returns:
    --------
    log_ratio : np.array
        log ratio pvalue matrix
    """
    r, c = mat.shape
    log_ratio = np.zeros((c, c))
    log_mat = np.log(mat)
    cs = np.unique(cats)

    def func(x):
        return significance_test(*[x[cats == k] for k in cs])

    for i in range(c-1):
        ratio = (log_mat[:, i].T - log_mat[:, i+1:].T).T
        m, p = np.apply_along_axis(func,
                                   axis=0,
                                   arr=ratio)
        log_ratio[i, i+1:] = np.squeeze(np.array(p.T))
    return log_ratio


def _gram_schmidt_basis(n):
    """
    Builds clr transformed basis derived from
    gram schmidt orthogonalization

    Parameters
    ----------
    n : int
        Dimension of the Aitchison simplex
    """
    basis = np.zeros((n, n-1))
    for j in range(n-1):
        i = j + 1
        e = np.array([(1/i)]*i + [-1] +
                     [0]*(n-i-1))*np.sqrt(i/(i+1))
        basis[:, j] = e
    return basis.T


@experimental(as_of="0.5.5")
def sbp_basis(sbp):
    r"""
    Builds an orthogonal basis from a sequential binary partition (SBP). As
    explained in [1]_, the SBP is a hierarchical collection of binary
    divisions of compositional parts. The child groups are divided again until
    all groups contain a single part. The SBP can be encoded in a
    :math:`(D - 1) \times D` matrix where, for each row, parts can be grouped
    by -1 and +1 tags, and 0 for excluded parts. The `sbp_basis` method was
    originally derived from function `gsi.buildilrBase()` found in the R
    package `compositions` [2]_. The ith balance is computed as follows

    .. math::
        b_i = \sqrt{ \frac{r_i s_i}{r_i+s_i} }
        \ln \left( \frac{g(x_{r_i})}{g(x_{s_i})} \right)

    where :math:`b_i` is the ith balance corresponding to the ith row in the
    SBP, :math:`r_i` and :math:`s_i` and the number of respectively `+1` and
    `-1` labels in the ith row of the SBP and where :math:`g(x) =
    (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric mean of :math:`x`.

    Parameters
    ----------
    sbp: np.array, int
        A contrast matrix, also known as a sequential binary partition, where
        every row represents a partition between two groups of features. A part
        labelled `+1` would correspond to that feature being in the numerator
        of the given row partition, a part labelled `-1` would correspond to
        features being in the denominator of that given row partition, and `0`
        would correspond to features excluded in the row partition.

    Returns
    -------
    numpy.array
        An orthonormal basis in the Aitchison simplex

    Examples
    --------
    >>> import numpy as np
    >>> sbp = np.array([[1, 1,-1,-1,-1],
    ...                 [1,-1, 0, 0, 0],
    ...                 [0, 0, 1,-1,-1],
    ...                 [0, 0, 0, 1,-1]])
    ...
    >>> sbp_basis(sbp)
    array([[ 0.31209907,  0.31209907,  0.12526729,  0.12526729,  0.12526729],
           [ 0.36733337,  0.08930489,  0.18112058,  0.18112058,  0.18112058],
           [ 0.17882092,  0.17882092,  0.40459293,  0.11888261,  0.11888261],
           [ 0.18112058,  0.18112058,  0.18112058,  0.36733337,  0.08930489]])

    References
    ----------
    .. [1] Parent, S.É., Parent, L.E., Egozcue, J.J., Rozane, D.E.,
       Hernandes, A., Lapointe, L., Hébert-Gentile, V., Naess, K.,
       Marchand, S., Lafond, J., Mattos, D., Barlow, P., Natale, W., 2013.
       The plant ionome revisited by the nutrient balance concept.
       Front. Plant Sci. 4, 39, http://dx.doi.org/10.3389/fpls.2013.00039.
    .. [2] van den Boogaart, K. Gerald, Tolosana-Delgado, Raimon and Bren,
       Matevz, 2014. `compositions`: Compositional Data Analysis. R package
       version 1.40-1. https://CRAN.R-project.org/package=compositions.
    """

    n_pos = (sbp == 1).sum(axis=1)
    n_neg = (sbp == -1).sum(axis=1)
    psi = np.zeros(sbp.shape)
    for i in range(0, sbp.shape[0]):
        psi[i, :] = sbp[i, :] * np.sqrt((n_neg[i] / n_pos[i])**sbp[i, :] /
                                        np.sum(np.abs(sbp[i, :])))
    return clr_inv(psi)


def _check_orthogonality(basis):
    """
    Checks to see if basis is truly orthonormal in the
    Aitchison simplex

    Parameters
    ----------
    basis: numpy.ndarray
        basis in the Aitchison simplex
    """
    basis = np.atleast_2d(basis)
    if not np.allclose(inner(basis, basis), np.identity(len(basis)),
                       rtol=1e-4, atol=1e-6):
        raise ValueError("Aitchison basis is not orthonormal")
