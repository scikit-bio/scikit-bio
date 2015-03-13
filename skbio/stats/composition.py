r"""
Composition Statistics (:mod:`skbio.stats.composition`)
=======================================================

.. currentmodule:: skbio.stats.composition

This module provides functions for compositional data analysis.

Many 'omics datasets are inheriently compositional - meaning that they are best
interpreted as proportions or percentages rather than absolute counts.

Formally, :math:`x` is a composition if :math:`\sum_{i=0}^D x_{i} = c` and
:math:`x_{i} > 0`, :math:`1 \leq i \leq D`  and :math:`c` is a real valued
constant and there are :math:`D` components for each composition. In this
module :math:`c=1`. Compositional data can be analyzed using Aitchison
geometry [1]_

However, in this framework, standard real Euclidean operations such as addition
and multiplication no longer apply.  Only operations such as perturbation and
power can be used to manipulate this data [1]_

This module allows two styles of manipulation of compositional data.
Compositional data can be analyzed using perturbation and power operations,
which can be useful for simulation studies. The alternative strategy is to
transform compositional data into the real space.  Right now, the centre log
ratio transform (clr) [1]_ can be used to accomplish this.  This transform can
be useful for performing standard statistical tools such as parametric
hypothesis testing, regressions and more.

The major caveat of using this framework is dealing with zeros.
In the Aitchison geometry, only compositions with nonzero components can be
considered. The multiplicative replacement technique [2]_ can be used to
substitute these zeros with small pseudocounts without introducing major
distortions to the data.

Functions
---------

.. autosummary::
   :toctree: generated/

   multiplicative_replacement
   perturb
   perturb_inv
   power
   clr
   centralize

Examples
--------

>>> import numpy as np

Consider a very simple environment with only 3 species.
The species in the environment are equally distributed and their
proportions are equilvalent:

>>> otus = np.array([1./3, 1./3., 1./3])

Suppose that an antibiotic kills off half of the population for the
first two species, but doesn't harm the third species.
Then the perturbation vector would be as follows

>>> antibiotic = np.array([1./2, 1./2, 1])

And the resulting perturbation would be

>>> perturb(otus, antibiotic)
array([ 0.25,  0.25,  0.5 ])

Reference
---------
.. [1] V. Pawlowsky-Glahn. "Lecture Notes on Compositional Data Analysis"
.. [2] J. A. Martin-Fernandez. "Dealing With Zeros and Missing Values in
       Compositional Data Sets Using Nonparametric Imputation"

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
import scipy.stats as ss


def _closure(mat):
    """
    Performs closure to ensure that all elements add up to 1

    Parameters
    ----------
    mat : array_like
       a matrix of proportions where
       rows = compositions
       columns = components

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    """
    mat = np.atleast_2d(mat)
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    mat = mat / mat.sum(axis=1, keepdims=True)
    return mat.squeeze()


def multiplicative_replacement(mat, delta=None):
    """
    Performs multiplicative replacement strategy to replace
    all of the zeros with small non-zero values.  A closure
    operation is applied so that the compositions still
    add up to 1

    Parameters
    ----------
    mat: array_like
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
    >>> from skbio.stats.composition import multiplicative_replacement
    >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
    >>> multiplicative_replacement(X)
    array([[ 0.1875,  0.375 ,  0.375 ,  0.0625],
           [ 0.0625,  0.4375,  0.4375,  0.0625]])

    """
    mat = np.atleast_2d(mat)
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    z_mat = (mat == 0)

    num_samps, num_feats = mat.shape
    tot = z_mat.sum(axis=1, keepdims=True)

    if delta is None:
        delta = (1. / num_feats)**2

    zcnts = 1 - tot * delta
    mat = np.where(z_mat, delta, zcnts * mat)
    return mat.squeeze()


def perturb(x, y):
    r"""
    Performs the perturbation operation

    This operation is defined as
    :math:`x \oplus y = C[x_1 y_1, ..., x_D y_D]`

    :math:`C[x]` is the closure operation defined as
    :math:`C[x] = [\frac{x_1}{\sum x},...,\frac{x_D}{\sum x}]`
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

    Notes
    -----
    - Each row must add up to 1

    - All of the values in x and y must be greater than zero

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb
    >>> x = np.array([.1,.3,.4, .2])
    >>> y = np.array([1./6,1./6,1./3,1./3])
    >>> perturb(x,y)
    array([ 0.0625,  0.1875,  0.5   ,  0.25  ])

    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    return _closure(x * y)


def perturb_inv(x, y):
    r"""
    Performs the inverse perturbation operation

    This operation is defined as
    :math:`x \ominus y = C[x_1 y_1^{-1}, ..., x_D y_D^{-1}]`

    :math:`C[x]` is the closure operation defined as
    :math:`C[x] = [\frac{x_1}{\sum x},...,\frac{x_D}{\sum x}]`
    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : numpy.ndarray
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : numpy.ndarray
        rows = compositions and
        columns = components

    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Notes
    -----
    - Each row must add up to 1 for x.

    - y doesn't neccessary have to be a matrix of compositions

    - All of the values in x and y must be greater than zero

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb_inv
    >>> x = np.array([.1,.3,.4, .2])
    >>> y = np.array([1./6,1./6,1./3,1./3])
    >>> perturb_inv(x,y)
    array([ 0.14285714,  0.42857143,  0.28571429,  0.14285714])

    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    return _closure(x / y)


def power(x, a):
    r"""
    Performs the power operation

    This operation is defined as follows
    :math:`x \odot a = C[x_1^a, ..., x_D^a]`

    :math:`C[x]` is the closure operation defined as
    :math:`C[x] = [\frac{x_1}{\sum x},...,\frac{x_D}{\sum x}]`
    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.


    Parameters
    ----------
    x : numpy.ndarray, float
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

    Notes
    -----
    - Each row must add up to 1 for x

    - All of the values in x must be greater than zero

    >>> import numpy as np
    >>> from skbio.stats.composition import power
    >>> x = np.array([.1,.3,.4, .2])
    >>> power(x, .1)
    array([ 0.23059566,  0.25737316,  0.26488486,  0.24714631])

    """
    x = np.asarray(x, dtype=np.float64)
    return _closure(x**a)


def clr(mat):
    r"""
    Performs centre log ratio transformation that transforms
    compositions from Aitchison geometry to the real space.
    This transformation is an isometry, but not an isomorphism.

    This transformation is defined for a composition :math:`x` as follows

    :math:`clr(x) = ln[\frac{x_1}{g_m(x)}, ..., \frac{x_D}{g_m(x)}]`
    where :math:`g_m(x) = (\prod_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.

    Parameters
    ----------
    mat : numpy.ndarray, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    -------
    numpy.ndarray
         clr transformed matrix

    Notes
    -----
    - Each row must add up to 1

    - All of the values must be greater than zero for mat

    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1,.3,.4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])

    """
    lmat = np.log(np.atleast_2d(mat))
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")

    num_samps, num_feats = lmat.shape
    gm = lmat.mean(axis=1, keepdims=True)

    return (lmat - gm).squeeze()


def centralize(mat):
    """
    This calculates the average sample and centers the data
    around this sample.

    Parameters
    ----------
    mat : numpy.ndarray, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    -------
    numpy.ndarray
         centered composition matrix

    Notes
    -----
    - Each row must add up to 1 for mat

    - All of the values must be greater than zero

    >>> import numpy as np
    >>> from skbio.stats.composition import centralize
    >>> X = np.array([[.1,.3,.4, .2],[.2,.2,.2,.4]])
    >>> centralize(X)
    array([[ 0.17445763,  0.30216948,  0.34891526,  0.17445763],
           [ 0.32495488,  0.18761279,  0.16247744,  0.32495488]])

    """
    mat = np.atleast_2d(mat)
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    cen = ss.gmean(mat, axis=0)
    return perturb_inv(mat, cen)
