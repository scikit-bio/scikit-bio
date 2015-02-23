r"""
Composition Statistics (:mod:`skbio.stats.composition`)
=======================================================

.. currentmodule:: skbio.stats.composition

This module provides functions for compositional data analysis.

Many 'omics datasets are inheriently compositional - meaning that they are best
interpreted as proportions or percentages rather than absolute counts.

Formally, :math:`x` is a composition if :math:`\sum_{i=0}^D x_{i} = c` and
:math:`x_{i} > 0`, :math:`1 \leq i \leq D`  and :math:`c` is a real valued constant
and there are :math:`D` components for each composition. In this module
:math:`c=1`. Compositional data can be analyzed using Aitchison geometry [1]_

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

Reference
---------
.. [1] V. Pawlowsky-Glahn. "Lecture Notes on Compositional Data Analysis"
.. [2] J. A. Martin-Fernandez. "Dealing With Zeros and Missing Values in
       Compositional Data Sets Using Nonparametric Imputation"

"""
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

import scipy.stats as ss


def _closure(mat):
    """
    Performs closure to ensure that all elements add up to 1
    mat : numpy.ndarray
       a matrix of proportions where
       rows = compositions
       columns = components

    Returns
    -------
    mat_ : numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    """
    if not isinstance(mat.dtype, float):
        mat = mat.astype(np.float64)

    if len(mat.shape) == 1:
        num_samps = len(mat)
        total = mat.sum()
    else:
        num_samps, num_feats = mat.shape
        total = np.reshape(mat.sum(axis=1), (num_samps, 1))
    mat_ = np.divide(mat, total)
    return mat_


def multiplicative_replacement(mat, delta=None):
    """
    Performs multiplicative replacement strategy to replace
    all of the zeros with small non-zero values.  A closure
    operation is applied so that the compositions still
    add up to 1

    mat: numpy.ndarray
       a matrix of proportions where
       rows = compositions and 
       columns = components

    Returns:
    --------
    mat_ : numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1
    """
    num_samps, num_feats = mat.shape
    if delta is None:
        delta = (1. / num_feats)**2
    z_mat = (mat == 0).astype(np.float32)
    zcnts = 1 - np.reshape(z_mat.sum(axis=1) * delta, (num_samps, 1))
    mat_ = _closure(z_mat*delta + np.multiply((1-z_mat),
                                              np.multiply(zcnts, mat)))
    return mat_


def perturb(x, y):
    r"""
    Performs the perturbation operation

    This operation is defined as
    :math:`x \oplus y = C[x_1 y_1, ..., x_D y_D]`
    where :math:`C[x]` is the closure operation on the composition
    :math:`x` and :math:`D` is the number of components for every
    composition.

    Parameters
    ----------
    x : numpy.ndarray, float
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : numpy.ndarray, float
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    mat_ : numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Notes
    -----
    - Each row must add up to 1

    - All of the values in x and y must be greater than zero

    """
    mat_ = _closure(np.multiply(x, y))
    return mat_


def perturb_inv(x, y):
    r"""
    Performs the inverse perturbation operation

    This operation is defined as
    :math:`x \ominus y = C[x_1 y_1^{-1}, ..., x_D y_D^{-1}]`
    where :math:`C[x]` is the closure operation on the composition
    :math:`x` and :math:`D` is the number of components for every
    composition.

    Parameters
    ----------
    x : numpy.ndarray
        a matrix of proportions where
        rows = compositions and
        columns = components
    y : numpy.ndarray
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    mat_ : numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Notes
    -----
    - Each row must add up to 1

    - All of the values in x and y must be greater than zero
    """
    _y = power(y, -1)
    mat = np.multiply(x, _y)
    return _closure(mat)


def power(x, a):
    r"""
    Performs the power operation

    This operation is defined as follows
    :math:`x \odot y = C[x_1^a, ..., x_D^a]`
    Where :math:`C[x]` is the closure operation on the composition
    :math:`x` and :math:`D` is the number of components for every
    composition.

    Parameters
    ----------
    x : numpy.ndarray, float
        a matrix of proportions where
        rows = compositions and
        columns = components
    a : numpy.ndarray, float
        a matrix of proportions where
        rows = compositions and
        columns = components

    Returns
    -------
    mat_ : numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Notes
    -----
    - Each row must add up to 1

    - All of the values in x must be greater than zero

    """
    mat = np.multiply(np.log(x), a)
    return _closure(np.exp(mat))


def clr(mat):
    r"""
    Performs centre log ratio transformation that transforms
    compositions from Aitchison geometry to the real space.
    This transformation is an isometry, but not an isomorphism.

    This transformation is defined for a composition :math:`x` as follows

    :math:`clr(x) = ln[\frac{x_1}{g_m(x)}, ..., \frac{x_D}{g_m(x)}]`
    where :math:`g_m(x) = (\prod_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.

    mat : numpy.ndarray, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    =======
    _clr : numpy.ndarray
         clr transformed matrix

    Notes
    -----
    - Each row must add up to 1

    - All of the values must be greater than zero

    """
    lmat = np.log(mat)
    if len(mat.shape) == 1:
        num_samps = len(mat)
        gm = lmat.mean()
    else:
        num_samps, num_feats = mat.shape
        gm = lmat.mean(axis=1)
        gm = np.reshape(gm, (num_samps, 1))
    _clr = lmat - gm
    return _clr


def centralize(mat):
    """
    This calculates the average sample and centers the data
    around this sample.

    mat : numpy.ndarray, float
       a matrix of proportions where
       rows = compositions and
       columns = components

    Returns
    =======
    centered : numpy.ndarray
         centered proportion matrix

    Notes
    -----
    - Each row must add up to 1

    - All of the values must be greater than zero

    """
    assert len(mat.shape) > 1, "Need more than 1 row"
    r, c = mat.shape
    cen = ss.gmean(mat, axis=0)
    cen = np.tile(cen, (r, 1))
    centered = perturb_inv(mat, cen)
    return centered
