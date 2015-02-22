r"""
Composition Statistics (:mod:`skbio.stats.composition`)
===============================================

.. currentmodule:: skbio.stats.composition

This module provides functions for compositional data analysis.

Functions
---------

.. autosummary::
   :toctree: generated/

   zero_imputation
   perturb
   perturb_inv
   power
   clr
   centralize
   

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
    Performs _closure to ensure that all elements add up to 1
    mat: numpy.ndarray
       columns = features
       rows = samples
    """
    if mat.dtype != type(0.0):
        mat = mat.astype(np.float64)

    if len(mat.shape) == 1:
        num_samps = len(mat)
        total = mat.sum()
    else:
        num_samps, num_feats = mat.shape
        total = np.reshape(mat.sum(axis=1), (num_samps, 1))
    return np.divide(mat, total)


def multiplicative_replacement(mat, delta=None):
    """
    Performs multiplicative replacement strategy
    mat: numpy.ndarray
       columns = features
       rows = samples
    Returns:
    --------
    mat: numpy.ndarray
    """
    num_samps, num_feats = mat.shape
    if delta==None:
        delta = (1. / num_feats)**2
    z_mat = (mat == 0).astype(np.float32)
    zcnts = 1 - np.reshape(z_mat.sum(axis=1) * delta, (num_samps, 1) )
    mat = z_mat*delta + np.multiply((1-z_mat), np.multiply(zcnts,mat))
    return _closure(mat)

def perturb(x, y):
    """
    Performs the perturbation operation
    x: numpy.ndarray
    y: numpy.ndarray
    """    
    mat = np.multiply(x, y)
    return _closure(mat)

def perturb_inv(x, y):
    """
    Performs the inverse perturbation operation
    x: numpy.ndarray
    y: numpy.ndarray
    """
    _y = power(y,-1)
    mat = np.multiply(x, _y)
    return _closure(mat)

def power(x, y):
    """
    Performs the perturbation operation
    x: numpy.ndarray
    y: numpy.ndarray
    """
    mat = np.multiply(np.log(x), y)
    return _closure(np.exp(mat))

def clr(mat):
    """
    Performs centre log ratio transformation
    Returns
    =======
    clr: numpy.ndarray
    
    """
    lmat = np.log(mat) # Need to account for zeros
    if len(mat.shape) == 1:
        num_samps = len(mat)
        gm = lmat.mean()
    else:
        num_samps, num_feats = mat.shape
        gm = lmat.mean(axis = 1)
        gm = np.reshape(gm, (num_samps, 1))
    _clr = lmat - gm
    return _clr
    
def centralize(mat):
    """
    This calculates the average sample and centers the data
    around this sample.

    mat: numpy.ndarray
       A contingency table of normalized proportions
       columns = features
       rows = samples

    Notes
    -----
    - Each row must add up to 1
    
    - All of the values must be greater than zero
    
    Reference
    ---------
    .. [1] J.J. Ezogcue http://www.sediment.uni-goettingen.de/staff/tolosana/extra/CoDa.pdf

    """
    assert len(mat.shape)>1, "Need more than 1 row"
    r,c = mat.shape
    cen = ss.gmean(mat, axis=0)
    cen = np.tile(cen, (r,1))
    centered = perturb_inv(mat, cen)
    return centered

