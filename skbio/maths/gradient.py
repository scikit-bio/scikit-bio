#!/usr/bin/env python
r"""
Gradient analyses (:mod:`skbio.maths.gradient`)
===============================================

.. currentmodule:: skbio.maths.gradient

This module provides functionality for performing gradient analyses.

Functions
---------

.. autosummary::
   :toctree: generated/

   gradient_analysis

"""
from __future__ import division

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from copy import deepcopy
import numpy as np
from collections import namedtuple

from skbio.util.sort import signed_natsort


def make_groups(ord_res, metamap, vector_category, sort_category=None):
    """Groups the sample ids in metamap by the values in vector_category
    Parameters
    ----------
    ord_res : skbio.maths.stats.ordination.OrdinationResults
        The ordination results namedtuple
    metamap : pandas.DataFrame
        The metadata map where index are samples ids and columns are the
        metadata categories
    vector_category : str
        The category from metamap to use to create the groups
    sort_category : str, optional
        The category from metamap to use to sort groups

    Returns
    -------
    dictionary
        A dictionary in which the keys represent the group label and values are
        ordered lists of (sort_value, sample id) tuples
    """
    # If sort_category is provided, we used the value of such category to sort
    # otherwise we use the sample id
    if sort_category:
        sort_val = lambda sid: metamap[sid][sort_category]
    else:
        sort_val = lambda sid: sid

    # Group by vector_category
    metamap_t = metamap.transpose()
    gb = metamap_t.groupby(vector_category)
    groups = {}
    for g, df in gb:
        groups[g] = signed_natsort([(sort_val(sid), sid) for sid in df.index])

    return groups


def weight_by_vector(vector, w_vector):
    """weights the values of 'vector' given a weighting vector 'w_vector'.

    Each value in 'vector' will be weighted by the 'rate of change'
    to 'optimal rate of change' ratio, meaning that when calling this function
    over evenly spaced 'w_vector' values, no change will be reflected on the
    output.

    Parameters
    ----------
    vector: numpy array
        Values to weight
    w_vector: numpy array
        Values used to weight 'vector'

    Returns
    -------
    numpy array
        A weighted version of 'vector'.

    Raises
    ------
    ValueError
        If vector and w_vector don't have equal lengths
        If w_vector is not a gradient
    TypeError
        If vector and w_vector are not iterables
    """
    try:
        if len(vector) != len(w_vector):
            raise ValueError("vector (%d) & w_vector (%d) must be equal "
                             "lengths" % (len(vector), len(w_vector)))
    except TypeError:
        raise TypeError("vector and w_vector must be iterables")

    # check no repeated values are passed in the weighting vector
    if len(list(set(w_vector))) != len(w_vector):
        raise ValueError("The weighting vector must be a gradient")

    # no need to weight in case of a one element vector
    if len(w_vector) == 1:
        return vector

    op_vector = []

    # Cast to float so divisions have a floating point resolution
    total_length = float(max(w_vector) - min(w_vector))

    # reflects the expected gradient between subsequent values in w_vector
    # the first value isn't weighted so subtract one from the number of
    # elements
    optimal_gradient = total_length/(len(w_vector)-1)

    for n, vector_value in enumerate(vector):
        # for all elements apply the weighting function
        if n != 0:
            op_vector.append(vector_value * (optimal_gradient) /
                             np.abs((w_vector[n]-w_vector[n-1])))
        # if it's the first element, just return it as is, no weighting to do
        else:
            op_vector.append(vector_value)

    return op_vector


def windowed_diff(vector, window_size):
    """Perform the first difference algorithm between windows of values in a
    vector and each value.

    Parameters
    ----------
    vector: numpy array
        Values to calculate the windowed_diff
    window_size: int or long
        Size of the window

    Returns
    -------
    numpy array
        Array where the Nth value is the difference between the mean of
        vector[N+1:N+1+window_size] and vector[N]. By definition this array
        will have 'window_size' less elements than 'vector'.

    Raises
    ------
    ValueError
        If the window_size is not a positive integer
        If the window_size is greater than the vector length
    """
    # check for consistency in window size and vector size
    if window_size < 1 or not isinstance(window_size, (long, int)):
        raise ValueError("The window_size must be a positive integer")

    if len(vector) <= window_size:
        raise ValueError("The window_size must be smaller than the vector")

    # replicate the last element as many times as required
    for index in range(0, window_size):
        vector = np.append(vector, vector[-1:], axis=0)

    op_vector = []

    for index in range(0, len(vector)-window_size):
        # mean has to be over axis 0 so it handles vectors of vectors
        element = np.mean(vector[(index+1):(index+1+window_size)], axis=0)
        op_vector.append(element - vector[index])

    return op_vector


class VectorResults(namedtuple('VectorResults', ('vector', 'calc',
                                                 'message'))):
    __slots__ = ()  # To avoid creating a dict, as a namedtuple doesn't have it

    def __new__(cls, vector, calc, message):
        return super(VectorResults, cls).__new__(cls, vector, calc, message)


class BaseVectors(object):
    """docstring for BaseVectors"""
    def __init__(self, coord_dict, eigvals, ids, vectors_axes=3,
                 custom_axes=None, weight=False, weighting_vector=None):
        self.message_buffer = []

        # We multiply the coord values with the value of
        # the eigvals represented
        if custom_axes:
            self._vectors = [coord_dict[id_][1:][:vectors_axes] *
                             eigvals[:vectors_axes] for id_ in ids
                             if id_ in coord_dict]
        else:
            self._vectors = [coord_dict[id_][:vectors_axes] *
                             eigvals[:vectors_axes] for id_ in ids
                             if id_ in coord_dict]
        if not self._vectors:
            raise TypeError("No samples to process, an empty list cannot "
                            "be processed")

        # the weighting can only be done over vectors with a length
        # greater than 1
        if weight and weighting_vector is not None and len(ids) > 1:
            try:
                vectors_copy = deepcopy(self._vectors)
                self._vectors = weight_by_vector(vectors_copy,
                                                 weighting_vector)
            except (FloatingPointError, ValueError):
                self.message_buffer.append("Could not weight group, no "
                                           "gradient in the the weighting "
                                           "vector.\n")
                self._vectors = vectors_copy
        elif weight and weighting_vector is None:
            raise ValueError("You must pass a weighting vector if you want to "
                             "weight your data")


class AverageVectors(BaseVectors):
    """docstring for AverageVectors"""

    def results(self):
        """"""
        center = np.average(self._vectors, axis=0)
        if len(self._vectors) == 1:
            vector = [np.linalg.norm(center)]
            calc = {'avg': vector}
        else:
            vector = [np.linalg.norm(i) for i in self._vectors - center]
            calc = {'avg': np.average(vector)}

        msg = ''.join(self.message_buffer) if self.message_buffer else None
        return VectorResults(vector, calc, msg)


class TrajectoryVectors(BaseVectors):
    """"""
    def results(self):
        """"""
        if len(self._vectors) == 1:
            vector = [np.linalg.norm(self._vectors)]
            calc = {'trajectory': vector}
        else:
            vector = [np.linalg.norm(self._vectors[i-1] - self._vectors[i])
                      for i in range(len(self._vectors)-1)]
            calc = {'trajectory': np.linalg.norm(vector)}

        msg = ''.join(self.message_buffer) if self.message_buffer else None
        return VectorResults(vector, calc, msg)


class DifferenceVectors(BaseVectors):
    """"""
    def results(self):
        """"""
        if len(self._vectors) == 1:
            vector = [np.linalg.norm(self._vectors)]
            calc = {'mean': vector, 'std': 0}
        elif len(self._vectors) == 2:
            vector = [np.linalg.norm(self._vectors[1] - self._vectors[0])]
            calc = {'mean': vector, 'std': 0}
        else:
            vec_norm = [np.linalg.norm(self._vectors[i-1] - self._vectors[i])
                        for i in range(len(self._vectors) - 1)]
            vector = np.diff(vec_norm)
            calc = {'mean': np.mean(vector), 'std': np.std(vector)}

        msg = ''.join(self.message_buffer) if self.message_buffer else None
        return VectorResults(vector, calc, msg)


class WindowDifferenceVectors(BaseVectors):
    """docstring for WindowDifferenceVectors"""

    def __init__(self, coord_dict, eigvals, ids, window_size, **kwargs):
        super(WindowDifferenceVectors, self).__init__(coord_dict, eigvals,
                                                      ids, **kwargs)
        self.window_size = window_size

    def results(self):
        """"""
        if len(self._vectors) == 1:
            vector = [np.linalg.norm(self._vectors)]
            calc = {'mean': vector, 'std': 0}
        elif len(self._vectors) == 2:
            vector = [np.linalg.norm(self._vectors[1] - self._vectors[0])]
            calc = {'mean': vector, 'std': 0}
        else:
            vec_norm = [np.linalg.norm(self._vectors[i-1] - self._vectors[i])
                        for i in range(len(self._vectors) - 1)]

            # windowed first differences won't be able on every group,
            # specially given the variation of size that a vector tends to have
            try:
                vector = windowed_diff(vec_norm, self.window_size)
            except ValueError:
                vector = vec_norm
                self.message_buffer.append("Cannot calculate the first "
                                           "difference with a window of size "
                                           "(%d).\n" % self.window_size)
            calc = {'mean': np.mean(vector), 'std': np.std(vector)}

        msg = ''.join(self.message_buffer) if self.message_buffer else None
        return VectorResults(vector, calc, msg)