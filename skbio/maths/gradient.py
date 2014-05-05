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
from collections import namedtuple, defaultdict

from skbio.util.sort import signed_natsort


class VectorResults(namedtuple('VectorResults', ('vector', 'calc',
                                                 'message'))):
    __slots__ = ()  # To avoid creating a dict, as a namedtuple doesn't have it

    def __new__(cls, vector, calc, message):
        return super(VectorResults, cls).__new__(cls, vector, calc, message)


class BaseVectors(object):
    r"""Base class for the Vector algorithms

    Parameters
    ----------
    coords : pandas.DataFrame
        The coordinates for each sample id
    eigvals : numpy 1-D array
        The eigenvalues for coords
    metamap : pandas.DataFrame
        The metadata map, indexed by sample ids and columns are metadata
        categories
    vector_categories : list of str, optional
        A list of metadata categories to use to create the vectors. If None is
        passed, the vectors for all metadata categories are computed. Default:
        None, compute all of them
    sort_category : str, optional
        The metadata category to use to sort the vectors. Default: None
    axes : int, optional
        The number of axes to account while doing the vector specific
        calculations. Pass 0 to compute all of them. Default: 3
    weighted : bool, optional
        If true, the output is weighted by the space between samples in the
        `sort_category` column
    """
    def __init__(self, coords, eigvals, metamap, vector_categories=None,
                 sort_category=None, axes=3, weighted=False):
        """
        Raises
        ------
        ValueError
            If vector_categories or sort_category are not in metamap
            If the values under sort_category are not numeric
        """
        if not vector_categories:
            # If vector_categories is not provided, use all the categories
            # present in self._metamap
            vector_categories = metamap.keys()
        else:
            # Check that vector_categories are in metamap
            for category in vector_categories:
                if category not in metamap:
                    raise ValueError("Category %s not present in metadata."
                                     % category)

        # Check that sort_categories is in metamap
        if sort_category and sort_category not in metamap:
            raise ValueError("Sort category %s not present in metadata."
                             % sort_category)

        if axes == 0:
            # If axes == 0, we should compute the vectors for all axes
            axes = len(eigvals)
        elif axes > len(eigvals) or axes < 0:
            # Axes should be 0 <= axes <= len(eigvals)
            raise ValueError("axes should be between 0 and the max number of "
                             "axes available (%d), found: %d "
                             % (len(eigvals), axes))

        # Restrict coordinates to those axes that we actually need to compute
        self._coords = coords.ix[:, :axes-1]
        self._eigvals = eigvals[:axes]
        self._metamap = metamap
        self._weighted = weighted

        # Remove any samples from coords not present in mapping file
        # and remove any samples from metamap not present in coords
        self._normalize_samples()

        # Create groups
        self._make_groups(vector_categories, sort_category)

        # Compute the weighting_vector
        self._weighting_vector = None
        if weighted:
            if not sort_category:
                raise ValueError("You should provide a sort category if you "
                                 "want to weight the vectors")
            try:
                self._weighting_vector = self._metamap[sort_category].astype(
                    np.float64)
            except ValueError:
                    raise ValueError("The sorting category must be numeric")

        # Initialize the message buffer
        self._message_buffer = []

    def _normalize_samples(self):
        """Ensures that self._coords and self._metamap have the same sample_ids

        Raises
        ------
        ValueError
            If `coords` and `metamap` does not have samples in common
        """
        # Figure out the sample ids in common
        coords_sample_ids = set(self._coords.index)
        mm_sample_ids = set(self._metamap.index)
        sample_ids = coords_sample_ids.intersection(mm_sample_ids)

        # Check if they actually have sample ids in common
        if not sample_ids:
            raise ValueError("Coordinates and metadata map had no samples "
                             "in common")

        # Need to take a subset of coords
        if coords_sample_ids != sample_ids:
            self._coords = self._coords.ix[sample_ids]
        # Need to take a subset of metamap
        if mm_sample_ids != sample_ids:
            self._metamap = self._metamap.ix[sample_ids]

    def _make_groups(self, vector_categories, sort_category):
        """Groups the sample ids in metamap by the values in vector_categories

        Creates self._groups a dictionary in which the keys represent the group
        label and values are ordered lists of (sort_value, sample id) tuples

        Parameters
        ----------
        vector_categories : list of str
            A list of metadata categories to use to create the groups. If None
            is passed, the groups for all metadata categories are computed.
            Default: None, compute all of them
        sort_category : str
            The category from self._metamap to use to sort groups
        """
        # If sort_category is provided, we used the value of such category to
        # otherwise we use the sample id
        if sort_category:
            sort_val = lambda sid: self._metamap[sort_category][sid]
        else:
            sort_val = lambda sid: sid

        self._groups = defaultdict(dict)
        for cat in vector_categories:
            # Group samples by category
            gb = self._metamap.groupby(cat)
            for g, df in gb:
                sorted_list = signed_natsort([(sort_val(sid), sid)
                                              for sid in df.index])
                self._groups[cat][g] = [val[1] for val in sorted_list]

    def get_vectors(self):
        """"""
        vector_results = defaultdict(dict)
        # Loop through all the categories that we should compute the vectors
        for cat, cat_groups in self._groups.iteritems():
            # Loop through all the category values present
            # in the current category
            for cat_value, sample_ids in cat_groups.iteritems():
                # Compute the vector for the current category value
                vector_results[cat][cat_value] = self._get_subgroup_vectors(
                    sample_ids)
        return vector_results

    def _get_subgroup_vectors(self, sids):
        """"""
        # We multiply the coord values with the value of
        # the eigvals represented
        vectors = self._coords.ix[sids] * self._eigvals

        if vectors.empty:
            # Raising a RuntimeError since in a usual execution this should
            # never happens. The only way this can happen is if the user
            # directly calls this function, which shouldn't be done
            # (that's why the function is private)
            raise RuntimeError("No samples to process, an empty list cannot "
                               "be processed")

        if self._weighted and len(sids) > 1:
            # the weighting can only be done over vectors with a length
            # greater than 1
            vectors_copy = deepcopy(vectors)
            try:
                vectors = self._weight_by_vector(vectors_copy,
                                                 self._weighting_vector[sids])
            except (FloatingPointError, ValueError):
                self._message_buffer.append("Could not weight group, no "
                                            "gradient in the the "
                                            "weighting vector.\n")
                vectors = vectors_copy

        return self._compute_vector_results(vectors)

    def _compute_vector_results(self, vectors):
        raise NotImplementedError("No algorithm is implemented on the base "
                                  "class.")

    def _weight_by_vector(self, vector, w_vector):
        """weights the values of 'vector' given a weighting vector 'w_vector'.

        Each value in 'vector' will be weighted by the 'rate of change'
        to 'optimal rate of change' ratio, meaning that when calling this
        function over evenly spaced 'w_vector' values, no change will be
        reflected on the output.

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
        if len(set(w_vector)) != len(w_vector):
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

        # for all elements apply the weighting function
        for n, vector_value in enumerate(vector):
            if n == 0:
                # if it's the first element, just return it as is,
                # no weighting to do
                op_vector.append(vector_value)
            else:
                op_vector.append(vector_value * (optimal_gradient) /
                                 np.abs((w_vector[n]-w_vector[n-1])))
        return op_vector


class AverageVectors(BaseVectors):
    """docstring for AverageVectors"""

    def _compute_vector_results(self, vectors):
        """"""
        center = np.average(vectors, axis=0)
        if len(vectors) == 1:
            vector = [np.linalg.norm(center)]
            calc = {'avg': vector}
        else:
            vector = [np.linalg.norm(i) for i in vectors - center]
            calc = {'avg': np.average(vector)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return VectorResults(vector, calc, msg)


class TrajectoryVectors(BaseVectors):
    """"""
    def _compute_vector_results(self, vectors):
        """"""
        if len(vectors) == 1:
            vector = [np.linalg.norm(vectors)]
            calc = {'trajectory': vector}
        else:
            vector = [np.linalg.norm(vectors[i-1] - vectors[i])
                      for i in range(len(vectors)-1)]
            calc = {'trajectory': np.linalg.norm(vector)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return VectorResults(vector, calc, msg)


class DifferenceVectors(BaseVectors):
    """"""
    def _compute_vector_results(self, vectors):
        """"""
        if len(vectors) == 1:
            vector = [np.linalg.norm(vectors)]
            calc = {'mean': vector, 'std': 0}
        elif len(vectors) == 2:
            vector = [np.linalg.norm(vectors[1] - vectors[0])]
            calc = {'mean': vector, 'std': 0}
        else:
            vec_norm = [np.linalg.norm(vectors[i-1] - vectors[i])
                        for i in range(len(vectors) - 1)]
            vector = np.diff(vec_norm)
            calc = {'mean': np.mean(vector), 'std': np.std(vector)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return VectorResults(vector, calc, msg)


class WindowDifferenceVectors(BaseVectors):
    """docstring for WindowDifferenceVectors"""

    def __init__(self, ord_res, metamap, vector_categories, window_size,
                 **kwargs):
        super(WindowDifferenceVectors, self).__init__(ord_res, metamap,
                                                      vector_categories,
                                                      **kwargs)
        self.window_size = window_size

    def _compute_vector_results(self, vectors):
        """"""
        if len(vectors) == 1:
            vector = [np.linalg.norm(vectors)]
            calc = {'mean': vector, 'std': 0}
        elif len(vectors) == 2:
            vector = [np.linalg.norm(vectors[1] - vectors[0])]
            calc = {'mean': vector, 'std': 0}
        else:
            vec_norm = [np.linalg.norm(vectors[i-1] - vectors[i])
                        for i in range(len(vectors) - 1)]

            # windowed first differences won't be able on every group,
            # specially given the variation of size that a vector tends to have
            try:
                vector = self._windowed_diff(vec_norm, self.window_size)
            except ValueError:
                vector = vec_norm
                self._message_buffer.append("Cannot calculate the first "
                                           "difference with a window of size "
                                           "(%d).\n" % self.window_size)
            calc = {'mean': np.mean(vector), 'std': np.std(vector)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return VectorResults(vector, calc, msg)

    def _windowed_diff(self, vector, window_size):
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
