#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio import DistanceMatrix
from numpy import zeros
from scipy.spatial.distance import pdist

def pw_distances(table, ids, metric):
    """Compute distances between all pairs of samples in table

        Parameters
        ----------
        table : 2D np.array or list
        ids : 1D iterable
        metric : str
            The name of the scipy pairwise distance (``pdist``) function
            to use when generating pairwise distances.

        Returns
        -------
        skbio.core.distance.DistanceMatrix

    """
    num_samples = len(ids)
    if num_samples != len(table):
        raise ValueError(
            "Number of columns in table must be equal to number of provided"
            " ids.")

    # initialize the result object
    dm = zeros((num_samples, num_samples))
    for i, id1 in enumerate(ids):
        v1 = table[i]
        for j, id2 in enumerate(ids[:i]):
            v2 = table[j]
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, ids)
