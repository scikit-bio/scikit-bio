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
from scipy.spatial.distance import pdist


def pw_distances(matrix, ids, metric="braycurtis"):
    """Compute distances between all pairs of columns in a matrix

    Parameters
    ----------
    matrix : 2D np.array/list of ints or floats
        Matrix containing count/abundance data where each list/array
        contains counts of observations in a given sample.
    ids : np.array or list, optional
        Identifiers for each sample in ``table``.
    metric : str, optional
        The name of the pairwise distance (``pdist``) function to use when
        generating pairwise distances. See the scipy ``pdist`` docs for
        available choices: http://bit.ly/1nAfiWf

    Returns
    -------
    skbio.core.distance.DistanceMatrix

    Raises
    ------
    ValueError
        If ``len(ids) != len(table)``.

    See Also
    --------
    scipy.spatial.distance.pdist
    pw_distances_from_table

    """
    num_samples = len(ids)
    if num_samples != len(matrix):
        raise ValueError(
            "Number of top-level entries in table must be equal to number of"
            " provided ids.")

    # initialize the result object
    dm = np.zeros((num_samples, num_samples))
    for i, id1 in enumerate(ids):
        v1 = matrix[i]
        for j, id2 in enumerate(ids[:i]):
            v2 = matrix[j]
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, ids)


def pw_distances_from_table(table, metric="braycurtis"):
    """Compute distances between all pairs of columns in a biom_format.Table

    Parameters
    ----------
    table : biom_format.Table
        ``Table`` containing count/abundance data of observations across
        samples.
    metric : str, optional
        The name of the pairwise distance (``pdist``) function to use when
        generating pairwise distances. See the scipy ``pdist`` docs for
        available choices: http://bit.ly/1nAfiWf

    Returns
    -------
    skbio.core.distance.DistanceMatrix

    See Also
    --------
    scipy.spatial.distance.pdist
    biom_format.Table
    pw_distances

    """
    sample_ids = table.sample_ids
    num_samples = len(sample_ids)
    dm = zeros((num_samples, num_samples))
    for i, sid1 in enumerate(sample_ids):
        v1 = table.data(sid1)
        for j, sid2 in enumerate(sample_ids[:i]):
            v2 = table.data(sid2)
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, sample_ids)
