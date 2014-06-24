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


def pw_distances(table, ids=None, metric="braycurtis"):
    """Compute distances between all pairs of samples in table

        Parameters
        ----------
        table : 2D np.array/list of ints or floats, or biom_format.Table
            Table containing count/abundance data where each list/array
            contains counts of observations in a given sample.
        ids : np.array or list, optional
            Identifiers for each sample in ``table``.
        metric : str
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
        TypeError
            If ids was not provided, and ``table`` was not provided as a
            ``biom_format.Table`` object.

    """
    if ids is None:
        try:
            ids = table.sample_ids
        except AttributeError:
            raise TypeError(
                "If not passing ids, table must be a type that defines a "
                "sample_ids attribute, such as a biom_format.Table object.")

    num_samples = len(ids)
    if num_samples != len(table):
        raise ValueError(
            "Number of top-level entries in table must be equal to number of"
            " provided ids.")

    # initialize the result object
    dm = np.zeros((num_samples, num_samples))
    for i, id1 in enumerate(ids):
        v1 = table[i]
        for j, id2 in enumerate(ids[:i]):
            v2 = table[j]
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, ids)
