from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn

import numpy as np

from skbio import DistanceMatrix
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform


def pw_distances(counts, ids=None, metric="braycurtis"):
    """Compute distances between all pairs of columns in a counts matrix

    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``.
    metric : str, optional
        The name of the pairwise distance function to use when generating
        pairwise distances. See the scipy ``pdist`` docs, linked under *See
        Also*, for available metrics.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between all pairs of samples (i.e., rows). The number of
        row and columns will be equal to the number of rows in ``counts``.

    Raises
    ------
    ValueError
        If ``len(ids) != len(counts)``.

    See Also
    --------
    scipy.spatial.distance.pdist
    pw_distances_from_table

    """
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError(
            "Number of rows in counts must be equal to number of provided "
            "ids.")

    distances = pdist(counts, metric)
    return DistanceMatrix(
        squareform(distances, force='tomatrix', checks=False), ids)


def pw_distances_from_table(table, metric="braycurtis"):
    """Compute distances between all pairs of samples in table

    Parameters
    ----------
    table : biom.table.Table
        ``Table`` containing count/abundance data of observations across
        samples.
    metric : str, optional
        The name of the pairwise distance function to use when generating
        pairwise distances. See the scipy ``pdist`` docs, linked under *See
        Also*, for available metrics.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between all pairs of samples. The number of row and columns
        will be equal to the number of samples in ``table``.

    See Also
    --------
    scipy.spatial.distance.pdist
    biom.table.Table
    pw_distances

    """
    warn("pw_distances_from_table is deprecated. In the future (tentatively "
         "scikit-bio 0.2.0), pw_distance will take a biom.table.Table object "
         "and this function will be removed. You will need to update your "
         "code to call pw_distances at that time.")
    sample_ids = table.sample_ids
    num_samples = len(sample_ids)

    # initialize the result object
    dm = np.zeros((num_samples, num_samples))
    for i, sid1 in enumerate(sample_ids):
        v1 = table.data(sid1)
        for j, sid2 in enumerate(sample_ids[:i]):
            v2 = table.data(sid2)
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, sample_ids)
