# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from functools import partial

import numpy as np
from scipy.spatial.distance import pdist

from skbio.diversity.beta._unifrac import (
    unweighted_unifrac, weighted_unifrac,_unweighted_unifrac_pdist_f,
    _weighted_unifrac_pdist_f)

from skbio.stats.distance import DistanceMatrix
from skbio.util._decorator import experimental, deprecated


@experimental(as_of="0.4.0")
def beta_diversity(metric, counts, ids=None, **kwargs):
    """Compute distances between all pairs of columns in a counts matrix

    Parameters
    ----------
    metric : str, callable
        The pairwise distance function as a string or callable to use when
        generating pairwise distances. See the scipy ``pdist`` docs and the
        scikit-bio functions linked under *See Also* for available metrics.
        Passing metrics as a string is preferable as this often results in an
        optimized version of the metric being used.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``.

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
    unweighted_unifrac
    weighted_unifrac
    scipy.spatial.distance.pdist
    pw_distances_from_table

    Notes
    -----
    The value that you provide for for ``metric`` can be either a string (e.g.,
    "unweighted_unifrac") or a function
    (e.g., ``skbio.diversity.beta.unweighted_unifrac``). The metric should
    generally be passed as a string, as this often uses an optimized version
    of the metric. For example, passing  ``"unweighted_unifrac"`` (a string)
    will be hundreds of times faster than passing the function
    ``skbio.diversity.beta.unweighted_unifrac``. The latter is faster if
    computing only one or a few distances, but in these cases the difference in
    runtime is negligible, so it's safer to just err on the side of passing
    ``metric`` as a string.

    """
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError(
            "Number of rows in counts must be equal to number of provided "
            "ids.")

    if metric == 'unweighted_unifrac':
        metric, counts, _ = _unweighted_unifrac_pdist_f(
            counts, otu_ids=kwargs['otu_ids'], tree=kwargs['tree'])
    elif metric == 'weighted_unifrac':
        try:
            normalized = kwargs['normalized']
        except KeyError:
            normalized=False
        metric, counts, _ = _weighted_unifrac_pdist_f(
            counts, otu_ids=kwargs['otu_ids'], tree=kwargs['tree'],
            normalized=normalized)
    elif callable(metric):
        metric = partial(metric, **kwargs)
    else:
        pass

    distances = pdist(counts, metric)
    return DistanceMatrix(distances, ids)

pw_distances_from_table_deprecation_reason = (
    "In the future, pw_distance will take a biom.table.Table object "
    "and this function will be removed. You will need to update your "
    "code to call beta_diversity at that time.")


@deprecated(as_of="0.4.0", until="0.4.1",
            reason=pw_distances_from_table_deprecation_reason)
def pw_distances_from_table(table, metric='braycurtis'):
    """Compute distances between all pairs of samples in table

    Parameters
    ----------
    table : biom.table.Table
        ``Table`` containing count/abundance data of observations across
        samples.
    metric : str, callable, optional
        The name of the pairwise distance function to use when generating
        pairwise distances. See the scipy ``pdist`` docs and the scikit-bio
        functions linked under *See Also* for available metrics.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between all pairs of samples. The number of row and columns
        will be equal to the number of samples in ``table``.

    See Also
    --------
    unweighted_unifrac
    weighted_unifrac
    scipy.spatial.distance.pdist
    biom.table.Table
    beta_diversity

    """
    sample_ids = table.ids(axis="sample")
    num_samples = len(sample_ids)

    # initialize the result object
    dm = np.zeros((num_samples, num_samples))
    for i, sid1 in enumerate(sample_ids):
        v1 = table.data(sid1)
        for j, sid2 in enumerate(sample_ids[:i]):
            v2 = table.data(sid2)
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, sample_ids)
