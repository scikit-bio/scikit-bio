# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from functools import partial

from scipy.spatial.distance import pdist

from skbio.diversity.beta._unifrac import (_unweighted_unifrac_pdist_f,
                                           _weighted_unifrac_pdist_f)

from skbio.stats.distance import DistanceMatrix
from skbio.util._decorator import experimental


@experimental(as_of="0.4.0")
def beta_diversity(metric, counts, ids=None, **kwargs):
    """Compute distances between all pairs of columns in a counts matrix

    Parameters
    ----------
    metric : str, callable
        The pairwise distance function to apply. See the scipy ``pdist`` docs
        and the scikit-bio functions linked under *See Also* for available
        metrics. Passing metrics as a string is preferable as this often
        results in an optimized version of the metric being used.
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
    skbio.diversity.alpha.alpha_diversity

    Notes
    -----
    The value that you provide for for ``metric`` can be either a string (e.g.,
    ``"unweighted_unifrac"``) or a function
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
        normalized = kwargs.get('normalized', False)
        metric, counts, _ = _weighted_unifrac_pdist_f(
            counts, otu_ids=kwargs['otu_ids'], tree=kwargs['tree'],
            normalized=normalized)
    elif callable(metric):
        metric = partial(metric, **kwargs)
    else:
        pass

    distances = pdist(counts, metric)
    return DistanceMatrix(distances, ids)
