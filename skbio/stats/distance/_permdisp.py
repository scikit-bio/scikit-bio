# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from scipy.cluster.hierarchy import centroid
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cdist
from scipy.optimize import minimize

import hdmedians as hd 

from ._base import (_preprocess_input, _run_monte_carlo_stats, _build_results)

from skbio.stats.ordination import pcoa
from skbio.util._decorator import experimental

@experimental(as_of="0.5.1")
def permdisp(distance_matrix, grouping, column=None, test='centroid',
                                                    permutations=999):
    """Test for Homogeneity of Multivariate Groups Disperisons using Marti 
    Anderson's PERMDISP2 procedure. PERMDISP is a multivariate analogue of
    Levene's test for homogeneity of multivariate variances. Non-euclidean 
    distances are handled by reducing the original distances to principle 
    coordinates.This procedure has used as a means of assessing beta diversity

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    grouping : 1-D array_like or pandas.DataFrame
        Vector indicating the assignment of objects to groups. For example,
        these could be strings or integers denoting which group an object
        belongs to. If `grouping` is 1-D ``array_like``, it must be the same
        length and in the same order as the objects in `distance_matrix`. If
        `grouping` is a ``DataFrame``, the column specified by `column` will be
        used as the grouping vector. The ``DataFrame`` must be indexed by the
        IDs in `distance_matrix` (i.e., the row labels must be distance matrix
        IDs), but the order of IDs between `distance_matrix` and the
        ``DataFrame`` need not be the same. All IDs in the distance matrix must
        be present in the ``DataFrame``. Extra IDs in the ``DataFrame`` are
        allowed (they are ignored in the calculations).
    column : str, optional
        Column name to use as the grouping vector if `grouping` is a
        ``DataFrame``. Must be provided if `grouping` is a ``DataFrame``.
        Cannot be provided if `grouping` is 1-D ``array_like``.
    test : {'centroid', 'median'}
        determines whether the analysis is done using centroid or spaitial
        median.
    permutations : int, optional
        Number of permutations to use when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.    
    
    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    See Also
    --------
    permanova

    Notes
    -----
    See [1]_ for the original method reference, as well as
    ``vegan::betadisper``, available in R's vegan package [2]_.

    References
    ----------
    .. [1] Anderson, Marti J. "Distance-Based Tests for Homogeneity of 
        Multivariate Dispersions." Biometrics 62 (2006):245-253

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    ...

    """
    #first reduce the matrix to euclidean space using pcoa
    ordination = pcoa(distance_matrix)

    #run preproces_input to obtain variables
    sample_size, num_groups, grouping, tri_idxs, distances = _preprocess_input(
        distance_matrix, grouping, column)
    
    #determine which test type will be used
    if test=='centroid':
        test_stat_function = partial(cen_oneway, ordination)
    elif test=='median':
        test_stat_function = partial(med_oneway, ordination)
    else:
        raise ValueError('Test must be centroid or median')
    
    #compute the stat and pval through monte_carlo tests
    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping,
                                           permutations)

    return _build_results('PERMDISP', 'F-value', sample_size, num_groups,
                          stat, p_value, permutations)

def _compute_centroid_groups(ordination, grouping):

    groups = []
    #group samples in pandas dataframe
    ordination.samples['grouping'] = grouping
    #compute centroids, store in separate dataframe
    centroids = ordination.samples.groupby('grouping').aggregate(_centroid)
    
    grouped = ordination.samples.apply(eu_dist, axis=1,
                                       centroids=centroids).groupby('grouping')
    for _, group in grouped:
        groups.append(group['distance'].tolist())

    return groups

def _centroid(x):
    return x.sum()/len(x)

#returns series of euclidean distances from each point to its centroid
def eu_dist(x, centroids=None):
    return pd.Series([euclidean(x[:-1],
                     centroids.loc[x.grouping]), x.grouping],
                     index=['distance', 'grouping'])

#ship the centroid series to f_oneway for fstat
def cen_oneway(ordination, grouping):
    stat, _ = f_oneway(*(_compute_centroid_groups(ordination, grouping)))
    return stat

def _compute_median_groups(ordination, grouping):
    
    groups = []

    ordination.samples['grouping'] = grouping

    def config_med(x):
        #slice, retype and transpose group array for hd.geomedian
        X = x.values[:, :-1]
        X = X.astype(np.float32)
        return np.array(hd.geomedian(X.T))

    medians = ordination.samples.groupby('grouping').aggregate(config_med)
    
    #return series of euclid distances from each point to its geo-median
    def eu_dist(x):
        return pd.Series([euclidean(x[:-1],
                         medians.loc[x.grouping]), x.grouping],
                         index=['distance', 'grouping'])

    for _, group in ordination.samples.apply(eu_dist,
                                             axis=1).groupby('grouping'):

        groups.append(group['distance'].tolist())
    
    return groups

#ship the median series out to f_oneway
def med_oneway(ordination, grouping):
    stat, _ = f_oneway(*(_compute_median_groups(ordination, grouping)))
    return stat
    

#these functions are for the spatial median computation,
#currently not sure whether to use these or the hdmedians package
#either return the same result
def geometric_mean(points, options={}):
    
    if len(points.shape) == 1:
        # geometric_median((0, 0)) has too much potential for error.
        # Did the user intend a single 2D point or two scalars?
        # Use np.median if you meant the latter.
        raise ValueError("Expected 2D array")
    if points.shape[1] > 2:
        # weiszfeld tends to converge faster in higher dimensions
        method = 'weiszfeld'
    else:
        method = 'minimize'

    return _methods[method](points, options)

def minimize_method(points, options={}):
    """
    Geometric median as a convex optimization problem.
    """

    # objective function
    def aggregate_distance(x):
        return cdist([x], points).sum()

    # initial guess: centroid
    centroid = points.mean(axis=0)

    optimize_result = minimize(aggregate_distance, centroid, method='COBYLA')

    return optimize_result.x

def weiszfeld_method(points, options={}):
    """
    Weiszfeld's algorithm as described on Wikipedia.
    """

    default_options = {'maxiter': 1000, 'tol': 1e-7}
    default_options.update(options)
    options = default_options

    def distance_func(x):
        return cdist([x], points)

    # initial guess: centroid
    guess = points.mean(axis=0)

    iters = 0

    while iters < options['maxiter']:
        distances = distance_func(guess).T

        # catch divide by zero
        # TODO: Wikipedia cites how to deal with distance 0
        distances = np.where(distances == 0, 1, distances)

        guess_next = (points/distances).sum(axis=0) / (1./distances).sum(axis=0)

        guess_movement = np.sqrt(((guess - guess_next)**2).sum())

        guess = guess_next

        if guess_movement <= options['tol']:
            break

        iters += 1

    return guess
    
_methods = {
    'minimize': minimize_method,
    'weiszfeld': weiszfeld_method,
}

