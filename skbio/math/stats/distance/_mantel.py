# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip

from itertools import combinations

import numpy as np
import pandas as pd
import scipy.misc
from scipy.stats import pearsonr, spearmanr

from skbio.core.distance import DistanceMatrix
from .base import p_value_to_str


def mantel(x, y, method='pearson', permutations=999, alternative='two-sided'):
    """Compute correlation between distance matrices using the Mantel test.

    The Mantel test compares two distance matrices by computing the correlation
    between the distances in the lower (or upper) triangular portions of the
    symmetric distance matrices. Correlation can be computed using Pearson's
    product-moment correlation coefficient or Spearman's rank correlation
    coefficient.

    As defined in [1]_, the Mantel test computes a test statistic :math:`r_M`
    given two symmetric distance matrices :math:`D_X` and :math:`D_Y`.
    :math:`r_M` is defined as

    .. math::

       r_M=\\frac{1}{d-1}\\sum_{i=1}^{n-1}\\sum_{j=i+1}^{n}
       stand(D_X)_{ij}stand(D_Y)_{ij}

    where

    .. math::

       d=\\frac{n(n-1)}{2}

    and :math:`n` is the number of rows/columns in each of the distance
    matrices. :math:`stand(D_X)` and :math:`stand(D_Y)` are distance matrices
    with their upper triangles containing standardized distances. Note that
    since :math:`D_X` and :math:`D_Y` are symmetric, the lower triangular
    portions of the matrices could equivalently have been used instead of the
    upper triangular portions (the current function behaves in this manner).

    If ``method='spearman'``, the above equation operates on ranked distances
    instead of the original distances.

    Statistical significance is assessed via a permutation test. The rows and
    columns of the first distance matrix (`x`) are randomly permuted a
    number of times (controlled via `permutations`). A correlation coefficient
    is computed for each permutation and the p-value is the proportion of
    permuted correlation coefficients that are equal to or more extreme
    than the original (unpermuted) correlation coefficient. Whether a permuted
    correlation coefficient is "more extreme" than the original correlation
    coefficient depends on the alternative hypothesis (controlled via
    `alternative`).

    Parameters
    ----------
    x, y : array_like or DistanceMatrix
        Input distance matrices to compare. Both matrices must have the same
        shape and be at least 3x3 in size. If ``array_like``, will be cast to
        ``DistanceMatrix`` (thus the requirements of a valid ``DistanceMatrix``
        apply to both `x` and `y`, such as symmetry and hollowness). If inputs
        are already ``DistanceMatrix`` instances, the IDs do not need to match
        between them; they are assumed to both be in the same order regardless
        of their IDs (the underlying data matrix is the only thing considered
        by this function).
    method : {'pearson', 'spearman'}
        Method used to compute the correlation between distance matrices.
    permutations : int, optional
        Number of times to randomly permute `x` when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.
    alternative : {'two-sided', 'greater', 'less'}
        Alternative hypothesis to use when calculating statistical
        significance. The default ``'two-sided'`` alternative hypothesis
        calculates the proportion of permuted correlation coefficients whose
        magnitude (i.e. after taking the absolute value) is greater than or
        equal to the absolute value of the original correlation coefficient.
        ``'greater'`` calculates the proportion of permuted coefficients that
        are greater than or equal to the original coefficient. ``'less'``
        calculates the proportion of permuted coefficients that are less than
        or equal to the original coefficient.

    Returns
    -------
    tuple of floats
        Correlation coefficient and p-value of the test.

    Raises
    ------
    ValueError
        If `x` and `y` are not the same shape and at least 3x3 in size, or an
        invalid `method`, number of `permutations`, or `alternative` are
        provided.

    See Also
    --------
    DistanceMatrix
    scipy.stats.pearsonr
    scipy.stats.spearmanr

    Notes
    -----
    The Mantel test was first described in [2]_. The general algorithm and
    interface are similar to ``vegan::mantel``, available in R's vegan
    package [3]_.

    ``np.nan`` will be returned for the p-value if `permutations` is zero or if
    the correlation coefficient is ``np.nan``. The correlation coefficient will
    be ``np.nan`` if one or both of the inputs does not have any variation
    (i.e. the distances are all constant) and ``method='spearman'``.

    References
    ----------
    .. [1] Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English
       Edition. Elsevier.

    .. [2] Mantel, N. (1967). "The detection of disease clustering and a
       generalized regression approach". Cancer Research 27 (2): 209-220. PMID
       6018555.

    .. [3] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    Define two 3x3 distance matrices:

    >>> x = [[0, 1, 2],
    ...      [1, 0, 3],
    ...      [2, 3, 0]]
    >>> y = [[0, 2, 7],
    ...      [2, 0, 6],
    ...      [7, 6, 0]]

    Compute the Pearson correlation between them and assess significance using
    a two-sided test with 999 permutations:

    >>> coeff, p_value = mantel(x, y)
    >>> round(coeff, 4)
    0.7559

    Thus, we see a moderate-to-strong positive correlation (:math:`r_M=0.7559`)
    between the two matrices.

    """
    if method == 'pearson':
        corr_func = pearsonr
    elif method == 'spearman':
        corr_func = spearmanr
    else:
        raise ValueError("Invalid correlation method '%s'." % method)

    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if alternative not in ('two-sided', 'greater', 'less'):
        raise ValueError("Invalid alternative hypothesis '%s'." % alternative)

    x = DistanceMatrix(x)
    y = DistanceMatrix(y)

    if x.shape != y.shape:
        raise ValueError("Distance matrices must have the same shape.")
    if x.shape[0] < 3:
        raise ValueError("Distance matrices must be at least 3x3 in size.")

    x_flat = x.condensed_form()
    y_flat = y.condensed_form()

    orig_stat = corr_func(x_flat, y_flat)[0]

    if permutations == 0 or np.isnan(orig_stat):
        p_value = np.nan
    else:
        perm_gen = (corr_func(x.permute(condensed=True), y_flat)[0]
                    for _ in range(permutations))
        permuted_stats = np.fromiter(perm_gen, np.float, count=permutations)

        if alternative == 'two-sided':
            count_better = (np.absolute(permuted_stats) >=
                            np.absolute(orig_stat)).sum()
        elif alternative == 'greater':
            count_better = (permuted_stats >= orig_stat).sum()
        else:
            count_better = (permuted_stats <= orig_stat).sum()

        p_value = (count_better + 1) / (permutations + 1)

    return orig_stat, p_value


def pwmantel(dms, labels=None, method='pearson', permutations=999,
             alternative='two-sided'):
    num_dms = len(dms)

    if num_dms < 2:
        raise ValueError

    if labels is None:
        labels = range(num_dms)
    else:
        if num_dms != len(labels):
            raise ValueError
        if len(set(labels)) != len(labels):
            raise ValueError

    num_combs = scipy.misc.comb(num_dms, 2, exact=True)
    results = np.empty(
        num_combs, dtype=[('dm1', object), ('dm2', object),
                          ('statistic', float), ('p-value', object),
                          ('n', int), ('method', object),
                          ('permutations', int), ('alternative', object)])

    for i, pair in enumerate(combinations(zip(labels, dms), 2)):
        (xlabel, x), (ylabel, y) = pair

        x, y = make_compatible_distance_matrices(x, y)

        if x.ids != y.ids:
            raise ValueError

        stat, p_val = mantel(x, y, method=method, permutations=permutations,
                             alternative=alternative)

        p_val_str = p_value_to_str(p_val, permutations)
        results[i] = (xlabel, ylabel, stat, p_val_str, x.shape[0], method,
                      permutations, alternative)

    return pd.DataFrame.from_records(results, index=('dm1', 'dm2'))


def make_compatible_distance_matrices(dm1, dm2, lookup=None):
    """ Intersect distance matrices and sort the values """
    if lookup:
        try:
            dm1_ids = [lookup[e] for e in dm1.ids]
            dm2_ids = [lookup[e] for e in dm2.ids]
        except KeyError as e:
            raise KeyError("All entries in both DMs must be in lookup if a "
                           "lookup is provided. Missing: %s" % str(e))
        dm1.ids = dm1_ids
        dm2.ids = dm2_ids

    order = [e for e in dm1.ids if e in dm2.ids]

    if len(order) == 0:
        raise ValueError

    # store the intersected distance matrices here
    matrices = []

    # iterate over the distance matrices and identifiers to match the data
    # note that the order must be the same between the two matrices
    for dm in (dm1, dm2):
        # the order is kept by getting the indices from this list
        indices = [dm.ids.index(element) for element in order]

        # this matrix contains the matched up data
        out = dm[indices][:, indices]
        matrices.append(DistanceMatrix(out, order))

    return matrices[0], matrices[1]
