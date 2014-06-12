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


def pwmantel(dms, labels=None, strict=True, lookup=None, method='pearson',
             permutations=999, alternative='two-sided'):
    """Run Mantel tests for every pair of distance matrices.

    Runs a Mantel test for each pair of distance matrices and collates the
    results in a data frame. Distance matrices do not need to be in the same
    ID order (contrary to how the ``mantel`` function behaves). Distance
    matrices will be re-ordered prior to running each pairwise test, and if
    ``strict=False``, IDs that don't match between a pair of distance matrices
    will be dropped prior to running the test (otherwise a ``ValueError`` will
    be raised if there are non-matching IDs between any pair of distance
    matrices).

    Parameters
    ----------
    dms : iterable of DistanceMatrix objects
        DistanceMatrix instances to perform pairwise Mantel tests upon.
    labels : iterable of str or int, optional
        Labels for each ``DistanceMatrix`` in `dms`. These are
        used in the results data frame to identify the pair of distance
        matrices used in a pairwise Mantel test. If ``None``, defaults to
        monotonically-increasing integers starting at zero.
    strict : bool, optional
        If ``True``, raises a ``ValueError`` if IDs are found that do not exist
        in both distance matrices for the current pairwise test. If ``False``,
        these "extra" (nonmatching) IDs are discarded before running the
        pairwise Mantel test.
    lookup : dict, optional
        Maps each ID in the distance matrices to a new ID. Used to match up IDs
        across distance matrices prior to running the Mantel test. If the IDs
        already match between the distance matrices in `dms`, this parameter is
        not necessary.
    method : {'pearson', 'spearman'}
        Correlation method. See ``mantel`` function for more details.
    permutations : int, optional
        Number of permutations. See ``mantel`` function for more details.
    alternative : {'two-sided', 'greater', 'less'}
        Alternative hypothesis. See ``mantel`` function for more details.

    Returns
    -------
    pandas.DataFrame
        Data frame containing the results of each pairwise test (one per row).
        Includes the number of objects considered in each test as column ``n``
        (after applying `lookup` and filtering non-matching IDs if
        ``strict=False``). Column ``p-value`` has the p-values formatted as
        strings with the correct number of decimal places, or ``N/A`` if a
        p-value could not be computed.

    See Also
    --------
    mantel

    """
    num_dms = len(dms)

    if num_dms < 2:
        raise ValueError("Must provide at least two distance matrices.")

    for dm in dms:
        if not isinstance(dm, DistanceMatrix):
            raise TypeError("Must provide DistanceMatrix instances as input.")

    if labels is None:
        labels = range(num_dms)
    else:
        if num_dms != len(labels):
            raise ValueError("Number of labels must match the number of "
                             "distance matrices.")
        if len(set(labels)) != len(labels):
            raise ValueError("Labels must be unique.")

    num_combs = scipy.misc.comb(num_dms, 2, exact=True)
    results_dtype = [('dm1', object), ('dm2', object), ('statistic', float),
                     ('p-value', object), ('n', int), ('method', object),
                     ('permutations', int), ('alternative', object)]
    results = np.empty(num_combs, dtype=results_dtype)

    for i, pair in enumerate(combinations(zip(labels, dms), 2)):
        (xlabel, x), (ylabel, y) = pair

        x, y = _order_dms(x, y, strict=strict, lookup=lookup)

        stat, p_val = mantel(x, y, method=method, permutations=permutations,
                             alternative=alternative)

        p_val_str = p_value_to_str(p_val, permutations)
        results[i] = (xlabel, ylabel, stat, p_val_str, x.shape[0], method,
                      permutations, alternative)

    return pd.DataFrame.from_records(results, index=('dm1', 'dm2'))


def _order_dms(x, y, strict=True, lookup=None):
    """Intersect distance matrices and put them in the same order."""
    if lookup is not None:
        # Create copy as we'll be modifying the IDs in place.
        x = x.copy()
        y = y.copy()

        try:
            x_ids = [lookup[id_] for id_ in x.ids]
            y_ids = [lookup[id_] for id_ in y.ids]
        except KeyError as e:
            raise KeyError("All IDs in both distance matrices must be in the "
                           "lookup. Missing ID: %s" % str(e))
        x.ids = x_ids
        y.ids = y_ids

    id_order = [id_ for id_ in x.ids if id_ in y]
    num_matches = len(id_order)

    if strict and ((num_matches != len(x.ids)) or (num_matches != len(y.ids))):
        raise ValueError("IDs exist that are not in both distance matrices.")

    if num_matches < 1:
        raise ValueError("No matching IDs exist between the distance "
                         "matrices.")

    return x.filter(id_order), y.filter(id_order)
