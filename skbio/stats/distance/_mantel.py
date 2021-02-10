# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import combinations

import warnings
import numpy as np
import pandas as pd
import scipy.special
from scipy.stats import kendalltau
from scipy.stats import PearsonRConstantInputWarning
from scipy.stats import PearsonRNearConstantInputWarning
from scipy.stats import SpearmanRConstantInputWarning

from skbio.stats.distance import DistanceMatrix
from skbio.util._decorator import experimental

from ._cutils import mantel_perm_pearsonr_cy


@experimental(as_of="0.4.0")
def mantel(x, y, method='pearson', permutations=999, alternative='two-sided',
           strict=True, lookup=None):
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
    x, y : DistanceMatrix or array_like
        Input distance matrices to compare. If `x` and `y` are both
        ``DistanceMatrix`` instances, they will be reordered based on matching
        IDs (see `strict` and `lookup` below for handling matching/mismatching
        IDs); thus they are not required to be in the same ID order. If `x` and
        `y` are ``array_like``, no reordering is applied and both matrices must
        have the same shape. In either case, `x` and `y` must be at least 3x3
        in size *after* reordering and matching of IDs.
    method : {'pearson', 'spearman','kendalltau'}
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
    strict : bool, optional
        If ``True``, raises a ``ValueError`` if IDs are found that do not exist
        in both distance matrices. If ``False``, any nonmatching IDs are
        discarded before running the test. See `n` (in Returns section below)
        for the number of matching IDs that were used in the test. This
        parameter is ignored if `x` and `y` are ``array_like``.
    lookup : dict, optional
        Maps each ID in the distance matrices to a new ID. Used to match up IDs
        across distance matrices prior to running the Mantel test. If the IDs
        already match between the distance matrices, this parameter is not
        necessary. This parameter is disallowed if `x` and `y` are
        ``array_like``.

    Returns
    -------
    corr_coeff : float
        Correlation coefficient of the test (depends on `method`).
    p_value : float
        p-value of the test.
    n : int
        Number of rows/columns in each of the distance matrices, after any
        reordering/matching of IDs. If ``strict=False``, nonmatching IDs may
        have been discarded from one or both of the distance matrices prior to
        running the Mantel test, so this value may be important as it indicates
        the *actual* size of the matrices that were compared.

    Raises
    ------
    ValueError
        If `x` and `y` are not at least 3x3 in size after reordering/matching
        of IDs, or an invalid `method`, number of `permutations`, or
        `alternative` are provided.
    TypeError
        If `x` and `y` are not both ``DistanceMatrix`` instances or
        ``array_like``.

    See Also
    --------
    DistanceMatrix
    scipy.stats.pearsonr
    scipy.stats.spearmanr
    pwmantel

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
    Import the functionality we'll use in the following examples:

    >>> from skbio import DistanceMatrix
    >>> from skbio.stats.distance import mantel

    Define two 3x3 distance matrices:

    >>> x = DistanceMatrix([[0, 1, 2],
    ...                     [1, 0, 3],
    ...                     [2, 3, 0]])
    >>> y = DistanceMatrix([[0, 2, 7],
    ...                     [2, 0, 6],
    ...                     [7, 6, 0]])

    Compute the Pearson correlation between them and assess significance using
    a two-sided test with 999 permutations:

    >>> coeff, p_value, n = mantel(x, y)
    >>> print(round(coeff, 4))
    0.7559

    Thus, we see a moderate-to-strong positive correlation (:math:`r_M=0.7559`)
    between the two matrices.

    In the previous example, the distance matrices (``x`` and ``y``) have the
    same IDs, in the same order:

    >>> x.ids
    ('0', '1', '2')
    >>> y.ids
    ('0', '1', '2')

    If necessary, ``mantel`` will reorder the distance matrices prior to
    running the test. The function also supports a ``lookup`` dictionary that
    maps distance matrix IDs to new IDs, providing a way to match IDs between
    distance matrices prior to running the Mantel test.

    For example, let's reassign the distance matrices' IDs so that there are no
    matching IDs between them:

    >>> x.ids = ('a', 'b', 'c')
    >>> y.ids = ('d', 'e', 'f')

    If we rerun ``mantel``, we get the following error notifying us that there
    are nonmatching IDs (this is the default behavior with ``strict=True``):

    >>> mantel(x, y)
    Traceback (most recent call last):
        ...
    ValueError: IDs exist that are not in both distance matrices.

    If we pass ``strict=False`` to ignore/discard nonmatching IDs, we see that
    no matches exist between `x` and `y`, so the Mantel test still cannot be
    run:

    >>> mantel(x, y, strict=False)
    Traceback (most recent call last):
        ...
    ValueError: No matching IDs exist between the distance matrices.

    To work around this, we can define a ``lookup`` dictionary to specify how
    the IDs should be matched between distance matrices:

    >>> lookup = {'a': 'A', 'b': 'B', 'c': 'C',
    ...           'd': 'A', 'e': 'B', 'f': 'C'}

    ``lookup`` maps each ID to ``'A'``, ``'B'``, or ``'C'``. If we rerun
    ``mantel`` with ``lookup``, we get the same results as the original
    example where all distance matrix IDs matched:

    >>> coeff, p_value, n = mantel(x, y, lookup=lookup)
    >>> print(round(coeff, 4))
    0.7559

    ``mantel`` also accepts input that is ``array_like``. For example, if we
    redefine `x` and `y` as nested Python lists instead of ``DistanceMatrix``
    instances, we obtain the same result:

    >>> x = [[0, 1, 2],
    ...      [1, 0, 3],
    ...      [2, 3, 0]]
    >>> y = [[0, 2, 7],
    ...      [2, 0, 6],
    ...      [7, 6, 0]]
    >>> coeff, p_value, n = mantel(x, y)
    >>> print(round(coeff, 4))
    0.7559

    It is import to note that reordering/matching of IDs (and hence the
    ``strict`` and ``lookup`` parameters) do not apply when input is
    ``array_like`` because there is no notion of IDs.

    """
    special = False  # set to true, if we have a dedicated implementation
    if method == 'pearson':
        special = True
    elif method == 'spearman':
        special = True
    elif method == 'kendalltau':
        corr_func = kendalltau
    else:
        raise ValueError("Invalid correlation method '%s'." % method)

    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if alternative not in ('two-sided', 'greater', 'less'):
        raise ValueError("Invalid alternative hypothesis '%s'." % alternative)

    x, y = _order_dms(x, y, strict=strict, lookup=lookup)

    n = x.shape[0]
    if n < 3:
        raise ValueError("Distance matrices must have at least 3 matching IDs "
                         "between them (i.e., minimum 3x3 in size).")

    if special:
        if method == 'pearson':
            orig_stat, permuted_stats = _mantel_stats_pearson(x, y,
                                                              permutations)
        elif method == 'spearman':
            orig_stat, permuted_stats = _mantel_stats_spearman(x, y,
                                                               permutations)
        else:
            raise ValueError("Invalid correlation method '%s'." % method)
    else:
        x_flat = x.condensed_form()
        y_flat = y.condensed_form()

        orig_stat = corr_func(x_flat, y_flat)[0]
        del x_flat

        permuted_stats = []
        if not (permutations == 0 or np.isnan(orig_stat)):
            perm_gen = (corr_func(x.permute(condensed=True), y_flat)[0]
                        for _ in range(permutations))
            permuted_stats = np.fromiter(perm_gen, np.float,
                                         count=permutations)

        del y_flat

    if permutations == 0 or np.isnan(orig_stat):
        p_value = np.nan
    else:
        if alternative == 'two-sided':
            count_better = (np.absolute(permuted_stats) >=
                            np.absolute(orig_stat)).sum()
        elif alternative == 'greater':
            count_better = (permuted_stats >= orig_stat).sum()
        else:
            count_better = (permuted_stats <= orig_stat).sum()

        p_value = (count_better + 1) / (permutations + 1)

    return orig_stat, p_value, n


def _mantel_stats_pearson_flat(x, y_flat, permutations):
    """Compute original and permuted stats using pearsonr.

    Parameters
    ----------
    x : DistanceMatrix
        Input distance matrix.
    y_flat: 1D array
        Compact representation of a distance matrix.
    permutations : int
        Number of times to randomly permute `x` when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and
        permuted_stats will be an empty array.

    Returns
    -------
    orig_stat : 1D array_like
        Correlation coefficient of the test.
    permuted_stats : 1D array_like
        Permuted correlation coefficients of the test.
    """

    x_flat = x.condensed_form()

    # If an input is constant, the correlation coefficient is not defined.
    if (x_flat == x_flat[0]).all() or (y_flat == y_flat[0]).all():
        warnings.warn(PearsonRConstantInputWarning())
        return np.nan, []

    # inline pearsonr, condensed from scipy.stats.pearsonr
    xmean = x_flat.mean()
    xm = x_flat - xmean
    normxm = np.linalg.norm(xm)
    xm_normalized = xm/normxm
    del xm
    del x_flat

    ymean = y_flat.mean()
    ym = y_flat - ymean
    normym = np.linalg.norm(ym)
    ym_normalized = ym/normym
    del ym

    threshold = 1e-13
    if (((normxm < threshold*abs(xmean)) or
         (normym < threshold*abs(ymean)))):
        # If all the values in x (likewise y) are very close to the mean,
        # the loss of precision that occurs in the subtraction xm = x - xmean
        # might result in large errors in r.
        warnings.warn(PearsonRNearConstantInputWarning())

    orig_stat = np.dot(xm_normalized, ym_normalized)

    # Presumably, if abs(orig_stat) > 1, then it is only some small artifact of
    # floating point arithmetic.
    orig_stat = max(min(orig_stat, 1.0), -1.0)

    mat_n = x._data.shape[0]
    # note: xmean and normxm do not change with permutations
    permuted_stats = []
    if not (permutations == 0 or np.isnan(orig_stat)):
        # inline DistanceMatrix.permute, grouping them together
        x_data = x._data
        if not x_data.flags.c_contiguous:
            x_data = np.asarray(x_data, order='C')

        # compute all pearsonr permutations at once
        # create first the list of permutations
        perm_order = np.empty([permutations, mat_n], dtype=np.int)
        for row in range(permutations):
            perm_order[row, :] = np.random.permutation(mat_n)

        permuted_stats = np.empty([permutations], dtype=x_data.dtype)
        mantel_perm_pearsonr_cy(x_data, perm_order, xmean, normxm,
                                ym_normalized, permuted_stats)

    return orig_stat, permuted_stats


def _mantel_stats_pearson(x, y, permutations):
    """Compute original and permuted stats using pearsonr.

    Parameters
    ----------
    x, y : DistanceMatrix
        Input distance matrices to compare.
    permutations : int
        Number of times to randomly permute `x` when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and
        permuted_stats will be an empty array.

    Returns
    -------
    orig_stat : 1D array_like
        Correlation coefficient of the test.
    permuted_stats : 1D array_like
        Permuted correlation coefficients of the test.
    """

    y_flat = y.condensed_form()
    return _mantel_stats_pearson_flat(x, y_flat, permutations)


def _mantel_stats_spearman(x, y, permutations):
    """Compute original and permuted stats using spearmanr.

    Parameters
    ----------
    x, y : DistanceMatrix
        Input distance matrices to compare.
    permutations : int
        Number of times to randomly permute `x` when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and
        permuted_stats will be an empty array.

    Returns
    -------
    orig_stat : 1D array_like
        Correlation coefficient of the test.
    permuted_stats : 1D array_like
        Permuted correlation coefficients of the test.
    """

    x_flat = x.condensed_form()
    y_flat = y.condensed_form()

    # If an input is constant, the correlation coefficient is not defined.
    if (x_flat == x_flat[0]).all() or (y_flat == y_flat[0]).all():
        warnings.warn(SpearmanRConstantInputWarning())
        return np.nan, []

    y_rank = scipy.stats.rankdata(y_flat)
    del y_flat

    x_rank = scipy.stats.rankdata(x_flat)
    del x_flat

    x_rank_matrix = DistanceMatrix(x_rank, x.ids)
    del x_rank

    # for our purposes, spearman is just pearson on rankdata
    return _mantel_stats_pearson_flat(x_rank_matrix, y_rank, permutations)


@experimental(as_of="0.4.0")
def pwmantel(dms, labels=None, method='pearson', permutations=999,
             alternative='two-sided', strict=True, lookup=None):
    """Run Mantel tests for every pair of given distance matrices.

    Runs a Mantel test for each pair of distance matrices and collates the
    results in a ``DataFrame``. Distance matrices do not need to be in the same
    ID order if they are ``DistanceMatrix`` instances. Distance matrices will
    be re-ordered prior to running each pairwise test, and if ``strict=False``,
    IDs that don't match between a pair of distance matrices will be dropped
    prior to running the test (otherwise a ``ValueError`` will be raised if
    there are nonmatching IDs between any pair of distance matrices).

    Parameters
    ----------
    dms : iterable of DistanceMatrix objects, array_like objects, or filepaths
        to distance matrices. If they are ``array_like``, no reordering or
        matching of IDs will be performed.
    labels : iterable of str or int, optional
        Labels for each distance matrix in `dms`. These are used in the results
        ``DataFrame`` to identify the pair of distance matrices used in a
        pairwise Mantel test. If ``None``, defaults to monotonically-increasing
        integers starting at zero.
    method : {'pearson', 'spearman'}
        Correlation method. See ``mantel`` function for more details.
    permutations : int, optional
        Number of permutations. See ``mantel`` function for more details.
    alternative : {'two-sided', 'greater', 'less'}
        Alternative hypothesis. See ``mantel`` function for more details.
    strict : bool, optional
        Handling of nonmatching IDs. See ``mantel`` function for more details.
    lookup : dict, optional
        Map existing IDs to new IDs. See ``mantel`` function for more details.

    Returns
    -------
    pandas.DataFrame
        ``DataFrame`` containing the results of each pairwise test (one per
        row). Includes the number of objects considered in each test as column
        ``n`` (after applying `lookup` and filtering nonmatching IDs if
        ``strict=False``). Column ``p-value`` will display p-values as ``NaN``
        if p-values could not be computed (they are stored as ``np.nan`` within
        the ``DataFrame``; see ``mantel`` for more details).

    See Also
    --------
    mantel
    DistanceMatrix.read

    Notes
    --------
    Passing a list of filepaths can be useful as it allows for a smaller amount
    of memory consumption as it only loads two matrices at a time as opposed to
    loading all distance matrices into memory.

    Examples
    --------
    Import the functionality we'll use in the following examples:

    >>> from skbio import DistanceMatrix
    >>> from skbio.stats.distance import pwmantel

    Define three 3x3 distance matrices:

    >>> x = DistanceMatrix([[0, 1, 2],
    ...                     [1, 0, 3],
    ...                     [2, 3, 0]])
    >>> y = DistanceMatrix([[0, 2, 7],
    ...                     [2, 0, 6],
    ...                     [7, 6, 0]])
    >>> z = DistanceMatrix([[0, 5, 6],
    ...                     [5, 0, 1],
    ...                     [6, 1, 0]])

    Run Mantel tests for each pair of distance matrices (there are 3 possible
    pairs):

    >>> pwmantel((x, y, z), labels=('x', 'y', 'z'),
    ...          permutations=0) # doctest: +NORMALIZE_WHITESPACE
                 statistic p-value  n   method  permutations alternative
    dm1 dm2
    x   y     0.755929     NaN  3  pearson             0   two-sided
        z    -0.755929     NaN  3  pearson             0   two-sided
    y   z    -0.142857     NaN  3  pearson             0   two-sided

    Note that we passed ``permutations=0`` to suppress significance tests; the
    p-values in the output are labelled ``NaN``.

    """
    num_dms = len(dms)

    if num_dms < 2:
        raise ValueError("Must provide at least two distance matrices.")

    if labels is None:
        labels = range(num_dms)
    else:
        if num_dms != len(labels):
            raise ValueError("Number of labels must match the number of "
                             "distance matrices.")
        if len(set(labels)) != len(labels):
            raise ValueError("Labels must be unique.")

    num_combs = scipy.special.comb(num_dms, 2, exact=True)
    results_dtype = [('dm1', object), ('dm2', object), ('statistic', float),
                     ('p-value', float), ('n', int), ('method', object),
                     ('permutations', int), ('alternative', object)]
    results = np.empty(num_combs, dtype=results_dtype)

    for i, pair in enumerate(combinations(zip(labels, dms), 2)):
        (xlabel, x), (ylabel, y) = pair
        if isinstance(x, str):
            x = DistanceMatrix.read(x)
        if isinstance(y, str):
            y = DistanceMatrix.read(y)

        stat, p_val, n = mantel(x, y, method=method, permutations=permutations,
                                alternative=alternative, strict=strict,
                                lookup=lookup)

        results[i] = (xlabel, ylabel, stat, p_val, n, method, permutations,
                      alternative)

    return pd.DataFrame.from_records(results, index=('dm1', 'dm2'))


def _order_dms(x, y, strict=True, lookup=None):
    """Intersect distance matrices and put them in the same order."""
    x_is_dm = isinstance(x, DistanceMatrix)
    y_is_dm = isinstance(y, DistanceMatrix)

    if (x_is_dm and not y_is_dm) or (y_is_dm and not x_is_dm):
        raise TypeError(
            "Mixing DistanceMatrix and array_like input types is not "
            "supported. Both x and y must either be DistanceMatrix instances "
            "or array_like, but not mixed.")
    elif x_is_dm and y_is_dm:
        if lookup is not None:
            x = _remap_ids(x, lookup, 'x', 'first')
            y = _remap_ids(y, lookup, 'y', 'second')

        id_order = [id_ for id_ in x.ids if id_ in y]
        num_matches = len(id_order)

        if (strict and ((num_matches != len(x.ids)) or
                        (num_matches != len(y.ids)))):
            raise ValueError("IDs exist that are not in both distance "
                             "matrices.")

        if num_matches < 1:
            raise ValueError("No matching IDs exist between the distance "
                             "matrices.")

        return x.filter(id_order), y.filter(id_order)
    else:
        # Both x and y aren't DistanceMatrix instances.
        if lookup is not None:
            raise ValueError("ID lookup can only be provided if inputs are "
                             "DistanceMatrix instances.")

        x = DistanceMatrix(x)
        y = DistanceMatrix(y)

        if x.shape != y.shape:
            raise ValueError("Distance matrices must have the same shape.")

        return x, y


def _remap_ids(dm, lookup, label, order):
    "Return a copy of `dm` with its IDs remapped based on `lookup`."""
    try:
        remapped_ids = [lookup[id_] for id_ in dm.ids]
    except KeyError as e:
        raise KeyError("All IDs in the %s distance matrix (%s) must be in "
                       "the lookup. Missing ID: %s" % (order, label, str(e)))

    # Create a copy as we'll be modifying the IDs in place.
    dm_copy = dm.copy()
    dm_copy.ids = remapped_ids
    return dm_copy
