# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import copy
import numpy as np
import pandas as pd
from skbio.util._decorator import experimental
from skbio.util._misc import check_random_state
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance._mantel import _order_dms, _check_dm_labels
from scipy.spatial.distance import cityblock


@experimental(as_of="0.4.0")
def mrm(y, *args, **kwargs):
    """
    This module performs a multiple linear regression on distance
    matrices.  This is done by extracting the upper (or lower)
    triangular portions of the symmetric distance matrices,
    and collapsing them to compute a multiple linear regression
    on the distances.  The signficances are computed by
    permuting the residuals.

    Parameters
    ----------
    y : DistanceMatrix, array_like
        Response distance matrix
    x1, x2, ... : DistanceMatrix, array_like
        to distance matrices.  These distance matrices are predictors
        for the regression.  They are also known as covariates
    labels : iterable of str or int, optional
        Labels for each distance matrix in `args`. These are used in the
        results ``DataFrame`` to identify the distance matrices in
        `x1, x2, ...`. If ``None``, defaults to monotonically-increasing
        integers starting at zero.
    permutations : int, optional
        Number of times to randomly permute `x` when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.
    random_state: int or np.RandomState, optional
        Pseudo number generator state used for random sampling.
    strict : bool, optional
        If ``True``, raises a ``ValueError`` if IDs are found that do not exist
        in all of the distance matrices. If ``False``, any nonmatching IDs are
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
    coefs : pd.DataFrame
       Contains information for each of the coefficients
       including the coefficient, the t-statistic and the p-value
       the coefficient contributes to the models
    summary : pd.Series
       Contains summary statistics, such as F-statistic,
       lack-of-fit p-value and R^2

    See Also
    --------
    DistanceMatrix

    References
    ----------
    .. [1] Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English
       Edition. Elsevier.
    .. [2] Legendre, P. Lapointe, F., Casgrain P. (1994) Modeling Brain
       Evolution from Behavior: A Permutational Regression Approach
    .. [3] Lichstein, J. (2007) Multiple regression on distance matrices: a
           multivariate spatial analysis tool
    .. [4] https://cran.r-project.org/web/packages/ecodist/index.html

    Examples
    --------
    Import the functionality we'll use in the following examples

    >>> from skbio import DistanceMatrix
    >>> from skbio.stats.distance import mrm

    Define two 3x3 distance matrices

    >>> x = DistanceMatrix([[0, 1, 2.1],
    ...                     [1, 0, 2.9],
    ...                     [2.1, 2.9, 0]])
    >>> y = DistanceMatrix([[0, 2, 4],
    ...                     [2, 0, 6],
    ...                     [4, 6, 0]])

    Now perform a multiple regression fit
    >>> B, T, pvals, F, model_pval, R2 = mrm(y, x)
    >>> print(B)
    intercept   -0.175824
    0            2.087912
    dtype: float64
    """
    # Unpack kwargs
    params = {'permutations': 999,
              'random_state': 0,
              'strict': True,
              'lookup': None,
              'labels': None}
    for key in ('permutations', 'random_state',
                'strict', 'lookup', 'labels'):
        params[key] = kwargs.get(key, params[key])

    permutations = params['permutations']
    random_state = params['random_state']
    strict = params['strict']
    lookup = params['lookup']
    labels = params['labels']

    labels = _check_dm_labels(labels, len(args))
    random_state = check_random_state(random_state)
    xargs = list(args)
    # Conform all of the ids in the distance matrices to the same order
    # Intersect all covariates
    for i in range(len(xargs)):
        for j in range(i):
            xargs[i], xargs[j] = _order_dms(xargs[i],
                                            xargs[j],
                                            strict=strict,
                                            lookup=lookup)
    # Intersect response matrix against all covariates
    for i in range(len(xargs)):
        y, xargs[i] = _order_dms(y, xargs[i],
                                 strict=strict,
                                 lookup=lookup)

    # Linearize all predictor distance matrices into
    # a single matrix
    n = len(y.data)
    X = np.vstack([np.ones((1, n*(n-1)/2))] +
                  [k.data[np.triu_indices(n, 1)] for k in xargs]).T
    Y = np.atleast_2d(y[np.triu_indices(n, 1)]).T
    cY = copy.deepcopy(Y)
    n, p = X.shape
    # Define regression function
    XX1 = np.linalg.pinv(X.T.dot(X))
    dfe, dfr = n - p,  p - 1

    def regress(Y):
        B = XX1.dot(X.T.dot(Y))
        Yhat = X.dot(B)
        _E = Y - Yhat
        mY = Y.mean()
        SSE = (_E**2).sum()
        SST = ((Y - mY)**2).sum()
        SSR = SST - SSE
        MSR, MSE = SSR / dfr, SSE / dfe
        T = np.ravel(B) / np.sqrt(np.diag(XX1) * MSE)
        F = MSR / MSE
        R2 = SSR / SST
        return Yhat, B, _E, T, F, R2

    # Permutation on labels
    Yhat, B, E, T, F, R2 = regress(Y)

    Fs = np.zeros(permutations)
    Ts = np.zeros((permutations, p))
    for i in range(permutations):
        random_state.shuffle(cY)
        _, _, _, T_, F_,  _ = regress(cY)
        Ts[i, :], Fs[i] = T_, F_
    # Calculate result statistics
    pvals = ((abs(T) <= abs(Ts)).sum(axis=0) + 1) / (permutations + 1)
    model_pval = ((F <= Fs).sum() + 1) / (permutations + 1)
    labs = ['intercept'] + list(labels)
    data = np.array([np.ravel(B),
                     np.ravel(T),
                     np.ravel(pvals)]).T
    coefs = pd.DataFrame(data,
                         index=labs,
                         columns=['coefficient',
                                  't-statistic',
                                  'p-value'])
    summary = pd.Series([F, model_pval, R2],
                        index=['F-statistic',
                               'lack-of-fit p-value',
                               'R^2'])
    return (coefs, summary)


@experimental(as_of="0.4.0")
def make_categorical_dms(x, metric=cityblock, ignore_nans=True):
    """
    Creates multiple distance matrices from a categorical vector.
    If vector has 3 categories A, B, C, then 3 different
    categorical indicator vectors will be created. From these categorical
    vectors, distance matrices will be created for the distances between
    each pair of categorical indicator vectors.

    Parameters
    ----------
    x : pd.Series or array_like
       1-D categorical vector.
    metric: function, option
       The distance metric to use to compare two categorical vectors
       default: scipy.spatial.distance.cityblock

    Returns
    -------
    iterable
       Iterable of tuples that contains
       1. skbio.DistanceMatrix
       2. string identifier
    """
    x = pd.Series(x)
    if ignore_nans:
        y = x[~pd.isnull(x)]

    cats = np.unique(y)
    for i in range(len(cats)):
        for j in range(i):
            a, b = (x == cats[j]), (x == cats[i])
            dm = DistanceMatrix.from_iterable(np.vstack([a, b]).T,
                                              metric=metric,
                                              keys=x.index)
            yield (dm, '%s:%s' % (str(cats[i]), str(cats[j])))
