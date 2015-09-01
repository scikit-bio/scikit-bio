# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import six
import numpy as np
from skbio.util._decorator import experimental
from skbio.util._misc import check_random_state
from skbio.stats.distance._mantel import _order_dms

from skbio.stats.distance import DistanceMatrix

@experimental(as_of="0.4.0")
def linregress(y, *args, **kwargs):
    """
    This module performs a multiple linear regression on distance
    matrices.  This is done by extracting the upper (or lower)
    triangular portions of the symmetric distance matrices,
    and collapsing them to compute a multiple linear regression
    on the distances.  The signficances are computed by
    permuting the residuals.

    Parameters
    ----------
    y : DistanceMatrix, array_like, or filepath
        Response distance matrix
    x1, x2, ... : DistanceMatrix, array_like, or filepath
        to distance matrices.  These distance matrices are predictors
        for the regression.  They are also known as covariates
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
    B : np.array
        Array of coefficients
    T : np.array
        Array of pseudo t-statistics for each coefficient
    pvals : np.array
        Array of p-values for each coefficient calculated
        using pseudo T-tests
    F : float
        pseudo F-statistic for lack of fit
    model_pval : float
        pvalue from pseudo F-test for lack of fit
    R2: float
        Coefficient of determination squared

    See Also
    --------
    DistanceMatrix

    References
    ----------
    .. [1] Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English
       Edition. Elsevier.
    .. [2] Legendre, P. Lapointe, F., Casgrain P. (1994) Modeling Brain
       Evolution from Behavior: A Permutational Regression Approach
    .. [3] https://cran.r-project.org/web/packages/ecodist/index.html
    """

    # Unpack kwargs
    params = {'permutations': 10000,
              'random_state': 0,
              'strict': True,
              'lookup': {}}
    for key in ('permutations', 'random_state',
                'strict', 'lookup'):
        params[key] = kwargs.get(key, params[key])
    permutations = params['permutations']
    random_state = params['random_state']
    strict = params['strict']
    lookup = params['lookup']

    random_state = check_random_state(random_state)

    # Read all of the distance matrices from files
    # if necessary
    if isinstance(y, six.string_types):
        y = DistanceMatrix.read(y)
    xargs = []
    for x in args:
        if isinstance(x, six.string_types):
            x = DistanceMatrix.read(x)
        xargs.append(x)

    # Conform all of the ids in the distance matrices to the same order
    if strict:
        for i in range(len(xargs)):
            y, xargs[i] = _order_dms(y, xargs[i],
                                     strict=strict,
                                     lookup=lookup)

    # Linearize all predictor distance matrices into
    # a single matrix
    n = len(y.data)
    X = np.vstack([np.ones((1, n*(n-1)/2))] + \
                  [k.data[np.triu_indices(n, 1)] for k in xargs]).T
    Y = np.atleast_2d(y[np.triu_indices(n, 1)]).T
    n, p = X.shape
    J = np.ones((n, n))
    I = np.identity(n)

    # Define regression function
    def regress(Y, X, computeR=False):
        XX1 = np.linalg.pinv(X.T.dot(X))
        B = XX1.dot(X.T.dot(Y))
        H = X.dot(XX1).dot(X.T)
        Yhat = H.dot(Y)

        SSE = Y.T.dot(I - H).dot(Y)
        SSR = Y.T.dot(H - (1./n)*J).dot(Y)
        dfe = n - p
        dfr = p - 1
        MSR = SSR / dfr
        MSE = SSE / dfe

        T = np.ravel(B) / np.sqrt(np.diag(XX1) * MSE)
        F = MSR / MSE
        if computeR:
            SST = Y.T.dot(I - (1./n)*J).dot(Y)
            R2 = SSR / SST
        else:
            R2 = None
        return Yhat, B, T, F, R2

    # Permutation on residuals
    Yhat, B, T, F, R2 = regress(Y, X, computeR=True)
    E = Y - Yhat
    Fs = np.zeros(permutations)
    Ts = np.zeros((permutations, p))
    for i in range(permutations):
        random_state.shuffle(E)
        Ynew = Yhat + E
        Yhat_, B_, T_, F_, _ = regress(Ynew, X, computeR=False)
        Ts[i, :] = T_
        Fs[i] = F_

    pvals = ((abs(T) >= abs(Ts)).sum(axis=0) + 1) / (permutations + 1)
    model_pval = ((F >= Fs).sum() + 1) / (permutations + 1)

    return (np.ravel(B),
            np.ravel(T),
            pvals,
            np.asscalar(F),
            np.asscalar(model_pval),
            np.asscalar(R2))
