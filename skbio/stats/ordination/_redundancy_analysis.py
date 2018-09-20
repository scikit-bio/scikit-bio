# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.linalg import svd, lstsq

from skbio.util._decorator import experimental
from ._ordination_results import OrdinationResults
from ._utils import corr, svd_rank, scale


@experimental(as_of="0.4.0")
def rda(y, x, scale_Y=False, scaling=1):
    r"""Compute redundancy analysis, a type of canonical analysis.

    It is related to PCA and multiple regression because the explained
    variables `y` are fitted to the explanatory variables `x` and PCA
    is then performed on the fitted values. A similar process is
    performed on the residuals.

    RDA should be chosen if the studied gradient is small, and CCA
    when it's large, so that the contingency table is sparse.

    Parameters
    ----------
    y : pd.DataFrame
        :math:`n \times p` response matrix, where :math:`n` is the number
        of samples and :math:`p` is the number of features. Its columns
        need be dimensionally homogeneous (or you can set `scale_Y=True`).
        This matrix is also referred to as the community matrix that
        commonly stores information about species abundances
    x : pd.DataFrame
        :math:`n \times m, n \geq m` matrix of explanatory
        variables, where :math:`n` is the number of samples and
        :math:`m` is the number of metadata variables. Its columns
        need not be standardized, but doing so turns regression
        coefficients into standard regression coefficients.
    scale_Y : bool, optional
        Controls whether the response matrix columns are scaled to
        have unit standard deviation. Defaults to `False`.
    scaling : int
        Scaling type 1 produces a distance biplot. It focuses on
        the ordination of rows (samples) because their transformed
        distances approximate their original euclidean
        distances. Especially interesting when most explanatory
        variables are binary.

        Scaling type 2 produces a correlation biplot. It focuses
        on the relationships among explained variables (`y`). It
        is interpreted like scaling type 1, but taking into
        account that distances between objects don't approximate
        their euclidean distances.

        See more details about distance and correlation biplots in
        [1]_, \S 9.1.4.

    Returns
    -------
    OrdinationResults
        Object that stores the computed eigenvalues, the
        proportion explained by each of them (per unit),
        transformed coordinates for feature and samples, biplot
        scores, sample constraints, etc.

    See Also
    --------
    ca
    cca
    OrdinationResults

    Notes
    -----
    The algorithm is based on [1]_, \S 11.1, and is expected to
    give the same results as ``rda(y, x)`` in R's package vegan.
    The eigenvalues reported in vegan are re-normalized to
    :math:`\sqrt{\frac{s}{n-1}}` `n` is the number of samples,
    and `s` is the original eigenvalues. Here we will only return
    the original eigenvalues, as recommended in [1]_.

    References
    ----------
    .. [1] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.

    """
    Y = y.values
    X = x.values

    n, p = y.shape
    n_, m = x.shape
    if n != n_:
        raise ValueError(
            "Both data matrices must have the same number of rows.")
    if n < m:
        # Mmm actually vegan is able to do this case, too
        raise ValueError(
            "Explanatory variables cannot have less rows than columns.")

    sample_ids = y.index
    feature_ids = y.columns
    # Centre response variables (they must be dimensionally
    # homogeneous)
    Y = scale(Y, with_std=scale_Y)
    # Centre explanatory variables
    X = scale(X, with_std=False)

    # Distribution of variables should be examined and transformed
    # if necessary (see paragraph 4 in p. 580 L&L 1998)

    # Compute Y_hat (fitted values by multivariate linear
    # regression, that is, linear least squares). Formula 11.6 in
    # L&L 1998 involves solving the normal equations, but that fails
    # when cond(X) ~ eps**(-0.5). A more expensive but much more
    # stable solution (fails when cond(X) ~ eps**-1) is computed
    # using the QR decomposition of X = QR:
    # (11.6) Y_hat = X [X' X]^{-1} X' Y
    #              = QR [R'Q' QR]^{-1} R'Q' Y
    #              = QR [R' R]^{-1} R'Q' Y
    #              = QR R^{-1} R'^{-1} R' Q' Y
    #              = Q Q' Y
    # and B (matrix of regression coefficients)
    # (11.4) B = [X' X]^{-1} X' Y
    #          = R^{-1} R'^{-1} R' Q' Y
    #          = R^{-1} Q'
    # Q, R = np.linalg.qr(X)
    # Y_hat = Q.dot(Q.T).dot(Y)
    # B = scipy.linalg.solve_triangular(R, Q.T.dot(Y))
    # This works provided X has full rank. When not, you can still
    # fix it using R's pseudoinverse or partitioning R. To avoid any
    # issues, like the numerical instability when trying to
    # reproduce an example in L&L where X was rank-deficient, we'll
    # just use `np.linalg.lstsq`, which uses the SVD decomposition
    # under the hood and so it's also more expensive.
    B, _, rank_X, _ = lstsq(X, Y)
    Y_hat = X.dot(B)
    # Now let's perform PCA on the fitted values from the multiple
    # regression
    u, s, vt = svd(Y_hat, full_matrices=False)
    # vt are the right eigenvectors, which is what we need to
    # perform PCA. That is, we're changing points in Y_hat from the
    # canonical basis to the orthonormal basis given by the right
    # eigenvectors of Y_hat (or equivalently, the eigenvectors of
    # the covariance matrix Y_hat.T.dot(Y_hat))
    # See 3) in p. 583 in L&L 1998
    rank = svd_rank(Y_hat.shape, s)
    # Theoretically, there're at most min(p, m, n - 1) non-zero eigenvalues

    U = vt[:rank].T  # U as in Fig. 11.2

    # Ordination in the space of response variables. Its columns are
    # sample scores. (Eq. 11.12)
    F = Y.dot(U)
    # Ordination in the space of explanatory variables. Its columns
    # are fitted sample scores. (Eq. 11.13)
    Z = Y_hat.dot(U)

    # Canonical coefficients (formula 11.14)
    # C = B.dot(U)  # Not used

    Y_res = Y - Y_hat
    # PCA on the residuals
    u_res, s_res, vt_res = svd(Y_res, full_matrices=False)
    # See 9) in p. 587 in L&L 1998
    rank_res = svd_rank(Y_res.shape, s_res)
    # Theoretically, there're at most min(p, n - 1) non-zero eigenvalues as

    U_res = vt_res[:rank_res].T
    F_res = Y_res.dot(U_res)  # Ordination in the space of residuals

    eigenvalues = np.r_[s[:rank], s_res[:rank_res]]

    # Compute scores
    if scaling not in {1, 2}:
        raise NotImplementedError("Only scalings 1, 2 available for RDA.")
    # According to the vegan-FAQ.pdf, the scaling factor for scores
    # is (notice that L&L 1998 says in p. 586 that such scaling
    # doesn't affect the interpretation of a biplot):
    eigvals = pd.Series(
        eigenvalues, index=['RDA%d' % (i+1) for i in range(len(eigenvalues))])
    const = np.sum(eigenvalues**2)**0.25
    if scaling == 1:
        scaling_factor = const
    elif scaling == 2:
        scaling_factor = eigenvalues / const
    feature_scores = np.hstack((U, U_res)) * scaling_factor
    sample_scores = np.hstack((F, F_res)) / scaling_factor

    feature_scores = pd.DataFrame(
        feature_scores, index=feature_ids,
        columns=['RDA%d' % (i+1) for i in range(feature_scores.shape[1])])
    sample_scores = pd.DataFrame(
        sample_scores, index=sample_ids,
        columns=['RDA%d' % (i+1) for i in range(sample_scores.shape[1])])
    # TODO not yet used/displayed
    sample_constraints = np.hstack((Z, F_res)) / scaling_factor
    sample_constraints = pd.DataFrame(
        sample_constraints, index=sample_ids,
        columns=['RDA%d' % (i+1) for i in range(sample_constraints.shape[1])])
    # Vegan seems to compute them as corr(X[:, :rank_X],
    # u) but I don't think that's a good idea. In fact, if
    # you take the example shown in Figure 11.3 in L&L 1998 you
    # can see that there's an arrow for each of the 4
    # environmental variables (depth, coral, sand, other) even if
    # other = not(coral or sand)
    biplot_scores = corr(X, u)
    biplot_scores = pd.DataFrame(
        biplot_scores, index=x.columns,
        columns=['RDA%d' % (i+1) for i in range(biplot_scores.shape[1])])
    # The "Correlations of environmental variables with sample
    # scores" from table 11.4 are quite similar to vegan's biplot
    # scores, but they're computed like this:
    # corr(X, F))
    p_explained = pd.Series(
        eigenvalues / eigenvalues.sum(),
        index=['RDA%d' % (i+1) for i in range(len(eigenvalues))])
    return OrdinationResults('RDA', 'Redundancy Analysis',
                             eigvals=eigvals,
                             proportion_explained=p_explained,
                             features=feature_scores,
                             samples=sample_scores,
                             biplot_scores=biplot_scores,
                             sample_constraints=sample_constraints)
