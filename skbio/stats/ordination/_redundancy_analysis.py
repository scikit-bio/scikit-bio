# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from ._base import Ordination, OrdinationResults
from ._utils import corr, svd_rank, scale



def RDA(self, Y, X, site_ids, species_ids, scale_Y=False, scaling=1):
    r"""Compute redundancy analysis, a type of canonical analysis.

    It is related to PCA and multiple regression because the explained
    variables `Y` are fitted to the explanatory variables `X` and PCA
    is then performed on the fitted values. A similar process is
    performed on the residuals.

    RDA should be chosen if the studied gradient is small, and CCA
    when it's large, so that the contingency table is sparse.

    Parameters
    ----------
    Y : array_like
        :math:`n \times p` response matrix. Its columns need be
        dimensionally homogeneous (or you can set `scale_Y=True`).
    X : array_like
        :math:`n \times m, n \geq m` matrix of explanatory
        variables. Its columns need not be standardized, but doing so
        turns regression coefficients into standard regression
        coefficients.
    scale_Y : bool, optional
        Controls whether the response matrix columns are scaled to
        have unit standard deviation. Defaults to `False`.
    site_ids : array_like
        Array of ids for the sites (rows)
    species_ids : array_like
        Array of ids for the speciess (columns)
    scaling : int
        Scaling type 1 produces a distance biplot. It focuses on
        the ordination of rows (sites) because their transformed
        distances approximate their original euclidean
        distances. Especially interesting when most explanatory
        variables are binary.

        Scaling type 2 produces a correlation biplot. It focuses
        on the relationships among explained variables (`Y`). It
        is interpreted like scaling type 1, but taking into
        account that distances between objects don't approximate
        their euclidean distances.

        See more details about distance and correlation biplots in
        [1]_, \S 9.1.4.


    Notes
    -----
    The algorithm is based on [1]_, \S 11.1, and is expected to
    give the same results as ``rda(Y, X)`` in R's package vegan.

    See Also
    --------
    CCA

    References
    ----------
    .. [1] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.

    """

    Y = np.asarray(Y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    return _rda(scale_Y)

def _rda(self, scale_Y, scaling):
    n, p = self.Y.shape
    n_, m = self.X.shape
    if n != n_:
        raise ValueError(
            "Both data matrices must have the same number of rows.")
    if n < m:
        # Mmm actually vegan is able to do this case, too
        raise ValueError(
            "Explanatory variables cannot have less rows than columns.")

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
    B, _, rank_X, _ = np.linalg.lstsq(X, Y)
    Y_hat = X.dot(B)
    # Now let's perform PCA on the fitted values from the multiple
    # regression
    u, s, vt = np.linalg.svd(Y_hat, full_matrices=False)
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
    # site scores. (Eq. 11.12)
    F = Y.dot(U)
    # Ordination in the space of explanatory variables. Its columns
    # are fitted site scores. (Eq. 11.13)
    Z = Y_hat.dot(U)

    # Canonical coefficients (formula 11.14)
    # C = B.dot(U)  # Not used

    Y_res = Y - Y_hat
    # PCA on the residuals
    u_res, s_res, vt_res = np.linalg.svd(Y_res, full_matrices=False)
    # See 9) in p. 587 in L&L 1998
    rank_res = svd_rank(Y_res.shape, s_res)
    # Theoretically, there're at most min(p, n - 1) non-zero eigenvaluesas

    U_res = vt_res[:rank_res].T
    F_res = Y_res.dot(U_res)  # Ordination in the space of residuals

    eigenvalues = np.r_[s[:rank], s_res[:rank_res]]

    return _scores(U, U_res, F, F_res, Z, u, eigenvalues, scaling)

def _scores(self, U, U_res, F, F_res, Z, u, eigenvalues, scaling):
    """Compute site, species and biplot scores for different scalings.

    Parameters
    ----------
    scaling : int

        Scaling type 1 produces a distance biplot. It focuses on
        the ordination of rows (sites) because their transformed
        distances approximate their original euclidean
        distances. Especially interesting when most explanatory
        variables are binary.

        Scaling type 2 produces a correlation biplot. It focuses
        on the relationships among explained variables (`Y`). It
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
        transformed coordinates for species and sites, biplot
        scores, site constraints, etc.

    See Also
    --------
    OrdinationResults

    References
    ----------

    .. [1] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.

    """
    if scaling not in {1, 2}:
        raise NotImplementedError("Only scalings 1, 2 available for RDA.")
    # According to the vegan-FAQ.pdf, the scaling factor for scores
    # is (notice that L&L 1998 says in p. 586 that such scaling
    # doesn't affect the interpretation of a biplot):
    eigvals = eigenvalues
    const = np.sum(eigvals**2)**0.25
    if scaling == 1:
        scaling_factor = const
    elif scaling == 2:
        scaling_factor = eigvals / const
    species_scores = np.hstack((U, U_res)) * scaling_factor
    site_scores = np.hstack((F, F_res)) / scaling_factor
    # TODO not yet used/displayed
    site_constraints = np.hstack((Z, F_res)) / scaling_factor
    # vegan seems to compute them as corr(X[:, :rank_X],
    # u) but I don't think that's a good idea. In fact, if
    # you take the example shown in Figure 11.3 in L&L 1998 you
    # can see that there's an arrow for each of the 4
    # environmental variables (depth, coral, sand, other) even if
    # other = not(coral or sand)
    biplot_scores = corr(X, u)
    # The "Correlations of environmental variables with site
    # scores" from table 11.4 are quite similar to vegan's biplot
    # scores, but they're computed like this:
    # corr(X, F))
    return OrdinationResults(eigvals=eigvals,
                             proportion_explained=eigvals / eigvals.sum(),
                             species=species_scores,
                             site=site_scores,
                             biplot=biplot_scores,
                             site_constraints=site_constraints,
                             site_ids=site_ids,
                             species_ids=species_ids)
