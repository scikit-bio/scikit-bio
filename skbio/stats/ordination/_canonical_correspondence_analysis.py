# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.linalg import svd, lstsq

from ._ordination_results import OrdinationResults
from ._utils import corr, svd_rank, scale


def cca(y, x, scaling=1):
    r"""Compute canonical (also known as constrained) correspondence analysis.

    Canonical (or constrained) correspondence analysis is a
    multivariate ordination technique. It appeared in community
    ecology [1]_ and relates community composition to the variation in
    the environment (or in other factors). It works from data on
    abundances or counts of samples and constraints variables,
    and outputs ordination axes that maximize sample separation among species.

    It is better suited to extract the niches of taxa than linear
    multivariate methods because it assumes unimodal response curves
    (habitat preferences are often unimodal functions of habitat
    variables [2]_).

    As more environmental variables are added, the result gets more
    similar to unconstrained ordination, so only the variables that
    are deemed explanatory should be included in the analysis.

    Parameters
    ----------
    y : DataFrame
        Samples by features table (n, m)
    x : DataFrame
        Samples by constraints table (n, q)
    scaling : int, {1, 2}, optional
        Scaling type 1 maintains :math:`\chi^2` distances between rows.
        Scaling type 2 preserves :math:`\chi^2` distances between columns.
        For a more detailed explanation of the interpretation, check Legendre &
        Legendre 1998, section 9.4.3.

    Returns
    -------
    OrdinationResults
        Object that stores the cca results.

    Raises
    ------
    ValueError
        If `x` and `y` have different number of rows
        If `y` contains negative values
        If `y` contains a row of only 0's.
    NotImplementedError
        If scaling is not 1 or 2.

    See Also
    --------
    ca
    rda
    OrdinationResults

    Notes
    -----
    The algorithm is based on [3]_, \S 11.2, and is expected to give
    the same results as ``cca(y, x)`` in R's package vegan, except
    that this implementation won't drop constraining variables due to
    perfect collinearity: the user needs to choose which ones to
    input.

    Canonical *correspondence* analysis shouldn't be confused with
    canonical *correlation* analysis (CCorA, but sometimes called
    CCA), a different technique to search for multivariate
    relationships between two datasets. Canonical correlation analysis
    is a statistical tool that, given two vectors of random variables,
    finds linear combinations that have maximum correlation with each
    other. In some sense, it assumes linear responses of "species" to
    "environmental variables" and is not well suited to analyze
    ecological data.

    References
    ----------
    .. [1] Cajo J. F. Ter Braak, "Canonical Correspondence Analysis: A
        New Eigenvector Technique for Multivariate Direct Gradient
        Analysis", Ecology 67.5 (1986), pp. 1167-1179.

    .. [2] Cajo J.F. Braak and Piet F.M. Verdonschot, "Canonical
        correspondence analysis and related multivariate methods in
        aquatic ecology", Aquatic Sciences 57.3 (1995), pp. 255-289.

    .. [3] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.

    """
    Y = y.values
    X = x.values

    # Perform parameter sanity checks
    if X.shape[0] != Y.shape[0]:
        raise ValueError(
            "The samples by features table 'y' and the samples by"
            " constraints table 'x' must have the same number of "
            " rows. 'y': {0} 'x': {1}".format(X.shape[0], Y.shape[0])
        )
    if Y.min() < 0:
        raise ValueError("The samples by features table 'y' must be nonnegative")
    row_max = Y.max(axis=1)
    if np.any(row_max <= 0):
        # Or else the lstsq call to compute Y_hat breaks
        raise ValueError(
            "The samples by features table 'y' cannot contain a " "row with only 0's"
        )
    if scaling not in {1, 2}:
        raise NotImplementedError("Scaling {0} not implemented.".format(scaling))

    # Step 1 (similar to Pearson chi-square statistic)
    grand_total = Y.sum()
    Q = Y / grand_total  # Relative frequencies of Y (contingency table)

    # Features and sample weights (marginal totals)
    column_marginals = Q.sum(axis=0)
    row_marginals = Q.sum(axis=1)

    # Formula 9.32 in Lagrange & Lagrange (1998). Notice that it's an
    # scaled version of the contribution of each cell towards Pearson
    # chi-square statistic.
    expected = np.outer(row_marginals, column_marginals)
    Q_bar = (Q - expected) / np.sqrt(expected)

    # Step 2. Standardize columns of X with respect to sample weights,
    # using the maximum likelihood variance estimator (Legendre &
    # Legendre 1998, p. 595)
    X = scale(X, weights=row_marginals, ddof=0)

    # Step 3. Weighted multiple regression.
    X_weighted = row_marginals[:, None] ** 0.5 * X
    B, _, rank_lstsq, _ = lstsq(X_weighted, Q_bar)
    Y_hat = X_weighted.dot(B)
    Y_res = Q_bar - Y_hat

    # Step 4. Eigenvalue decomposition
    u, s, vt = svd(Y_hat, full_matrices=False)
    rank = svd_rank(Y_hat.shape, s)
    s = s[:rank]
    u = u[:, :rank]
    vt = vt[:rank]
    U = vt.T

    # Step 5. Eq. 9.38
    U_hat = Q_bar.dot(U) * s**-1

    # Residuals analysis
    u_res, s_res, vt_res = svd(Y_res, full_matrices=False)
    rank = svd_rank(Y_res.shape, s_res)
    s_res = s_res[:rank]
    u_res = u_res[:, :rank]
    vt_res = vt_res[:rank]

    U_res = vt_res.T
    U_hat_res = Y_res.dot(U_res) * s_res**-1

    eigenvalues = np.r_[s, s_res] ** 2

    # Scalings (p. 596 L&L 1998):
    # feature scores, scaling 1
    V = (column_marginals**-0.5)[:, None] * U

    # sample scores, scaling 2
    V_hat = (row_marginals**-0.5)[:, None] * U_hat

    # sample scores, scaling 1
    F = V_hat * s

    # feature scores, scaling 2
    F_hat = V * s

    # Sample scores which are linear combinations of constraint
    # variables
    Z_scaling1 = (row_marginals**-0.5)[:, None] * Y_hat.dot(U)
    Z_scaling2 = Z_scaling1 * s**-1

    # Feature residual scores, scaling 1
    V_res = (column_marginals**-0.5)[:, None] * U_res

    # Sample residual scores, scaling 2
    V_hat_res = (row_marginals**-0.5)[:, None] * U_hat_res

    # Sample residual scores, scaling 1
    F_res = V_hat_res * s_res

    # Feature residual scores, scaling 2
    F_hat_res = V_res * s_res

    eigvals = eigenvalues
    if scaling == 1:
        features_scores = np.hstack((V, V_res))
        sample_scores = np.hstack((F, F_res))
        sample_constraints = np.hstack((Z_scaling1, F_res))
    elif scaling == 2:
        features_scores = np.hstack((F_hat, F_hat_res))
        sample_scores = np.hstack((V_hat, V_hat_res))
        sample_constraints = np.hstack((Z_scaling2, V_hat_res))

    biplot_scores = corr(X_weighted, u)

    pc_ids = ["CCA%d" % (i + 1) for i in range(len(eigenvalues))]
    sample_ids = y.index
    feature_ids = y.columns
    eigvals = pd.Series(eigenvalues, index=pc_ids)
    samples = pd.DataFrame(sample_scores, columns=pc_ids, index=sample_ids)
    features = pd.DataFrame(features_scores, columns=pc_ids, index=feature_ids)

    biplot_scores = pd.DataFrame(
        biplot_scores, index=x.columns, columns=pc_ids[: biplot_scores.shape[1]]
    )
    sample_constraints = pd.DataFrame(
        sample_constraints, index=sample_ids, columns=pc_ids
    )

    return OrdinationResults(
        "CCA",
        "Canonical Correspondence Analysis",
        eigvals,
        samples,
        features=features,
        biplot_scores=biplot_scores,
        sample_constraints=sample_constraints,
        proportion_explained=eigvals / eigvals.sum(),
    )
