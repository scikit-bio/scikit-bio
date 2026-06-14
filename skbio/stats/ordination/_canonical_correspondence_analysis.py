# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------


from ._ordination_results import OrdinationResults
from ._utils import corr, svd_rank, scale
from skbio.table._tabular import _create_table, _create_table_1d, _ingest_table
from skbio.util._array import ingest_array


def _xp_svd(xp, mat, full_matrices=False):
    """Perform SVD using scipy, converting to/from xp arrays.

    Parameters
    ----------
    xp : module
        Array API compatible module.
    mat : array
        2-D array.
    full_matrices : bool, optional
        Whether to compute full or reduced SVD.

    Returns
    -------
    u : array
    s : array
    vt : array

    """
    return xp.linalg.svd(mat, full_matrices=full_matrices)


def _xp_lstsq(xp, a, b):
    """Perform least squares using scipy, converting to/from xp arrays.

    Parameters
    ----------
    xp : module
        Array API compatible module.
    a : array
        2-D array (m, n).
    b : array
        2-D array (m, k).

    Returns
    -------
    x : array
        Least squares solution.
    rank : int
        Effective rank of `a`.

    """
    res= xp.linalg.lstsq(a, b, rcond=None)
    result = res[0]
    rank = res[2] if len(res) > 2 else None
    return result, rank

def cca(
    y,
    x,
    scaling=1,
    sample_ids=None,
    feature_ids=None,
    constraint_ids=None,
    output_format=None,
):
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
    y : table_like
        Samples by features table (n, m). See :ref:`supported formats <table_like>`.
    x : table_like
        Samples by constraints table (n, q). See above.
    scaling : int, {1, 2}, optional
        Scaling type 1 maintains :math:`\chi^2` distances between rows.
        Scaling type 2 preserves :math:`\chi^2` distances between columns.
        For a more detailed explanation of the interpretation, check Legendre &
        Legendre 1998, section 9.4.3.
    constraint_ids : list of str, optional
        List of identifiers for metadata variables or constraints. If not provided
        implicitly by the input data structure or explicitly by the user, defaults
        to integers starting at zero.
    sample_ids, feature_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

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
    Y, y_sample_ids, feature_ids = _ingest_table(
        y, sample_ids=sample_ids, feature_ids=feature_ids
    )
    X, x_sample_ids, constraint_ids = _ingest_table(
        x, sample_ids=sample_ids, feature_ids=constraint_ids
    )

    xp, Y, X = ingest_array(Y, X)

    Y = xp.astype(Y, xp.float64)
    X = xp.astype(X, xp.float64)

    # Perform parameter sanity checks
    if X.shape[0] != Y.shape[0]:
        raise ValueError(
            "The samples by features table 'y' and the samples by"
            " constraints table 'x' must have the same number of "
            " rows. 'y': {0} 'x': {1}".format(X.shape[0], Y.shape[0])
        )
    if xp.any(Y < 0):
        raise ValueError(
            "The samples by features table 'y' must be nonnegative"
        )
    row_max = xp.max(Y, axis=1)
    if xp.any(row_max <= 0):
        raise ValueError(
            "The samples by features table 'y' cannot contain a row "
            "with only 0's"
        )
    if scaling not in {1, 2}:
        raise NotImplementedError(
            "Scaling {0} not implemented.".format(scaling)
        )

    # Step 1 (similar to Pearson chi-square statistic)
    grand_total = xp.sum(Y)
    Q = Y / grand_total

    # Features and sample weights (marginal totals)
    column_marginals = xp.sum(Q, axis=0)
    row_marginals = xp.sum(Q, axis=1)

    # Formula 9.32 in Lagrange & Lagrange (1998)
    expected = row_marginals[:, None] * column_marginals[None, :]
    #column marginals can be zero to prevent that
    mask = expected > 0
    safe = xp.where(mask, expected, 1)
    Q_bar = xp.where(mask, (Q - expected) / xp.sqrt(safe), 0)

    # Step 2. Standardize columns of X with respect to sample weights
    X = scale(X, weights=row_marginals, ddof=0)


    # Step 3. Weighted multiple regression
    X_weighted = row_marginals[:, None] ** 0.5 * X
    B, rank_lstsq = _xp_lstsq(xp, X_weighted, Q_bar)
    Y_hat = X_weighted @ B
    Y_res = Q_bar - Y_hat

    # Step 4. Eigenvalue decomposition
    u, s, vt = _xp_svd(xp, Y_hat, full_matrices=False)
    rank = svd_rank(Y_hat.shape, s)
    s = s[:rank]
    u = u[:, :rank]
    vt = vt[:rank]
    U = vt.T

    # Step 5. Eq. 9.38
    eps = xp.finfo(s.dtype).eps
    inv_s = xp.where(s > eps, 1.0 / s, 0.0)
    U_hat = (Q_bar @ U) * inv_s

    # Residuals analysis
    u_res, s_res, vt_res = _xp_svd(xp, Y_res, full_matrices=False)
    rank_res = svd_rank(Y_res.shape, s_res)

    if rank_res > 0:
        s_res = s_res[:rank_res]
        vt_res = vt_res[:rank_res]
        U_res = vt_res.T

        eps_res = xp.finfo(s_res.dtype).eps
        inv_s_res = xp.where(s_res > eps_res, 1.0 / s_res, 0.0)
        U_hat_res = (Y_res @ U_res) * inv_s_res

    else:
        U_res = xp.zeros((Y.shape[1], 0), dtype=Y.dtype)
        U_hat_res = xp.zeros((Y.shape[0], 0), dtype=Y.dtype)
        s_res = xp.zeros((0,), dtype=Y.dtype)

    eigenvalues = xp.concatenate([s, s_res]) ** 2

    # Scalings (p. 596 L&L 1998):
    # feature scores, scaling 1
    V = (column_marginals ** -0.5)[:, None] * U

    # sample scores, scaling 2
    V_hat = (row_marginals ** -0.5)[:, None] * U_hat

    # sample scores, scaling 1
    F = V_hat * s

    # feature scores, scaling 2
    F_hat = V * s

    # Sample scores which are linear combinations of constraint variables
    Z_scaling1 = (row_marginals ** -0.5)[:, None] * (Y_hat @ U)
    Z_scaling2 = Z_scaling1 * inv_s

    # Feature residual scores, scaling 1
    V_res = (column_marginals ** -0.5)[:, None] * U_res

    # Sample residual scores, scaling 2
    V_hat_res = (row_marginals ** -0.5)[:, None] * U_hat_res

    # Sample residual scores, scaling 1
    F_res = V_hat_res * s_res

    # Feature residual scores, scaling 2
    F_hat_res = V_res * s_res

    eigvals = eigenvalues
    if scaling == 1:
        features_scores = xp.concatenate([V, V_res], axis=1)
        sample_scores = xp.concatenate([F, F_res], axis=1)
        sample_constraints = xp.concatenate([Z_scaling1, F_res], axis=1)
    elif scaling == 2:
        features_scores = xp.concatenate([F_hat, F_hat_res], axis=1)
        sample_scores = xp.concatenate([V_hat, V_hat_res], axis=1)
        sample_constraints = xp.concatenate([Z_scaling2, V_hat_res], axis=1)

    biplot_scores = corr(X_weighted, u)

    pc_ids = ["CCA%d" % (i + 1) for i in range(eigenvalues.shape[0])]
    eigvals = _create_table_1d(eigvals, index=pc_ids, backend=output_format)
    samples = _create_table(
        sample_scores, columns=pc_ids, index=y_sample_ids,
        backend=output_format,
    )
    features = _create_table(
        features_scores, columns=pc_ids, index=feature_ids,
        backend=output_format,
    )
    biplot_scores = _create_table(
        biplot_scores,
        index=constraint_ids,
        columns=pc_ids[:biplot_scores.shape[1]],
        backend=output_format,
    )
    sample_constraints = _create_table(
        sample_constraints, index=y_sample_ids, columns=pc_ids,
        backend=output_format,
    )

    return OrdinationResults(
        "CCA",
        "Canonical Correspondence Analysis",
        eigvals,
        samples,
        sample_ids=y_sample_ids,
        features=features,
        feature_ids=feature_ids,
        biplot_scores=biplot_scores,
        sample_constraints=sample_constraints,
        constraint_ids=constraint_ids,
        proportion_explained=eigvals / eigvals.sum(),
    )
