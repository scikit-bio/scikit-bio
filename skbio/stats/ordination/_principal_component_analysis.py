# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from scipy.linalg import eigh, svd
from scipy.sparse.linalg import eigsh, svds

from ._ordination_results import OrdinationResults
from skbio.table._tabular import _create_table, _create_table_1d, _ingest_table


def pca(
    X,
    method="eigh",
    iterative=False,
    dimensions=None,
    sample_ids=None,
    feature_ids=None,
    output_format=None,
):
    r"""Perform Principal Component Analysis (PCA).

    Principal component analysis (PCA) is a dimensionality reduction technique
    that finds linear combinations of features which maximize the variance among
    samples. The vectors of feature weights resulting from this are the principal
    components of the data. The original samples are projected onto the principal
    components to obtain a lower-dimensional representation of the data that
    captures as much variance as possible.

    PCA operates on a set of `n` samples that are each associated with a set of
    `p` features. Each sample is then a row of the `n` :math:`\times` `p`
    data matrix :math:`\mathbf{X}`, whose columns are the features of the data.
    The goal of PCA is to find vectors :math:`\mathbf{w}_i` in feature space
    along which the variance of the samples is maximized, the principal
    components of the data matrix:

    .. math::
        \mathbf{w}_i^\ast = \arg\max_{\mathbf{w}_i}
        \operatorname{Var}(\mathbf{X_c} \mathbf{w}_i)

    where :math:`\mathbf{X_c}` is the mean-centered data matrix,
    centered by columns (features). The principal components
    of :math:`\mathbf{X}` are the unit eigenvectors of the
    covariance matrix :math:`\mathbf{\Sigma}`:

    .. math::
        \mathbf{\Sigma} = \frac{1}{n-1} \mathbf{X_c}^T
        \mathbf{X_c}

    Each entry :math:`\mathbf{\Sigma}_{ij}` is the covariance between
    features `i` and `j`:

    .. math::
        \mathbf{\Sigma}_{ij} = \operatorname{Cov}(\mathbf{X}_{\cdot i},
        \mathbf{X}_{\cdot j})

    The eigenvalue associated with each principal direction is the variance
    of the data along that direction:

    .. math::
        \sigma_i^2 = \operatorname{Var}(\mathbf{X_c} \mathbf{w}_i)

    Parameters
    ----------
    X : table_like
        Samples by features table (n, p). See :ref:`supported formats <table_like>`.
    method : str, optional
        Matrix decomposition method to use. Default is "eigh" (eigendecomposition),
        which computes exact eigenvectors and eigenvalues of the covariance matrix.
        The alternative is "svd" (singular value decomposition), which bypasses
        computing the full covariance matrix and instead computes the singular
        values of the input matrix.
    iterative : bool, optional
        Whether to use iterative algorithms via ARPACK to compute only the specified
        number of eigenvalues. Default is False, which uses dense algorithms via
        LAPACK to compute all eigenvalues. Only applied if ``dimensions`` is
        specified; otherwise, dense algorithms are used regardless.
    dimensions : int, optional
        Number of principal components to compute. Must be a positive integer less
        than or equal to min(n, p). If not provided, all principal components will
        be computed.
    sample_ids, feature_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

    Returns
    -------
    OrdinationResults
        Object that stores the PCA results, including eigenvalues, the
        proportion of variance explained by each of them, and transformed
        sample coordinates.

    Raises
    ------
    ValueError
        If ``dimensions`` is not a positive integer less than or equal to
        min(n_samples, n_features)
    ValueError
        If ``method`` is not one of "eigh" or "svd"

    See Also
    --------
    OrdinationResults

    Notes
    -----
    Principal Component Analysis (PCA) was first described in [1]_ and [2]_.

    Two methods are provided for performing the matrix decomposition: The
    default method, ``eigh``, constructs the covariance matrix first and then
    performs eigenvalue decomposition to obtain the principal components. The
    alternative method, ``svd``, performs singular value decomposition (SVD)
    on the centered data matrix to obtain the same result. Because SVD avoids
    explicitly computing the covariance matrix, it may be more accurate than
    eigenvalue decomposition for very small eigenvalues.

    The ``iterative`` parameter determines whether the methods are performed
    by dense algorithms, which compute all eigenvalues to full accuracy via
    LAPACK, or sparse algorithms, which compute only the specified number of
    eigenvalues iteratively by ARPACK. Dense algorithms are used by default,
    while iterative algorithms are preferred for very large data matrices
    when the desired number of principal components is much less than the
    size of the data matrix.

    The number of dimensions is determined by ``dimensions``, and if no
    dimension is specified, the default is to compute all eigenvalues.
    The iterative parameter is incompatible with computing all eigenvalues;
    therefore, if no dimension is specified, dense algorithms are used
    regardless.

    If ``iterative`` is specified as True, ``svd`` is computed by SciPy's
    dense "svd" function, and ``eigh`` is computed by SciPy's dense "eigh"
    function. If ``iterative`` is specified as False, the methods are computed
    by SciPy's sparse "svds" and "eigsh" functions respectively.

    References
    ----------
    .. [1] Pearson, K. (1901). On lines and planes of closest fit to systems of
       points in space. Philosophical Magazine, 2(11), 559–572.

    .. [2] Hotelling, H. (1933). Analysis of a complex of statistical variables
       into principal components. Journal of Educational Psychology, 24(6),
       417–441, 498–520.

    """

    # Ingestion of input data matrix
    X, row_ids, column_ids = _ingest_table(X, sample_ids, feature_ids)

    # Center the data matrix
    X_centered = X - np.mean(X, axis=0)
    n = X.shape[0]

    # Handle dimensions parameter
    if dimensions is None:
        sub_idx = None
    elif dimensions <= 0 or dimensions > min(X.shape):
        raise ValueError(
            "Dimensions must be a positive integer less than or equal to "
            "min(n_samples, n_features)"
        )
    elif type(dimensions) is not int:
        raise ValueError("Dimensions must be a positive integer")
    else:
        sub_idx = [n - dimensions, n - 1]

    # Perform eigendecomposition
    if method == "eigh":
        # Compute the covariance matrix
        covariance_matrix = np.cov(X_centered, rowvar=False)

        # Perform iterative sparse eigendecomposition only if "iterative"
        # is selected and a valid dimension is specified
        if not iterative or dimensions is None:
            # Perform full dense eigendecomposition
            variances, components = eigh(covariance_matrix, subset_by_index=sub_idx)

            # Sort output of eigh in descending order (it is ascending by default)
            components = components.T[::-1]
            variances = variances[::-1]
        else:
            # Perform iterative sparse eigendecomposition
            variances, components = eigsh(covariance_matrix, k=dimensions)
            components = components.T

        # Compute total variance using the fact that the trace
        # of a matrix is the sum of its eigenvalues --
        # We want to avoid this for SVD
        # because it involves computing the covariance matrix
        total_variance = np.trace(covariance_matrix)

    # Perform singular value decomposition
    elif method == "svd":
        # The eigenvalues of the covariance matrix are the right singular values
        # of the centered data matrix, scaled by the total degrees of freedom
        # (number of samples - 1)

        # Perform iterative sparse SVD only if "iterative" is selected and
        # a valid dimension is specified
        if not iterative or dimensions is None:
            # Perform full dense SVD
            _, singular_values, components = svd(X_centered)

            # Collect only the specified number of components if applicable
            if dimensions:
                singular_values = singular_values[:dimensions]
                components = components[:dimensions]
        else:
            # Perform iterative sparse SVD
            _, singular_values, components = svds(X_centered, k=dimensions)

        # Compute explained variances
        variances = np.square(singular_values) / (n - 1)

        # The total variance, the trace of the covariance matrix, is equal to
        # the squared Frobenius norm of the centered data matrix scaled by
        # the total degrees of freedom [need source?]
        total_variance = np.sum(np.square(X_centered)) / (n - 1)

    # Handle method parameter
    else:
        raise ValueError("Method must be 'eigh' or 'svd'")

    # Iterative algorithms may return an unsorted output
    if iterative:
        idx = variances.argsort()[::-1]
        variances = variances[idx]
        components = components[idx, :]

    # Project samples onto principal components
    projected_samples = np.dot(X_centered, components.T)

    # Build the OrdinationResults object
    pc_ids = ["%s%d" % ("PC", i + 1) for i in range(variances.shape[0])]
    eigvals = _create_table_1d(variances, index=pc_ids, backend=output_format)
    samples = _create_table(
        projected_samples, index=row_ids, columns=pc_ids, backend=output_format
    )
    features = _create_table(
        components, index=pc_ids, columns=column_ids, backend=output_format
    )

    return OrdinationResults(
        short_method_name="PCA",
        long_method_name="Principal Component Analysis",
        eigvals=eigvals,
        samples=samples,
        sample_ids=row_ids,
        features=features,
        feature_ids=pc_ids,
        proportion_explained=eigvals / total_variance,
    )
