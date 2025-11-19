# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from scipy.linalg import eigh, svd
from scipy.sparse.linalg import svds
from ._ordination_results import OrdinationResults
from skbio.table._tabular import _create_table, _create_table_1d, _ingest_table


def pca(
    X,
    method="svd",
    dimensions=None,
    sample_ids=None,
    feature_ids=None,
    output_format=None,
):
    r"""Perform Principal Component Analysis (PCA).

    Principal component analysis (PCA) is a dimensionality reduction technique
    that finds the directions (principal components) in feature space which
    maximize the variance among samples. The original samples are projected
    onto these directions to obtain a lower-dimensional representation of the
    data that captures as much variance as possible.

    PCA operates on a set of *n* samples that are each associated with a set of
    *p* features. Each sample is then a row of the *n*x*p* data matrix *X*,
    whose columns are the features of the data.
    The goal of PCA is to find directions **w**_i in feature space along which
    the variance of the samples is maximized, the so-called principle components
    of the data matrix:
    max_**w**_i var(X_centered**w**_i)

    The principle components of *X* are the unit eigenvectors of the
    covariance matrix:
    C = 1/(n-1) X_centered^T X_centered    < can use sigma_XX instead
    for which each entry C_ij is the covariance between features i and j
    C_ij = cov(**x**_i, **x**_j)
    The corresponding eigenvalue of each principal direction is the variance
    of the data along that direction.
    sigma_i^2 = var(X_centered**w**_i)

    Parameters
    ----------
    X : table_like
        Samples by features table (n, p). See :ref:`supported formats <table_like>`.
    method : str, optional
        Matrix decomposition method to use. Default is
        ####################################################          for review \/
        Unlike our PCoA, which uses a highly specialized (I think Cython) implementation
        of FSVD, we attempt to use SciPy's sparse SVD "svds" to iteratively compute
        only the specified number of eigenvalues using ARPACK. If dimensions is not
        specified, and the user still wishes to use SVD, we fall back on SciPy's dense
        SVD "svd" which computes all eigenvalues using LAPACK.

        I've heard that SVD is preferred for PCA because it is generally faster and
        more numerically stable than eigendecomposition because it does not compute
        the full covariance matrix, but this needs a source. I've made SVD the
        default method because of this, but I can change it to eigendecomposition
        if that is preferred.
        ####################################################
    dimensions : int, optional
        Number of principal components to compute. Must be a positive integer less
        than or equal to min(n_samples, n_features). If not provided, all principal
        components will be computed.
    sample_ids, feature_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

    Returns
    -------
    OrdinationResults
        Object that stores the PCA results.

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
    Principal Component Analysis (PCA) was first described in [1]_.

    *This section is currently WIP

    References
    ----------

    \/ This is the very first reference I could find for PCA

    .. [1] Pearson, K. (1901). On lines and planes of closest fit to systems of
       points in space. Philosophical Magazine, 2(11), 559–572.


    \/ This is a second reference that generalized PCA to statistical applications,
    although it might be not be relevant or practical to cite such historical sources.
    I can remove them both if that's best.

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

    # Handle method parameter
    if method not in ["eigh", "svd"]:
        raise ValueError("Method must be 'eigh' or 'svd'")

    # Perform eigendecomposition
    elif method == "eigh":
        # Compute the covariance matrix
        covariance_matrix = np.cov(X_centered, rowvar=False)

        # Perform eigendecomposition
        variances, components = eigh(covariance_matrix, subset_by_index=sub_idx)

        # Sort output of eigh in descending order (it is ascending by default)
        components = components.T[::-1]
        variances = variances[::-1]

        # Compute total variance using the fact that the trace
        # of a matrix is the sum of its eigenvalues --
        # We want to avoid this for SVD
        # because it involves computing the covariance matrix
        total_variance = np.trace(covariance_matrix)

    # Perform singular value decomposition
    elif method == "svd":
        # The eigenvalues of the covariance matrix are the right singular values
        # of the centered data matrix scaled by the total degrees of freedom
        # (number of samples - 1)
        if dimensions is None:
            # Perform dense SVD with LAPACK if dimension is not specified
            _, singular_values, components = svd(X_centered)
        else:
            # Perform sparse SVD with ARPACK if dimension is specified
            _, singular_values, components = svds(X_centered, k=dimensions)

            # Sort output of sparse SVD (which may be unordered)
            idx = singular_values.argsort()[::-1]
            singular_values = singular_values[idx]
            components = components[idx, :]

        # Compute explained variances
        variances = np.square(singular_values) / (n - 1)  # limit dimensions

        # The total variance, the trace of the covariance matrix, is equal to
        # the squared Frobenius norm of the centered data matrix scaled by
        # the total degrees of freedom [need source?]
        total_variance = np.sum(np.square(X_centered)) / (n - 1)

    # Project samples onto principal components
    projected_samples = np.dot(X_centered, components.T)

    # Build the OrdinationResults object
    dimensions = n if dimensions is None else dimensions  # [maybe change this]
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
