# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Robust Principal Component Analysis (RPCA).

This module implements Robust PCA for sparse compositional data analysis. RPCA utilizes
matrix completion with OptSpace as a pre-processing step to handle missing entries
before dimensionality reduction. This may be used to analyze compositional data which
has already been processed using the robust centered log-ratio (rclr) transformation.

References
----------
.. [1] Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
   Knight, R., & Zengler, K. (2019). A novel sparse compositional technique reveals
   microbial perturbations. mSystems, 4(1), 10-1128.

"""

import numpy as np

from ._ordination_results import OrdinationResults
from ._optspace import optspace
from ._principal_component_analysis import _pca
from skbio.table._tabular import (
    _create_table,
    _create_table_1d,
    _ingest_table,
)


def rpca(
    X,
    dimensions=3,
    max_iter=10000,
    sample_ids=None,
    feature_ids=None,
    output_format=None,
):
    r"""Perform Robust Principal Component Analysis.

    Robust PCA (RPCA) is an ordination method for data with missing entries.
    In the context of compositional data, RPCA is used to analyze data that has
    been transformed via the robust centered log-ratio (rclr) transformation,
    in which zeros become missing entries. In RPCA, OptSpace matrix completion is
    applied to impute missing entries, and PCA is performed on the completed matrix.

    Parameters
    ----------
    X : table_like of shape (n_samples, n_features)
        Input data table. See :ref:`supported formats <table_like>`.
    dimensions : int, optional
        Number of principal components to compute. Must be a positive integer less
        than or equal to min(n, p). Default is 3.
    max_iter : int, optional
        Maximum iterations for OptSpace algorithm. Default is 10000.

        .. versionchanged:: 0.7.4
            A non-positive or non-integer value now raises ``ValueError``
            instead of failing with an unrelated error partway through.
    sample_ids, feature_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

    Returns
    -------
    OrdinationResults
        The ordination results including sample coordinates, feature loadings,
        eigenvalues, and proportion explained.

    Raises
    ------
    ValueError
        If ``dimensions`` is not a positive integer less than or equal to
        min(n_samples, n_features)
    ValueError
        If ``max_iter`` is not a positive integer.
    ValueError
        If any row or column is fully unobserved.

    See Also
    --------
    pca
    pcoa
    optspace
    skbio.stats.composition.rclr

    Notes
    -----
    For compositional data with zero values, RPCA is most often preceded by the robust
    centered log-ratio (RCLR) transformation. See :func:`~skbio.stats.composition.rclr`
    for details. The entire procedure is referred to as robust Aitchison PCA.

    References
    ----------
    .. [1] Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
       Knight, R., & Zengler, K. (2019). A novel sparse compositional technique reveals
       microbial perturbations. mSystems, 4(1), 10-1128.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from skbio.stats.ordination import rpca

    Create a simple count table:

    >>> rng = np.random.default_rng(42)
    >>> counts = rng.poisson(5, size=(10, 20))
    >>> counts[counts < 2] = 0  # Add some zeros
    >>> table = pd.DataFrame(
    ...     counts,
    ...     index=[f'S{i}' for i in range(10)],
    ...     columns=[f'F{i}' for i in range(20)],
    ... )

    Perform RPCA:

    >>> ordination = rpca(table, dimensions=3)
    >>> print(ordination.proportion_explained[:3])  # doctest: +SKIP
    PC1    0.497464
    PC2    0.435214
    PC3    0.067322
    dtype: float64

    """

    # Ingestion of input data matrix
    X, row_ids, column_ids = _ingest_table(X, sample_ids, feature_ids)

    # Check observed entries
    obs_mask = ~np.isnan(X)
    obs_rows, obs_cols = np.any(obs_mask, axis=1), np.any(obs_mask, axis=0)

    if np.any(~obs_rows) or np.any(~obs_cols):
        raise ValueError(
            "Input to RPCA must not contain any fully unobserved rows or columns"
        )

    # Apply OptSpace for matrix completion
    X = optspace(X, dimensions=dimensions, max_iter=max_iter)

    # Perform PCA on the completed matrix
    output = _pca(X, method="svd", dimensions=dimensions)

    # Build the OrdinationResults object
    pc_ids = ["%s%d" % ("PC", i + 1) for i in range(output["variances"].shape[0])]
    eigvals = _create_table_1d(
        output["variances"],
        index=pc_ids,
        backend=output_format,
    )
    samples = _create_table(
        output["projected_samples"],
        index=row_ids,
        columns=pc_ids,
        backend=output_format,
    )
    features = _create_table(
        output["components"],
        index=column_ids,
        columns=pc_ids,
        backend=output_format,
    )

    return OrdinationResults(
        short_method_name="RPCA",
        long_method_name="Robust Principal Component Analysis",
        eigvals=eigvals,
        samples=samples,
        sample_ids=row_ids,
        features=features,
        feature_ids=column_ids,
        proportion_explained=eigvals / output["total_variance"],
    )
