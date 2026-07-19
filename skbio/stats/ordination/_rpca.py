# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Robust Principal Component Analysis (RPCA).

This module implements Robust PCA for sparse compositional data analysis.
RPCA combines the robust centered log-ratio (rclr) transformation with
OptSpace matrix completion to perform dimensionality reduction on
compositional data with many zeros.

References
----------
.. [1] Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A,
   Knight R, Zengler K. 2019. A Novel Sparse Compositional Technique
   Reveals Microbial Perturbations. mSystems 4:e00016-19.

"""

import numpy as np
import pandas as pd
from scipy.linalg import svd

from ._ordination_results import OrdinationResults
from skbio.stats.composition._base import rclr
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
    sample_ids=None,
    feature_ids=None,
    max_iter=10000,
    output_format=None,
):
    r"""Perform Robust Principal Component Analysis.

    \/ Rewrite this, and mention the application of rclr for preprocessing
    link to rclr

    Robust PCA (RPCA) is an ordination method for sparse compositional data.
    It applies the robust centered log-ratio (rclr) transformation followed
    by OptSpace matrix completion to handle the zeros introduced by the
    log transformation.

    Parameters
    ----------
    X : table_like
        Samples by features table (n, p). See :ref:`supported formats <table_like>`.
        Values should be non-negative counts or abundances.
    dimensions : int, optional
            dimensions : int, optional
        Number of principal components to compute. Must be a positive integer less
        than or equal to min(n, p). Default is 3.
    max_iter : int, optional
        Maximum iterations for OptSpace algorithm. Default is 5.
    sample_ids, feature_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

    Returns
    -------
    OrdinationResults
        The ordination results including sample coordinates, feature
        loadings, eigenvalues, and proportion explained.

    Raises
    ------
    ValueError
        If input is not a DataFrame, contains negative values,
        or has insufficient samples/features after filtering.

    See Also
    --------
    pca
    ctf
    pcoa
    rclr
    OptSpace

    Notes
    -----

    Show how to preprocess data with rclr, and show how to use with

    RPCA is designed for cross-sectional studies where each sample
    represents an independent observation. For repeated-measures or
    longitudinal data, consider using CTF instead.

    The method proceeds as follows:

    1. Filter the table by minimum counts and frequency (optional)
    2. Apply robust CLR transformation (log-ratio with geometric mean)
    3. Use OptSpace to complete the matrix (fill NaN values)
    4. Perform SVD on the completed matrix to extract principal components

    References
    ----------
    \/ Add citation for RPCA in general

    .. [1] Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A,
       Knight R, Zengler K. 2019. A Novel Sparse Compositional Technique
       Reveals Microbial Perturbations. mSystems 4:e00016-19.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from skbio.stats.ordination import rpca

    Create a simple count table:

    >>> np.random.seed(42)
    >>> counts = np.random.poisson(5, size=(10, 20))
    >>> counts[counts < 2] = 0  # Add some zeros
    >>> table = pd.DataFrame(counts,
    ...                      index=['sample_%d' % i for i in range(10)],
    ...                      columns=['feature_%d' % i for i in range(20)])

    Perform RPCA:

    >>> ordination = rpca(table, dimensions=3)
    >>> print(ordination.proportion_explained[:3])  # doctest: +SKIP
    PC1    0.35...
    PC2    0.20...
    PC3    0.15...
    dtype: float64

    """

    # Ingestion of input data matrix
    X, row_ids, column_ids = _ingest_table(X, sample_ids, feature_ids)

    # Note:
    # Input validation may be redundant here, since optspace already performs this

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
