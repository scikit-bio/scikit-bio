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
from ._rclr import matrix_rclr
from ._optspace import OptSpace


def _filter_table(table, min_sample_count=0, min_feature_count=0,
                  min_feature_frequency=0.0):
    """Filter a feature table by minimum counts and frequency.

    Parameters
    ----------
    table : pd.DataFrame
        Feature table (samples x features).
    min_sample_count : int
        Minimum total count per sample.
    min_feature_count : int
        Minimum total count per feature.
    min_feature_frequency : float
        Minimum proportion of samples a feature must appear in.

    Returns
    -------
    pd.DataFrame
        Filtered table.
    """
    table = table.copy()

    # Filter by minimum sample count
    if min_sample_count > 0:
        sample_sums = table.sum(axis=1)
        table = table.loc[sample_sums >= min_sample_count]

    # Filter by minimum feature count
    if min_feature_count > 0:
        feature_sums = table.sum(axis=0)
        table = table.loc[:, feature_sums >= min_feature_count]

    # Filter by minimum feature frequency
    if min_feature_frequency > 0:
        n_samples = table.shape[0]
        feature_freq = (table > 0).sum(axis=0) / n_samples
        table = table.loc[:, feature_freq >= min_feature_frequency]

    return table


def rpca(table, n_components=3, min_sample_count=0, min_feature_count=0,
         min_feature_frequency=0.0, max_iterations=5):
    r"""Perform Robust Principal Component Analysis.

    Robust PCA (RPCA) is an ordination method for sparse compositional data.
    It applies the robust centered log-ratio (rclr) transformation followed
    by OptSpace matrix completion to handle the zeros introduced by the
    log transformation.

    Parameters
    ----------
    table : pd.DataFrame
        A feature table with samples as rows and features as columns.
        Values should be non-negative counts or abundances.
    n_components : int, optional
        Number of principal components to compute. Default is 3.
    min_sample_count : int, optional
        Minimum total count per sample. Samples with lower counts
        are removed. Default is 0 (no filtering).
    min_feature_count : int, optional
        Minimum total count per feature. Features with lower counts
        are removed. Default is 0 (no filtering).
    min_feature_frequency : float, optional
        Minimum proportion of samples a feature must appear in.
        Features appearing in fewer samples are removed. Default is 0.
    max_iterations : int, optional
        Maximum iterations for OptSpace algorithm. Default is 5.

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
    ctf
    pcoa
    matrix_rclr
    OptSpace

    Notes
    -----
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

    >>> ordination = rpca(table, n_components=3)
    >>> print(ordination.proportion_explained[:3])  # doctest: +SKIP
    PC1    0.35...
    PC2    0.20...
    PC3    0.15...
    dtype: float64

    """
    # Validate input
    if not isinstance(table, pd.DataFrame):
        raise ValueError("Input table must be a pandas DataFrame.")

    if (table.values < 0).any():
        raise ValueError("Input table contains negative values.")

    # Filter table
    table = _filter_table(table, min_sample_count, min_feature_count,
                          min_feature_frequency)

    if table.shape[0] < 3:
        raise ValueError(
            f"After filtering, only {table.shape[0]} samples remain. "
            "At least 3 samples are required."
        )

    if table.shape[1] < n_components:
        raise ValueError(
            f"After filtering, only {table.shape[1]} features remain. "
            f"This is less than n_components ({n_components})."
        )

    # Get sample and feature IDs
    sample_ids = table.index.tolist()
    feature_ids = table.columns.tolist()

    # Apply rclr transformation
    rclr_table = matrix_rclr(table.values)

    # Apply OptSpace for matrix completion
    optspace = OptSpace(n_components=n_components,
                        max_iterations=max_iterations)
    optspace.fit(rclr_table)

    # Get the completed matrix
    completed = optspace.transform()

    # Center the completed matrix
    completed_centered = completed - np.nanmean(completed, axis=0)

    # Perform SVD on the completed centered matrix
    U, s, Vt = svd(completed_centered, full_matrices=False)

    # Truncate to n_components
    U = U[:, :n_components]
    s = s[:n_components]
    Vt = Vt[:n_components, :]

    # Compute eigenvalues (squared singular values)
    eigvals = s ** 2

    # Compute proportion explained
    total_variance = np.sum(svd(completed_centered, compute_uv=False) ** 2)
    proportion_explained = eigvals / total_variance

    # Compute sample coordinates (scores)
    sample_coordinates = U * s

    # Compute feature loadings
    feature_loadings = Vt.T * s

    # Create axis labels
    axis_labels = ['PC%d' % i for i in range(1, n_components + 1)]

    # Create OrdinationResults
    ordination_results = OrdinationResults(
        short_method_name='RPCA',
        long_method_name='Robust Principal Component Analysis',
        eigvals=pd.Series(eigvals, index=axis_labels),
        samples=pd.DataFrame(sample_coordinates, index=sample_ids,
                             columns=axis_labels),
        features=pd.DataFrame(feature_loadings, index=feature_ids,
                              columns=axis_labels),
        proportion_explained=pd.Series(proportion_explained,
                                        index=axis_labels)
    )

    return ordination_results
