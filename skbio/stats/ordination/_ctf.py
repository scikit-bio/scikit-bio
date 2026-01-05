# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Compositional Tensor Factorization (CTF).

This module implements CTF for analyzing repeated-measures compositional
data. CTF extends RPCA to handle longitudinal or multi-condition studies
by constructing a three-way tensor and performing tensor factorization.

References
----------
.. [1] Martino C, Shenhav L, Marotz CA, et al. 2020. Context-aware
   dimensionality reduction deconvolutes gut microbial community
   dynamics. Nature Biotechnology 39:165-168.

"""

import numpy as np
import pandas as pd

from ._ordination_results import OrdinationResults
from ._rclr import tensor_rclr
from ._tf import TensorFactorization


def _build_tensor(table, sample_metadata, individual_id_column, state_column):
    """Build a 3-way tensor from a feature table and metadata.

    Parameters
    ----------
    table : pd.DataFrame
        Feature table (samples x features).
    sample_metadata : pd.DataFrame
        Sample metadata with individual_id and state columns.
    individual_id_column : str
        Column name for individual/subject identifiers.
    state_column : str
        Column name for condition/state identifiers.

    Returns
    -------
    tuple
        (tensor, individual_ids, state_ids, feature_ids) where tensor
        has shape (n_individuals, n_states, n_features).
    """
    # Ensure sample_metadata has table samples
    common_samples = table.index.intersection(sample_metadata.index)
    if len(common_samples) == 0:
        raise ValueError("No common samples between table and metadata.")

    table = table.loc[common_samples]
    sample_metadata = sample_metadata.loc[common_samples]

    # Get unique individuals and states
    individual_ids = sample_metadata[individual_id_column].unique().tolist()
    state_ids = sample_metadata[state_column].unique().tolist()
    feature_ids = table.columns.tolist()

    n_individuals = len(individual_ids)
    n_states = len(state_ids)
    n_features = len(feature_ids)

    # Build tensor with NaN for missing entries
    tensor = np.full((n_individuals, n_states, n_features), np.nan)

    # Create lookup dictionaries
    ind_to_idx = {ind: i for i, ind in enumerate(individual_ids)}
    state_to_idx = {state: i for i, state in enumerate(state_ids)}

    # Fill tensor
    for sample_id in common_samples:
        ind = sample_metadata.loc[sample_id, individual_id_column]
        state = sample_metadata.loc[sample_id, state_column]

        i = ind_to_idx[ind]
        j = state_to_idx[state]

        # Handle potential duplicates by averaging
        current_val = tensor[i, j, :]
        new_val = table.loc[sample_id].values

        if np.all(np.isnan(current_val)):
            tensor[i, j, :] = new_val
        else:
            # Average with existing value
            tensor[i, j, :] = (current_val + new_val) / 2

    return tensor, individual_ids, state_ids, feature_ids


def _filter_tensor(tensor, individual_ids, state_ids, feature_ids,
                   min_sample_count=0, min_feature_count=0,
                   min_feature_frequency=0.0):
    """Filter tensor by minimum counts and frequency.

    Parameters
    ----------
    tensor : ndarray
        3D tensor (individuals x states x features).
    individual_ids : list
        Individual identifiers.
    state_ids : list
        State/condition identifiers.
    feature_ids : list
        Feature identifiers.
    min_sample_count : int
        Minimum count per sample.
    min_feature_count : int
        Minimum total count per feature.
    min_feature_frequency : float
        Minimum frequency of feature occurrence.

    Returns
    -------
    tuple
        Filtered (tensor, individual_ids, state_ids, feature_ids).
    """
    tensor = tensor.copy()

    # Filter features by minimum count
    if min_feature_count > 0:
        feature_sums = np.nansum(tensor, axis=(0, 1))
        keep_features = feature_sums >= min_feature_count
        tensor = tensor[:, :, keep_features]
        feature_ids = [f for f, k in zip(feature_ids, keep_features) if k]

    # Filter features by frequency
    if min_feature_frequency > 0:
        n_samples = tensor.shape[0] * tensor.shape[1]
        feature_presence = np.sum(~np.isnan(tensor) & (tensor > 0), axis=(0, 1))
        feature_freq = feature_presence / n_samples
        keep_features = feature_freq >= min_feature_frequency
        tensor = tensor[:, :, keep_features]
        feature_ids = [f for f, k in zip(feature_ids, keep_features) if k]

    # Filter individuals with no observed data
    ind_has_data = np.any(~np.isnan(tensor), axis=(1, 2))
    tensor = tensor[ind_has_data, :, :]
    individual_ids = [i for i, k in zip(individual_ids, ind_has_data) if k]

    # Filter states with no observed data
    state_has_data = np.any(~np.isnan(tensor), axis=(0, 2))
    tensor = tensor[:, state_has_data, :]
    state_ids = [s for s, k in zip(state_ids, state_has_data) if k]

    return tensor, individual_ids, state_ids, feature_ids


def ctf(table, sample_metadata, individual_id_column, state_column,
        n_components=3, min_sample_count=0, min_feature_count=0,
        min_feature_frequency=0.0, max_als_iterations=25,
        max_rtpm_iterations=25, n_initializations=5):
    r"""Perform Compositional Tensor Factorization.

    CTF is an ordination method for repeated-measures compositional data.
    It extends RPCA to handle longitudinal or multi-condition studies by
    constructing a three-way tensor (subjects x conditions x features)
    and performing tensor factorization.

    Parameters
    ----------
    table : pd.DataFrame
        A feature table with samples as rows and features as columns.
        Values should be non-negative counts or abundances.
    sample_metadata : pd.DataFrame
        Metadata for each sample, indexed by sample ID. Must contain
        the individual_id_column and state_column.
    individual_id_column : str
        Name of the column in sample_metadata that identifies
        individual subjects.
    state_column : str
        Name of the column in sample_metadata that identifies
        conditions or time points.
    n_components : int, optional
        Number of components to compute. Default is 3.
    min_sample_count : int, optional
        Minimum total count per sample. Default is 0.
    min_feature_count : int, optional
        Minimum total count per feature. Default is 0.
    min_feature_frequency : float, optional
        Minimum proportion of samples a feature must appear in.
        Default is 0.
    max_als_iterations : int, optional
        Maximum ALS iterations. Default is 25.
    max_rtpm_iterations : int, optional
        Maximum RTPM iterations. Default is 25.
    n_initializations : int, optional
        Number of random initializations. Default is 5.

    Returns
    -------
    OrdinationResults
        Subject ordination results with subject coordinates and
        feature loadings.
    OrdinationResults
        State ordination results with state coordinates and
        feature loadings.

    Raises
    ------
    ValueError
        If input validation fails or insufficient data after filtering.

    See Also
    --------
    rpca
    TensorFactorization

    Notes
    -----
    CTF is designed for repeated-measures studies where each subject
    is sampled multiple times across different conditions or time
    points. The method:

    1. Constructs a 3-way tensor (subjects x conditions x features)
    2. Applies tensor rclr transformation
    3. Performs CP tensor decomposition to find low-rank structure
    4. Returns separate ordinations for subjects and conditions

    References
    ----------
    .. [1] Martino C, Shenhav L, Marotz CA, et al. 2020. Context-aware
       dimensionality reduction deconvolutes gut microbial community
       dynamics. Nature Biotechnology 39:165-168.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from skbio.stats.ordination import ctf

    Create a feature table and metadata for a longitudinal study:

    >>> np.random.seed(42)
    >>> n_samples = 20
    >>> n_features = 15
    >>> counts = np.random.poisson(10, size=(n_samples, n_features))
    >>> counts[counts < 2] = 0
    >>> sample_ids = ['s%d' % i for i in range(n_samples)]
    >>> table = pd.DataFrame(counts, index=sample_ids,
    ...                      columns=['f%d' % i for i in range(n_features)])

    Create metadata with 5 subjects measured at 4 time points:

    >>> metadata = pd.DataFrame({
    ...     'subject': ['subj_%d' % (i // 4) for i in range(n_samples)],
    ...     'timepoint': ['t%d' % (i % 4) for i in range(n_samples)]
    ... }, index=sample_ids)

    Perform CTF:

    >>> subject_ord, state_ord = ctf(
    ...     table, metadata, 'subject', 'timepoint', n_components=2
    ... )  # doctest: +SKIP

    """
    # Validate inputs
    if not isinstance(table, pd.DataFrame):
        raise ValueError("Input table must be a pandas DataFrame.")

    if not isinstance(sample_metadata, pd.DataFrame):
        raise ValueError("sample_metadata must be a pandas DataFrame.")

    if individual_id_column not in sample_metadata.columns:
        raise ValueError(
            f"Column '{individual_id_column}' not found in sample_metadata."
        )

    if state_column not in sample_metadata.columns:
        raise ValueError(
            f"Column '{state_column}' not found in sample_metadata."
        )

    if (table.values < 0).any():
        raise ValueError("Input table contains negative values.")

    # Build tensor from table and metadata
    tensor, individual_ids, state_ids, feature_ids = _build_tensor(
        table, sample_metadata, individual_id_column, state_column
    )

    # Filter tensor
    tensor, individual_ids, state_ids, feature_ids = _filter_tensor(
        tensor, individual_ids, state_ids, feature_ids,
        min_sample_count, min_feature_count, min_feature_frequency
    )

    # Validate dimensions
    if len(individual_ids) < 3:
        raise ValueError(
            f"Only {len(individual_ids)} individuals remain after filtering. "
            "At least 3 are required."
        )

    if len(state_ids) < 2:
        raise ValueError(
            f"Only {len(state_ids)} states remain after filtering. "
            "At least 2 are required."
        )

    if len(feature_ids) < n_components:
        raise ValueError(
            f"Only {len(feature_ids)} features remain after filtering. "
            f"This is less than n_components ({n_components})."
        )

    # Apply tensor rclr transformation
    rclr_tensor = tensor_rclr(tensor)

    # Perform tensor factorization
    tf = TensorFactorization(
        n_components=n_components,
        max_als_iterations=max_als_iterations,
        max_rtpm_iterations=max_rtpm_iterations,
        n_initializations=n_initializations
    )
    tf.fit(rclr_tensor)

    # Extract factor matrices
    subject_loadings = tf.get_subject_loadings()
    state_loadings = tf.get_condition_loadings()
    feature_loadings = tf.get_feature_loadings()

    # Create axis labels
    axis_labels = ['PC%d' % i for i in range(1, n_components + 1)]

    # Compute eigenvalue-like quantities from factor norms
    # These represent the variance captured by each component
    factor_norms = np.array([
        np.linalg.norm(subject_loadings[:, i]) *
        np.linalg.norm(state_loadings[:, i]) *
        np.linalg.norm(feature_loadings[:, i])
        for i in range(n_components)
    ])

    # Normalize to get proportion explained (approximation)
    proportion_explained = factor_norms / np.sum(factor_norms)
    eigvals = factor_norms ** 2

    # Create subject OrdinationResults
    subject_ordination = OrdinationResults(
        short_method_name='CTF',
        long_method_name='Compositional Tensor Factorization (Subjects)',
        eigvals=pd.Series(eigvals, index=axis_labels),
        samples=pd.DataFrame(subject_loadings, index=individual_ids,
                             columns=axis_labels),
        features=pd.DataFrame(feature_loadings, index=feature_ids,
                              columns=axis_labels),
        proportion_explained=pd.Series(proportion_explained,
                                        index=axis_labels)
    )

    # Create state OrdinationResults
    state_ordination = OrdinationResults(
        short_method_name='CTF',
        long_method_name='Compositional Tensor Factorization (States)',
        eigvals=pd.Series(eigvals, index=axis_labels),
        samples=pd.DataFrame(state_loadings, index=state_ids,
                             columns=axis_labels),
        features=pd.DataFrame(feature_loadings, index=feature_ids,
                              columns=axis_labels),
        proportion_explained=pd.Series(proportion_explained,
                                        index=axis_labels)
    )

    return subject_ordination, state_ordination
