# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Robust Centered Log-Ratio (rclr) Transform.

This module provides the robust centered log-ratio transformation for
compositional data analysis. The rclr transform is designed to handle
sparse data (i.e., data with many zeros) by only operating on observed
(non-zero) values.

References
----------
.. [1] Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A,
   Knight R, Zengler K. 2019. A Novel Sparse Compositional Technique
   Reveals Microbial Perturbations. mSystems 4:e00016-19.

"""

import numpy as np


def _matrix_closure(mat):
    """Perform matrix closure (normalize rows to sum to 1).

    Parameters
    ----------
    mat : ndarray
        A 2D array of non-negative values. May contain NaN for missing entries.

    Returns
    -------
    ndarray
        A 2D array where each row sums to 1. NaN values remain NaN.

    Notes
    -----
    Rows that sum to zero are left unchanged (as zeros).
    NaN values are excluded from the sum and remain NaN in the output.
    """
    mat = np.atleast_2d(mat).astype(np.float64)
    # Use nansum to handle NaN values (treating them as missing)
    row_sums = np.nansum(mat, axis=1, keepdims=True)
    # Avoid division by zero for empty rows
    row_sums[row_sums == 0] = 1.0
    return mat / row_sums


def matrix_rclr(mat):
    r"""Perform robust centered log-ratio transformation on a matrix.

    The robust centered log-ratio (rclr) transformation is similar to the
    standard CLR transformation, but it only operates on observed (non-zero)
    values. This makes it suitable for sparse compositional data such as
    microbiome count data.

    For each sample (row), the transformation computes:

    .. math::

        rclr(x_i) = \log(x_i) - \frac{1}{|S_i|} \sum_{j \in S_i} \log(x_j)

    where :math:`S_i` is the set of indices with non-zero values in sample
    :math:`i`, and :math:`|S_i|` is the number of non-zero values.

    Parameters
    ----------
    mat : array_like
        A 2D array of non-negative values (samples x features).
        Typically a count or abundance table.

    Returns
    -------
    ndarray
        A 2D array of the same shape with rclr-transformed values.
        Zero values in the input are replaced with NaN in the output.

    Raises
    ------
    ValueError
        If the input is not a 2D array, contains negative values,
        or contains non-finite values.

    See Also
    --------
    rpca
    ctf

    Notes
    -----
    The rclr transformation has several advantages for sparse compositional
    data:

    1. It does not require pseudocount addition, which can bias results
    2. It preserves the zero/non-zero structure of the data
    3. It allows for matrix completion methods (like OptSpace) to be applied

    The geometric mean is computed only over non-zero values in each row,
    making it "robust" to the presence of zeros.

    References
    ----------
    .. [1] Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A,
       Knight R, Zengler K. 2019. A Novel Sparse Compositional Technique
       Reveals Microbial Perturbations. mSystems 4:e00016-19.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import matrix_rclr
    >>> data = np.array([[1, 2, 0, 4],
    ...                  [0, 3, 3, 0],
    ...                  [2, 2, 2, 2]])
    >>> result = matrix_rclr(data)
    >>> print(np.round(result, 3))
    [[-0.693  0.     nan  0.693]
     [   nan  0.     0.     nan]
     [ 0.     0.     0.     0.   ]]

    """
    mat = np.atleast_2d(mat).astype(np.float64)

    if mat.ndim != 2:
        raise ValueError("Input must be a 2D array. Got array with "
                         f"{mat.ndim} dimensions.")

    # Check for negative values (excluding NaN which represents missing data)
    non_nan_mask = ~np.isnan(mat)
    if np.any(mat[non_nan_mask] < 0):
        raise ValueError("Input contains negative values. The rclr "
                         "transformation requires non-negative values.")

    # Check for infinities (NaN is allowed for missing entries)
    if np.any(np.isinf(mat)):
        raise ValueError("Input contains infinite values.")

    # Close the matrix (normalize rows to sum to 1)
    closed = _matrix_closure(mat)

    # Take natural log (zeros become -inf)
    with np.errstate(divide='ignore'):
        log_closed = np.log(closed)

    # Create output array initialized with NaN
    result = np.full_like(mat, np.nan)

    # For each row, compute rclr only for observed (non-zero) values
    for i in range(mat.shape[0]):
        # Find non-zero positions in original data
        observed_mask = mat[i] > 0
        if np.any(observed_mask):
            # Compute geometric mean over observed values only
            log_vals = log_closed[i, observed_mask]
            geo_mean_log = np.mean(log_vals)
            # Center by geometric mean
            result[i, observed_mask] = log_vals - geo_mean_log

    return result


def tensor_rclr(tensor):
    r"""Perform robust centered log-ratio transformation on a tensor.

    Extends the rclr transformation to N-dimensional arrays by treating
    the data as a collection of 2D matrices.

    Parameters
    ----------
    tensor : array_like
        An N-dimensional array of non-negative values where N >= 2.
        The first axis is treated as samples, and the last axis as features.
        Middle axes are treated as conditions/states.

    Returns
    -------
    ndarray
        An N-dimensional array of the same shape with rclr-transformed values.

    Raises
    ------
    ValueError
        If the input has fewer than 2 dimensions, contains negative values,
        or contains non-finite values.

    See Also
    --------
    matrix_rclr
    ctf

    Notes
    -----
    For a 3D tensor of shape (n_samples, n_conditions, n_features), the
    transformation is applied by first reshaping to 2D
    (n_samples * n_conditions, n_features), applying matrix_rclr,
    then reshaping back to 3D.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import tensor_rclr
    >>> # 3D tensor: 2 samples, 2 conditions, 3 features
    >>> data = np.array([[[1, 2, 3],
    ...                   [4, 5, 6]],
    ...                  [[7, 8, 9],
    ...                   [10, 11, 12]]])
    >>> result = tensor_rclr(data)

    """
    tensor = np.atleast_2d(tensor).astype(np.float64)

    if tensor.ndim < 2:
        raise ValueError("Input must have at least 2 dimensions. "
                         f"Got array with {tensor.ndim} dimensions.")

    # Check for negative values (excluding NaN)
    non_nan_mask = ~np.isnan(tensor)
    if np.any(tensor[non_nan_mask] < 0):
        raise ValueError("Input contains negative values. The rclr "
                         "transformation requires non-negative values.")

    # Check for infinities (NaN is allowed for missing entries)
    if np.any(np.isinf(tensor)):
        raise ValueError("Input contains infinite values.")

    # Store original shape
    original_shape = tensor.shape

    # Reshape to 2D: (n_samples * middle_dims, n_features)
    n_features = original_shape[-1]
    tensor_2d = tensor.reshape(-1, n_features)

    # Apply matrix_rclr
    result_2d = matrix_rclr(tensor_2d)

    # Reshape back to original shape
    return result_2d.reshape(original_shape)
