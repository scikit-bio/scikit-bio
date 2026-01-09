# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Tensor Factorization for Compositional Data.

This module provides tensor factorization methods for analyzing
multi-way compositional data, such as longitudinal microbiome studies.
The implementation uses CANDECOMP/PARAFAC (CP) decomposition with
Alternating Least Squares (ALS) optimization.

References
----------
.. [1] Martino C, Shenhav L, Marotz CA, et al. 2021. Context-aware
   dimensionality reduction deconvolutes gut microbial community
   dynamics. Nature Biotechnology 39:165-168.

.. [2] Kolda TG, Bader BW. 2009. Tensor Decompositions and Applications.
   SIAM Review 51(3):455-500.

"""

import numpy as np
from scipy.linalg import svd


def _khatri_rao(matrices):
    """Compute the Khatri-Rao product of a list of matrices.

    The Khatri-Rao product is the column-wise Kronecker product.

    Parameters
    ----------
    matrices : list of ndarray
        List of 2D arrays with the same number of columns.

    Returns
    -------
    ndarray
        Khatri-Rao product with shape (prod(n_rows), n_cols).
    """
    if len(matrices) == 0:
        raise ValueError("At least one matrix is required.")

    result = matrices[0]
    for mat in matrices[1:]:
        n1, r = result.shape
        n2 = mat.shape[0]
        new_result = np.zeros((n1 * n2, r))
        for j in range(r):
            new_result[:, j] = np.outer(result[:, j], mat[:, j]).flatten()
        result = new_result

    return result


def _unfold_tensor(tensor, mode):
    """Unfold (matricize) a tensor along a specified mode.

    Parameters
    ----------
    tensor : ndarray
        N-dimensional array.
    mode : int
        The mode along which to unfold (0-indexed).

    Returns
    -------
    ndarray
        2D matrix unfolding of the tensor.
    """
    ndim = tensor.ndim
    dims = list(range(ndim))
    dims.remove(mode)
    dims = [mode] + dims

    # Permute and reshape
    permuted = np.transpose(tensor, dims)
    return permuted.reshape(tensor.shape[mode], -1)


def _refold_tensor(matrix, mode, shape):
    """Refold a matrix back into a tensor.

    Parameters
    ----------
    matrix : ndarray
        2D matrix (unfolded tensor).
    mode : int
        The mode along which the tensor was unfolded.
    shape : tuple
        Target tensor shape.

    Returns
    -------
    ndarray
        Refolded tensor.
    """
    ndim = len(shape)
    dims = list(range(ndim))
    dims.remove(mode)
    dims = [mode] + dims

    # Compute the permuted shape
    permuted_shape = [shape[mode]] + [shape[d] for d in dims[1:]]

    # Reshape and inverse permute
    tensor_permuted = matrix.reshape(permuted_shape)

    # Compute inverse permutation
    inverse_dims = [0] * ndim
    for i, d in enumerate(dims):
        inverse_dims[d] = i

    return np.transpose(tensor_permuted, inverse_dims)


def _construct_tensor_from_cp(factors):
    """Construct a tensor from CP decomposition factors.

    Parameters
    ----------
    factors : list of ndarray
        List of factor matrices, one for each mode.

    Returns
    -------
    ndarray
        Reconstructed tensor.
    """
    n_modes = len(factors)
    rank = factors[0].shape[1]
    shape = tuple(f.shape[0] for f in factors)

    # Initialize tensor
    tensor = np.zeros(shape)

    # Sum outer products
    for r in range(rank):
        component = factors[0][:, r]
        for mode in range(1, n_modes):
            component = np.outer(component, factors[mode][:, r]).reshape(
                component.shape + (factors[mode].shape[0],)
            )
        tensor += component

    return tensor


def _rtpm_update(tensor, factors, mode, observed_mask):
    """Perform one Robust Tensor Power Method update for a mode.

    Parameters
    ----------
    tensor : ndarray
        The observed tensor with NaN for missing entries.
    factors : list of ndarray
        Current factor matrices for all modes.
    mode : int
        Mode to update.
    observed_mask : ndarray
        Binary mask indicating observed entries.

    Returns
    -------
    ndarray
        Updated factor matrix for the specified mode.
    """
    n_modes = len(factors)
    rank = factors[0].shape[1]

    # Get the unfolded tensor
    tensor_unfolded = _unfold_tensor(tensor, mode)
    mask_unfolded = _unfold_tensor(observed_mask, mode)

    # Replace NaN with 0 for computation
    tensor_unfolded = np.nan_to_num(tensor_unfolded, nan=0.0)

    # Compute Khatri-Rao product of all other factors
    other_factors = [factors[m] for m in range(n_modes) if m != mode]
    if len(other_factors) > 1:
        kr_product = _khatri_rao(other_factors[::-1])  # Reverse order for proper KR
    else:
        kr_product = other_factors[0]

    # Update factor with least squares
    new_factor = np.zeros((tensor.shape[mode], rank))

    for i in range(tensor.shape[mode]):
        # Get row of unfolded tensor and mask
        y_i = tensor_unfolded[i, :]
        m_i = mask_unfolded[i, :]

        # Only use observed entries
        observed_idx = m_i > 0

        if np.sum(observed_idx) > 0:
            A = kr_product[observed_idx, :]
            b = y_i[observed_idx]

            # Regularized least squares
            AtA = A.T.dot(A)
            Atb = A.T.dot(b)

            # Add small regularization
            reg = 1e-6 * np.eye(rank)
            try:
                new_factor[i, :] = np.linalg.solve(AtA + reg, Atb)
            except np.linalg.LinAlgError:
                # Fallback to least squares
                new_factor[i, :] = np.linalg.lstsq(A, b, rcond=None)[0]

    return new_factor


def _initialize_factors(tensor_shape, rank, observed_mask, tensor, method='svd'):
    """Initialize factor matrices for tensor factorization.

    Parameters
    ----------
    tensor_shape : tuple
        Shape of the tensor.
    rank : int
        Rank of the decomposition.
    observed_mask : ndarray
        Binary mask for observed entries.
    tensor : ndarray
        The observed tensor (with NaN for missing).
    method : str
        Initialization method: 'svd' or 'random'.

    Returns
    -------
    list of ndarray
        Initialized factor matrices.
    """
    n_modes = len(tensor_shape)
    factors = []

    # Fill NaN with zeros for initialization
    tensor_filled = np.nan_to_num(tensor, nan=0.0)

    for mode in range(n_modes):
        if method == 'svd' and tensor_shape[mode] > rank:
            # SVD-based initialization
            unfolded = _unfold_tensor(tensor_filled, mode)
            try:
                U, s, Vt = svd(unfolded, full_matrices=False)
                factor = U[:, :rank] * np.sqrt(s[:rank])
            except Exception:
                factor = np.random.randn(tensor_shape[mode], rank)
        else:
            # Random initialization
            factor = np.random.randn(tensor_shape[mode], rank)

        # Normalize columns
        norms = np.linalg.norm(factor, axis=0, keepdims=True)
        norms[norms == 0] = 1
        factor = factor / norms

        factors.append(factor)

    return factors


class TensorFactorization:
    r"""Tensor factorization using CP decomposition with ALS.

    This class implements CANDECOMP/PARAFAC (CP) tensor decomposition
    using Alternating Least Squares (ALS) optimization. It is designed
    to handle tensors with missing values, making it suitable for
    analyzing incomplete multi-way data.

    Parameters
    ----------
    n_components : int, optional
        The rank of the decomposition. Default is 3.
    max_als_iterations : int, optional
        Maximum number of ALS iterations. Default is 25.
    max_rtpm_iterations : int, optional
        Maximum iterations for Robust Tensor Power Method
        initialization. Default is 25.
    n_initializations : int, optional
        Number of random initializations to try. Default is 5.
    tol : float, optional
        Convergence tolerance. Default is 1e-7.
    center : bool, optional
        Whether to center the factors after fitting. Default is True.

    Attributes
    ----------
    factors : list of ndarray
        The factor matrices after fitting, one per tensor mode.
    converged : bool
        Whether the algorithm converged.
    reconstruction_error : float
        Final reconstruction error.

    See Also
    --------
    ctf

    Notes
    -----
    CP decomposition approximates a tensor as a sum of rank-1 tensors:

    .. math::

        \mathcal{X} \approx \sum_{r=1}^{R} \mathbf{a}_r \circ \mathbf{b}_r
        \circ \mathbf{c}_r

    where :math:`\circ` denotes the outer product.

    References
    ----------
    .. [1] Kolda TG, Bader BW. 2009. Tensor Decompositions and Applications.
       SIAM Review 51(3):455-500.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import TensorFactorization
    >>> # Create a 3-way tensor
    >>> np.random.seed(42)
    >>> tensor = np.random.randn(5, 4, 3)
    >>> # Add some missing values
    >>> tensor[::2, ::2, 0] = np.nan
    >>> # Fit the model
    >>> tf = TensorFactorization(n_components=2, max_als_iterations=50)
    >>> tf.fit(tensor)  # doctest: +SKIP
    >>> # Get factor matrices
    >>> subjects, conditions, features = tf.factors  # doctest: +SKIP

    """

    def __init__(self, n_components=3, max_als_iterations=25,
                 max_rtpm_iterations=25, n_initializations=5,
                 tol=1e-7, center=True):
        self.n_components = n_components
        self.max_als_iterations = max_als_iterations
        self.max_rtpm_iterations = max_rtpm_iterations
        self.n_initializations = n_initializations
        self.tol = tol
        self.center = center

        self.factors = None
        self.converged = False
        self.reconstruction_error = np.inf

    def fit(self, X):
        """Fit the tensor factorization model.

        Parameters
        ----------
        X : ndarray
            An N-dimensional array (N >= 3) with NaN for missing entries.

        Returns
        -------
        self
            The fitted TensorFactorization instance.

        Raises
        ------
        ValueError
            If input has fewer than 3 dimensions or n_components
            exceeds tensor dimensions.
        """
        X = np.asarray(X, dtype=np.float64)

        if X.ndim < 3:
            raise ValueError(
                f"Input must have at least 3 dimensions, got {X.ndim}."
            )

        if not np.any(np.isnan(X)):
            raise ValueError(
                "Input must contain missing values (NaN). For complete "
                "tensors, use standard CP decomposition methods."
            )

        tensor_shape = X.shape
        rank = self.n_components

        if rank > min(tensor_shape):
            raise ValueError(
                f"n_components ({rank}) cannot exceed minimum tensor "
                f"dimension ({min(tensor_shape)})."
            )

        # Create observed mask
        observed_mask = ~np.isnan(X)
        observed_mask = observed_mask.astype(np.float64)

        best_factors = None
        best_error = np.inf

        # Try multiple initializations
        for init_idx in range(self.n_initializations):
            # Initialize factors
            factors = _initialize_factors(
                tensor_shape, rank, observed_mask, X,
                method='svd' if init_idx == 0 else 'random'
            )

            # ALS iterations
            prev_error = np.inf

            for als_iter in range(self.max_als_iterations):
                # Update each mode
                for mode in range(X.ndim):
                    factors[mode] = _rtpm_update(
                        X, factors, mode, observed_mask
                    )

                # Compute reconstruction error
                reconstructed = _construct_tensor_from_cp(factors)
                error_tensor = (X - reconstructed) * observed_mask
                error_tensor = np.nan_to_num(error_tensor, nan=0.0)
                current_error = np.sqrt(np.sum(error_tensor ** 2))

                # Check convergence
                if abs(prev_error - current_error) < self.tol:
                    break

                prev_error = current_error

            if current_error < best_error:
                best_error = current_error
                best_factors = [f.copy() for f in factors]

        self.factors = best_factors
        self.reconstruction_error = best_error
        self.converged = True

        # Center factors if requested
        if self.center:
            for i in range(len(self.factors)):
                self.factors[i] = self.factors[i] - np.mean(
                    self.factors[i], axis=0, keepdims=True
                )

        return self

    def transform(self, X=None):
        """Reconstruct the complete tensor.

        Parameters
        ----------
        X : ndarray, optional
            Not used, present for API compatibility.

        Returns
        -------
        ndarray
            The reconstructed tensor.

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.factors is None:
            raise ValueError("Model has not been fitted. Call fit() first.")

        return _construct_tensor_from_cp(self.factors)

    def fit_transform(self, X):
        """Fit the model and return the reconstructed tensor.

        Parameters
        ----------
        X : ndarray
            An N-dimensional array with NaN for missing entries.

        Returns
        -------
        ndarray
            The reconstructed tensor.
        """
        self.fit(X)
        return self.transform()

    def get_subject_loadings(self):
        """Get subject (first mode) loadings.

        Returns
        -------
        ndarray
            Subject loadings of shape (n_subjects, n_components).

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.factors is None:
            raise ValueError("Model has not been fitted. Call fit() first.")
        return self.factors[0]

    def get_condition_loadings(self):
        """Get condition/state (second mode) loadings.

        Returns
        -------
        ndarray
            Condition loadings of shape (n_conditions, n_components).

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.factors is None:
            raise ValueError("Model has not been fitted. Call fit() first.")
        return self.factors[1]

    def get_feature_loadings(self):
        """Get feature (last mode) loadings.

        Returns
        -------
        ndarray
            Feature loadings of shape (n_features, n_components).

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.factors is None:
            raise ValueError("Model has not been fitted. Call fit() first.")
        return self.factors[-1]
