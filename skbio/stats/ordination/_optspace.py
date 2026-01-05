# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""OptSpace Matrix Completion Algorithm.

This module provides the OptSpace algorithm for low-rank matrix completion
from partially observed entries. It is used by the Robust PCA (RPCA)
ordination method.

The algorithm minimizes the objective:

.. math::

    \min_{U, V, S} \|P_\Omega(Y - USV^T)\|_F^2

where :math:`Y` is the observed matrix, :math:`U` and :math:`V` are the
left and right singular vector matrices, :math:`S` is the diagonal matrix
of singular values, and :math:`P_\Omega` projects onto the observed entries.

References
----------
.. [1] Keshavan RH, Montanari A, Oh S. 2010. Matrix Completion from a
   Few Entries. IEEE Transactions on Information Theory 56(6):2980-2998.

.. [2] Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A,
   Knight R, Zengler K. 2019. A Novel Sparse Compositional Technique
   Reveals Microbial Perturbations. mSystems 4:e00016-19.

"""

import numpy as np
from scipy.sparse.linalg import svds
from scipy.linalg import svd


def _svd_sort(U, S, V):
    """Sort SVD components by descending singular values.

    Also applies sign correction for reproducibility.

    Parameters
    ----------
    U : ndarray
        Left singular vectors.
    S : ndarray
        Singular values.
    V : ndarray
        Right singular vectors.

    Returns
    -------
    tuple
        Sorted (U, S, V) with sign correction applied.
    """
    # Sort by descending singular values
    idx = np.argsort(S)[::-1]
    S = S[idx]
    U = U[:, idx]
    V = V[:, idx]

    # Apply sign correction for reproducibility
    # Make the largest element in each column of U positive
    for i in range(U.shape[1]):
        if np.abs(U[:, i].min()) > np.abs(U[:, i].max()):
            U[:, i] *= -1
            V[:, i] *= -1

    return U, S, V


def _compute_gradient(M, E, U, V, n, m, observed_mask):
    """Compute gradient for OptSpace optimization.

    Parameters
    ----------
    M : ndarray
        Reconstructed matrix.
    E : ndarray
        Error matrix.
    U : ndarray
        Left singular vectors.
    V : ndarray
        Right singular vectors.
    n : int
        Number of rows.
    m : int
        Number of columns.
    observed_mask : ndarray
        Binary mask indicating observed entries.

    Returns
    -------
    tuple
        Gradients for U and V.
    """
    # Apply mask to error
    E_masked = E * observed_mask

    # Compute gradients
    grad_U = E_masked.dot(V) / n
    grad_V = E_masked.T.dot(U) / m

    return grad_U, grad_V


def _line_search(U, V, S, grad_U, grad_V, M_observed, observed_mask,
                 n, m, step_size=1.0, max_iter=50, tol=1e-6):
    """Perform line search for step size optimization.

    Parameters
    ----------
    U, V, S : ndarray
        Current SVD components.
    grad_U, grad_V : ndarray
        Gradients for U and V.
    M_observed : ndarray
        Original observed matrix (with NaN for unobserved).
    observed_mask : ndarray
        Binary mask for observed entries.
    n, m : int
        Matrix dimensions.
    step_size : float
        Initial step size.
    max_iter : int
        Maximum line search iterations.
    tol : float
        Convergence tolerance.

    Returns
    -------
    tuple
        Updated (U, V, step_size).
    """
    # Current objective
    M = U.dot(S).dot(V.T)
    E = (M_observed - M) * observed_mask
    E = np.nan_to_num(E, nan=0.0)
    current_obj = np.sum(E ** 2)

    for _ in range(max_iter):
        # Trial update
        U_new = U - step_size * grad_U
        V_new = V - step_size * grad_V

        # Orthonormalize via QR
        U_new, _ = np.linalg.qr(U_new)
        V_new, _ = np.linalg.qr(V_new)

        # Compute new objective
        # First compute optimal S for new U, V
        S_new = _compute_singular_values(U_new, V_new, M_observed, observed_mask)

        M_new = U_new.dot(S_new).dot(V_new.T)
        E_new = (M_observed - M_new) * observed_mask
        E_new = np.nan_to_num(E_new, nan=0.0)
        new_obj = np.sum(E_new ** 2)

        if new_obj < current_obj:
            return U_new, V_new, S_new, step_size
        else:
            step_size *= 0.5

        if step_size < tol:
            break

    return U, V, S, step_size


def _compute_singular_values(U, V, M_observed, observed_mask):
    """Compute optimal singular values given U and V.

    Solves the least squares problem to find optimal S.

    Parameters
    ----------
    U : ndarray
        Left singular vectors (n x r).
    V : ndarray
        Right singular vectors (m x r).
    M_observed : ndarray
        Observed matrix values.
    observed_mask : ndarray
        Binary mask for observed entries.

    Returns
    -------
    ndarray
        Optimal diagonal singular value matrix (r x r).
    """
    r = U.shape[1]
    n, m = M_observed.shape

    # Build and solve the linear system for S
    # For each observed entry (i,j): M[i,j] = sum_k U[i,k] * S[k,k] * V[j,k]
    # This is a least squares problem

    # For diagonal S, we solve independently for each singular value
    S = np.zeros((r, r))

    for k in range(r):
        # Build the linear system for s_k
        uv_k = np.outer(U[:, k], V[:, k])
        uv_k_masked = uv_k * observed_mask

        # Flatten and solve
        a = uv_k_masked.flatten()
        b = (M_observed * observed_mask).flatten()
        b = np.nan_to_num(b, nan=0.0)

        # Filter to observed entries
        observed_flat = observed_mask.flatten() > 0
        if np.sum(observed_flat) > 0:
            # Least squares solution
            a_obs = a[observed_flat]
            b_obs = b[observed_flat]

            if np.sum(a_obs ** 2) > 1e-10:
                S[k, k] = np.dot(a_obs, b_obs) / np.dot(a_obs, a_obs)

    return S


class OptSpace:
    r"""Matrix completion using the OptSpace algorithm.

    OptSpace is an algorithm for recovering a low-rank matrix from a
    subset of observed entries. It uses gradient descent on the
    Grassmann manifold to find the optimal low-rank approximation.

    Parameters
    ----------
    n_components : int, optional
        The rank of the matrix to recover. Default is 3.
    max_iterations : int, optional
        Maximum number of iterations. Default is 5.
    tol : float, optional
        Convergence tolerance. Default is 1e-5.

    Attributes
    ----------
    U : ndarray
        Left singular vectors after fitting.
    S : ndarray
        Singular values (as diagonal matrix) after fitting.
    V : ndarray
        Right singular vectors after fitting.
    converged : bool
        Whether the algorithm converged.

    See Also
    --------
    rpca

    Notes
    -----
    The algorithm proceeds as follows:

    1. Initialize U, V using trimmed SVD of the observed matrix
    2. Iteratively:
       a. Compute optimal S given current U, V
       b. Compute gradients with respect to U, V
       c. Update U, V using gradient descent with line search
       d. Project U, V back to Grassmann manifold

    References
    ----------
    .. [1] Keshavan RH, Montanari A, Oh S. 2010. Matrix Completion from a
       Few Entries. IEEE Transactions on Information Theory 56(6):2980-2998.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import OptSpace
    >>> # Create a low-rank matrix with missing entries
    >>> np.random.seed(42)
    >>> true_U = np.random.randn(10, 2)
    >>> true_V = np.random.randn(8, 2)
    >>> true_M = true_U.dot(true_V.T)
    >>> # Mask some entries
    >>> M_observed = true_M.copy()
    >>> M_observed[::2, ::2] = np.nan  # Hide some entries
    >>> # Recover the matrix
    >>> opt = OptSpace(n_components=2, max_iterations=10)
    >>> M_recovered = opt.fit_transform(M_observed)

    """

    def __init__(self, n_components=3, max_iterations=5, tol=1e-5):
        self.n_components = n_components
        self.max_iterations = max_iterations
        self.tol = tol
        self.U = None
        self.S = None
        self.V = None
        self.converged = False

    def fit(self, X):
        """Fit the OptSpace model to the observed matrix.

        Parameters
        ----------
        X : ndarray
            A 2D array with observed values and NaN for missing entries.

        Returns
        -------
        self
            The fitted OptSpace instance.

        Raises
        ------
        ValueError
            If input is not 2D or n_components exceeds matrix dimensions.
        """
        X = np.asarray(X, dtype=np.float64)

        if X.ndim != 2:
            raise ValueError(f"Input must be 2D, got {X.ndim}D array.")

        n, m = X.shape
        r = self.n_components

        if r > min(n, m):
            raise ValueError(
                f"n_components ({r}) cannot exceed min matrix dimension "
                f"({min(n, m)})."
            )

        # Create observed mask (1 for observed, 0 for missing)
        observed_mask = ~np.isnan(X)
        observed_mask = observed_mask.astype(np.float64)

        # Replace NaN with 0 for computation
        X_filled = np.nan_to_num(X, nan=0.0)

        # Compute sparsity for rescaling
        n_observed = np.sum(observed_mask)
        sparsity = n_observed / (n * m)

        # Rescale observed values for sparse initialization
        X_scaled = X_filled / max(sparsity, 1e-10)

        # Initialize with truncated SVD
        try:
            if r < min(n, m) - 1:
                U, s, Vt = svds(X_scaled, k=r)
                V = Vt.T
            else:
                U, s, Vt = svd(X_scaled, full_matrices=False)
                U = U[:, :r]
                s = s[:r]
                V = Vt[:r, :].T
        except Exception:
            # Fallback to random initialization
            U = np.random.randn(n, r)
            V = np.random.randn(m, r)
            U, _ = np.linalg.qr(U)
            V, _ = np.linalg.qr(V)
            s = np.ones(r)

        # Sort by singular values
        U, s, V = _svd_sort(U, s, V)

        # Initialize S as diagonal matrix
        S = np.diag(s) * sparsity  # Scale back

        # Optimization loop
        step_size = 1.0
        prev_obj = np.inf

        for iteration in range(self.max_iterations):
            # Compute current reconstruction and error
            M = U.dot(S).dot(V.T)
            E = (X - M) * observed_mask
            E = np.nan_to_num(E, nan=0.0)

            # Current objective (Frobenius norm of error)
            obj = np.sum(E ** 2)

            # Check convergence
            if abs(prev_obj - obj) < self.tol:
                self.converged = True
                break

            prev_obj = obj

            # Compute gradients
            grad_U, grad_V = _compute_gradient(M, E, U, V, n, m, observed_mask)

            # Update with line search
            U, V, S, step_size = _line_search(
                U, V, S, grad_U, grad_V, X, observed_mask,
                n, m, step_size=step_size
            )

        # Recompute final S
        S = _compute_singular_values(U, V, X, observed_mask)

        # Final sort
        s_diag = np.diag(S)
        U, s_diag, V = _svd_sort(U, s_diag, V)
        S = np.diag(s_diag)

        self.U = U
        self.S = S
        self.V = V

        return self

    def transform(self, X=None):
        """Reconstruct the complete matrix.

        Parameters
        ----------
        X : ndarray, optional
            Not used, present for API compatibility.

        Returns
        -------
        ndarray
            The reconstructed low-rank matrix.

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.U is None:
            raise ValueError("Model has not been fitted. Call fit() first.")

        return self.U.dot(self.S).dot(self.V.T)

    def fit_transform(self, X):
        """Fit the model and return the reconstructed matrix.

        Parameters
        ----------
        X : ndarray
            A 2D array with observed values and NaN for missing entries.

        Returns
        -------
        ndarray
            The reconstructed low-rank matrix.
        """
        self.fit(X)
        return self.transform()

    def get_loadings(self):
        """Get sample and feature loadings.

        Returns
        -------
        tuple
            (sample_loadings, feature_loadings) where sample_loadings
            has shape (n_samples, n_components) and feature_loadings
            has shape (n_features, n_components).

        Raises
        ------
        ValueError
            If the model has not been fitted.
        """
        if self.U is None:
            raise ValueError("Model has not been fitted. Call fit() first.")

        s = np.sqrt(np.diag(self.S))
        sample_loadings = self.U * s
        feature_loadings = self.V * s

        return sample_loadings, feature_loadings
