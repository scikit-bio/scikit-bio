# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""OptSpace Matrix Completion Algorithm.

This module provides the OptSpace algorithm for low-rank matrix completion
from partially observed entries. It is used by the Robust PCA (RPCA)
ordination method.

The algorithm minimizes the objective:

.. math::

    \min_{U, V, S} \|P_\Omega(M_{obs} - USV^T)\|_F^2

where :math:`M_{obs}` is the observed matrix, :math:`U` and :math:`V` are the
left and right singular vector matrices, :math:`S` is the diagonal matrix
of singular values, and :math:`P_\Omega` projects onto the observed entries.

References
----------
.. [1] Keshavan, R. H., Montanari, A., & Oh, S. (2010). Matrix completion from a
   few entries. IEEE transactions on information theory, 56(6), 2980-2998.

.. [2] Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
   Knight, R., & Zengler, K. (2019). A novel sparse compositional technique reveals
   microbial perturbations. mSystems, 4(1), 10-1128.

"""

from warnings import warn

import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import svds, lsmr, LinearOperator


def _check_unobserved(observed_mask):

    # Check for fully unobserved rows or columns
    keep_rows = np.any(observed_mask, axis=1)
    keep_cols = np.any(observed_mask, axis=0)

    # Check that input is not completely unobserved
    if np.all(~keep_rows):
        raise ValueError(
            "Matrix completion requires at least one observed row "
            "and one observed column."
        )

    # Warn if any row or column is fully unobserved
    if np.any(~keep_rows) or np.any(~keep_cols):
        warn(
            "Input contains rows or columns which are fully unobserved. "
            "These rows or columns will be ignored for matrix completion. "
            "Fully unobserved rows or columns will remain as NaN in the output.",
            RuntimeWarning,
        )

    return keep_rows, keep_cols


def _trim(X, observed_mask, m, n, n_observed):
    """Trim over-represented rows and columns.

    Any row or column with more than half the average observed entries per
    row or column respectively is set to zero per Keshavan et al. (2010).
    This makes the low-rank structure of the observed data more
    prominent."""

    n_observed_rows = np.sum(observed_mask, axis=1)
    n_observed_cols = np.sum(observed_mask, axis=0)

    row_threshold = 2 * n_observed / m
    col_threshold = 2 * n_observed / n

    valid_rows = n_observed_rows <= row_threshold
    valid_cols = n_observed_cols <= col_threshold

    trim_mask = np.outer(valid_rows, valid_cols)

    return np.where(trim_mask & observed_mask, X, 0.0)


def _svd_init(X_trimmed, r):
    """Initialize U and V with an SVD of the trimmed observation matrix."""

    # Attempt sparse SVDS
    try:
        # Sparse SVDS
        U, s, Vt = svds(X_trimmed, k=r, solver="propack")
        V = Vt.T

        # Sort unsorted singular values from SVDS
        idx = np.argsort(s)[::-1]
        U = U[:, idx]
        V = Vt[idx, :].T

    # Sparse SVDS may fail to converge for sparse matrices, in which
    # case we fall back to dense SVD
    except Exception as e:
        warn(
            "Sparse SVD failed. Falling back to dense SVD instead.",
            RuntimeWarning,
        )

        # Dense SVD
        U, _, Vt = svd(X_trimmed, full_matrices=False)
        U = U[:, :r]
        V = Vt[:r, :].T

    return U, V


def jacobian_S(U, V, S, rows, cols):
    """Compute J_S(dS).

    The Jacobian is

    J(dU, dV, dS) = J_S(dS) + J_U(dU) + J_V(dV).

    The S component of the Jacobian determines how changes in S (dS) contribute to
    the reconstruction error over the observed entries.

    J_S(dS) = P_\\Omega(U dS V^T).
    """

    W = U @ S @ V.T
    return W[rows, cols]


def jacobian_S_adj(U, V, w, observed_mask):
    """Compute J_S*(W).

    The Jacobian adjoint is defined with respect to the inner product by

    <J(dU, dV, dS), W> = <(dU, dV, dS), J*(W)>

    Thus,

    J_S*(W) = U^T P_\\Omega(W) V.
    """

    W = np.zeros_like(observed_mask, dtype=U.dtype)
    W[observed_mask] = w
    ds = U.T @ W @ V
    return ds.ravel()


def _solve_S(U, V, b, observed_mask, tol):
    """Compute optimal S given U and V.

    Solves the least squares problem to find the optimal S that
    minimizes the reconstruction difference on the observed entries:

    arg min_S ||P_\\Omega(U V S^T - M_observed)||_F^2

    where P_\\Omega is the projection onto the observed entries.
    This is the least-squares solution to the system

    J_S(dS) = -R

    This is solved via lsmr.
    """

    r = U.shape[1]
    n_observed = np.sum(observed_mask)
    rows, cols = np.where(observed_mask)

    def matvec(s):
        return jacobian_S(U, V, s.reshape(r, r), rows, cols)

    def rmatvec(w):
        return jacobian_S_adj(U, V, w, observed_mask)

    J_S = LinearOperator(
        shape=(n_observed, r**2),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=U.dtype,
    )

    s = lsmr(J_S, b, atol=tol, btol=tol)[0]

    return s.reshape(r, r)


def jacobian_UV(U, V, S, dU, dV, rows, cols):
    """Compute J_UV(dU, dV).

    The Jacobian is

    J(dU, dV, dS) = J_S(dS) + J_U(dU) + J_V(dV),

    where J_UV(dU, dV) = J_U(dU) + J_V(dV). The UV component of the Jacobian
    determines how changes in U and V (dU, dV) contribute to the reconstruction
    error over the observed entries.

    J_UV(dU, dV) = P_\\Omega(dU S V^T + U S dV^T)

    This is pre-composed with projection of the pair (dU, dV) onto the tangent
    space of (U, V).
    """

    dU_t = dU - U @ (U.T @ dU)
    dV_t = dV - V @ (V.T @ dV)
    W = dU_t @ S @ V.T
    W += U @ S @ dV_t.T

    return W[rows, cols]


def jacobian_UV_adj(U, V, S, w, observed_mask):
    """Compute J*(W).

    The Jacobian adjoint is defined with respect to the inner product by

    <J(dU, dV, dS), W> = <(dU, dV, dS), J*(W)>

    Thus,

    J_UV*(W) = (P_\\Omega(W) V S^T, P_\\Omega(W)^T U S)

    This is projected back to the tangent space of (U, V).
    """

    W = np.zeros_like(observed_mask, dtype=U.dtype)
    W[observed_mask] = w

    dU = W @ V @ S.T
    dV = W.T @ U @ S

    dU -= U @ (U.T @ dU)
    dV -= V @ (V.T @ dV)

    return dU, dV


def pack(dU, dV):
    """Pack dU and dV to a single vector dx"""
    return np.concatenate([dU.ravel(), dV.ravel()])


def unpack(x, U_shape, V_shape):
    """Unpack the vector dx back to dU and dV."""
    nu = np.prod(U_shape)

    dU = x[:nu].reshape(U_shape)
    dV = x[nu:].reshape(V_shape)

    return dU, dV


def solve_gauss_newton_step(U, V, S, observed_mask, R, tol, damp):
    """Solve (J_UV* J_UV)dx = -J_UV* R.

    The Gauss-Newton step is the vector dx = (dU, dV), where dU and dV are
    tangent vectors in their respective Grassmann manifolds. The step is the
    least-squares solution of the system J_UV dx = -R, and it is computed using
    the LSMR algorithm.
    """

    nvars = U.size + V.size
    rows, cols = np.where(observed_mask)

    def matvec(x):
        dU, dV = unpack(x, U.shape, V.shape)
        return jacobian_UV(U, V, S, dU, dV, rows, cols)

    def rmatvec(y):
        dU, dV = jacobian_UV_adj(U, V, S, y, observed_mask)
        return pack(dU, dV)

    J_UV = LinearOperator(
        shape=(np.sum(observed_mask), nvars),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=U.dtype,
    )

    step = lsmr(J_UV, -R.ravel(), atol=tol, btol=tol, damp=damp)[0]

    return unpack(step, U.shape, V.shape)


def line_search(U, V, dU, dV, obj0, obj_fn, alpha0, tau, c):
    """Backtracking line search."""

    # Maximum depth of line search (no greater than 1e-16)
    max_ls = 16

    # Directional derivative of objective along search direction
    deriv = np.sum(dU**2) + np.sum(dV**2)

    # Initialize step size one update prior
    alpha = alpha0 / tau

    # Best pair (U, V)
    best = (U, V, obj0, alpha0)
    converged = False

    for _ in range(max_ls):
        # Retract step to Grassmann manifold
        U_try = retract_grassmann(U, alpha * dU)
        V_try = retract_grassmann(V, alpha * dV)

        # Recompute objective
        obj_try = obj_fn(U_try, V_try)

        # Update best pair (U, V)
        if obj_try < best[2]:
            best = (U_try, V_try, obj_try, alpha)

        # Armijo-Goldstein condition for sufficient decrease
        if obj_try < obj0 - c * alpha * deriv:
            converged = True
            break

        # Update step
        alpha *= tau

    # Check convergence
    if not converged:
        warn(
            "Sufficient decrease was not satisfied in line search even though the"
            f"depth limit was reached (1e-{max_ls}).",
            RuntimeWarning,
        )

    return best


def retract_grassmann(X, dX):
    """Retract the updated matrix X + dX back to the Grassmann manifold."""
    return np.linalg.qr(X + dX, mode="reduced")[0]


def optspace(X, dimensions=3, max_iter=10000, tol=1e-5, method="GD"):
    r"""Matrix completion using the OptSpace algorithm.

    OptSpace is an algorithm for recovering a low-rank matrix from a subset of observed
    entries. It uses gradient descent on the Grassmann manifold to find the optimal
    low-rank approximation.

    Parameters
    ----------
    X : array_like of shape (n_samples, n_features)
        A matrix with observed values and NaN for missing entries.
    dimensions : int, optional
        The rank of the matrix to recover. Default is 3.
    max_iter : int, optional
        Maximum number of iterations. Default is 10000.
    tol : float, optional
        Convergence tolerance. Default is 1e-5.
    method : {'GD', 'GN'}, optional
        The optimization method to use. Options are gradient descent ("GD", default)
        and Gauss-Newton ("GN").

    Returns
    -------
    ndarray of shape (n_samples, n_dimensions)
        The reconstructed optimal low-rank matrix.

    Raises
    ------
    ValueError
        If input is not 2D.
    ValueError
        If ``dimensions`` is not a positive integer less than or equal to
        ``min(n_samples, n_features)``.
    ValueError
        If ``method`` is not one of "GN" or "GD".
    ValueError
        If input is fully unobserved.

    See Also
    --------
    rpca

    Notes
    -----
    OptSpace was first described in [1]_.

    The algorithm proceeds as follows:

    1. Initialize `U`, `V` using trimmed SVD of the observed matrix.
    2. Iteratively:
       a. Compute optimal `S` given current `U`, `V`.
       b. Update `U`, `V` with the Gauss-Newton step `dU`, `dV`.
       c. Project `U`, `V` back to Grassmann manifold.

    References
    ----------
    .. [1] Keshavan, R. H., Montanari, A., & Oh, S. (2010). Matrix completion from a
       few entries. IEEE transactions on information theory, 56(6), 2980-2998.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import optspace
    >>> # Create a low-rank matrix
    >>> m, n, r = 600, 600, 5
    >>> rng = np.random.default_rng(0)
    >>> U_true = rng.normal(size=(m,r))
    >>> V_true = rng.normal(size=(n,r))
    >>> M_true = U_true @ V_true.T
    >>> # Mask some entries
    >>> M_obs = M_true.copy()
    >>> p_observe = 0.4  # 40% observed
    >>> mask = rng.random((m, n)) < p_observe
    >>> M_obs[~mask] = np.nan
    >>> # Recover the matrix
    >>> M_hat = optspace(M_obs, dimensions=r)
    """

    X = np.asarray(X, dtype=np.float64)

    # Validate input

    if X.ndim != 2:
        raise ValueError(f"Input must be 2D, got {X.ndim}D array.")

    elif type(dimensions) is not int:
        raise ValueError("Dimensions must be a positive integer")

    elif dimensions < 1 or dimensions > min(X.shape):
        raise ValueError(
            "Dimensions must be a positive integer less than or equal to "
            "min(n_samples, n_features)"
        )

    elif method not in ("GN", "GD"):
        raise ValueError("Method must be 'GN' or 'GD'")

    # Create observed mask
    observed_mask = ~np.isnan(X)

    # Check for unobserved rows or columns
    keep_rows, keep_cols = _check_unobserved(observed_mask)

    # Update matrix
    X_keep = X[np.ix_(keep_rows, keep_cols)]
    observed_mask = observed_mask[np.ix_(keep_rows, keep_cols)]

    n_observed = np.sum(observed_mask)
    m, n = X_keep.shape
    r = dimensions

    # Trim over-represented rows and columns
    X_trimmed = _trim(X_keep, observed_mask, m, n, n_observed)

    # Compute density for rescaling
    density = n_observed / (n * m)

    # Rescale observed values for sparse initialization
    X_trimmed /= density

    # Initialize with truncated SVD
    U, V = _svd_init(X_trimmed, r)

    # Vectorize data matrix over observed entries
    rows, cols = np.where(observed_mask)
    b = X_keep[rows, cols]

    # Iteratively solve for U, V, and S by minimizing the objective

    # Convergence parameters
    prev_obj = np.inf
    converged = False
    damp = 0
    alpha, tau, c = 1, 0.1, 1e-4

    # Objective function
    def _compute_obj(U, V):
        # Compute optimal S given current U, V
        S_curr = _solve_S(U, V, b, observed_mask, tol)

        # Compute current error
        R_curr = jacobian_S(U, V, S_curr, rows, cols) - b

        # Current objective (Frobenius norm of error over observed entries)
        obj_curr = np.sum(R_curr**2)

        return obj_curr

    for i in range(max_iter):
        # Compute optimal S given current U, V
        S = _solve_S(U, V, b, observed_mask, tol)

        # Compute current error
        R = jacobian_S(U, V, S, rows, cols) - b

        # Current objective (Frobenius norm of error over observed entries)
        obj = np.sum(R**2)

        # Check convergence
        if np.abs(prev_obj - obj) / obj < tol:
            # Gradient descent has a hierarchical convergence criterion
            # which incrementally sharpens the resolution of line search
            if method == "GD":
                if c < 1:
                    c *= 100
                else:
                    converged = True
                    break

            # Gauss-Newton exits immediately upon convergence
            if method == "GN":
                converged = True
                break

        prev_obj = obj

        # Update via gradient descent with line search
        if method == "GD":
            # Compute gradient directions
            dU, dV = jacobian_UV_adj(U, V, S, -R, observed_mask)

            # Perform backtracking line search
            U, V, obj, alpha = line_search(
                U, V, dU, dV, obj, _compute_obj, alpha, tau, c
            )

        # Update via Gauss-Newton
        if method == "GN":
            # Compute Gauss-Newton step
            dU, dV = solve_gauss_newton_step(U, V, S, observed_mask, R, tol, damp)

            # Retract updates back to Grassmann manifold
            U = retract_grassmann(U, dU)
            V = retract_grassmann(V, dV)

    # Check convergence
    if not converged:
        warn(
            f"OptSpace did not converge after 'max_iter' ({max_iter}) iterations.",
            RuntimeWarning,
        )

    # Form reconstructed matrix, and expand back to original shape
    X_hat = np.full(X.shape, np.nan)
    X_hat[np.ix_(keep_rows, keep_cols)] = U @ S @ V.T

    return X_hat
