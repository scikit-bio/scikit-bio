# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""OptSpace matrix completion algorithm.

This module provides the OptSpace algorithm for low-rank matrix completion from
partially observed entries. It is used by the Robust PCA (RPCA) ordination method.

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
from scipy.linalg import svd, qr, solve_triangular
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


def _obs_basis(Ui, Vj, tol=1e-12):
    """Factorize the S-Jacobian restricted to the observed entries.

    ``J_S`` is the linear map ``dS -> P_\\Omega(U dS V^T)``. Written as a matrix
    acting on ``dS.ravel()`` it has one row per observed entry ``k = (i, j)``,

    A[k, a * r + b] = U[i, a] * V[j, b],

    i.e. ``A[k] = kron(U[i], V[j])``, which matches the row-major
    ``s.reshape(r, r)`` convention used throughout this module. ``A`` is
    ``(n_observed, r ** 2)``, so it is thin (only 9 columns at the default rank
    of 3) and can be factorized once per outer iteration and reused by every
    matrix-vector product.

    The economy QR factorization ``A[:, perm] = Q R`` yields both operations
    needed downstream: ``Q`` is an orthonormal basis of ``range(J_S)``, and
    ``R`` back-substitutes the least-squares solve for ``S``.

    Column pivoting is used so that rank deficiency (possible when the observed
    entries do not determine all ``r ** 2`` degrees of freedom of ``S``) is
    detected from the magnitudes of the diagonal of ``R``, which pivoting orders
    non-increasingly. Trailing columns whose pivot falls below ``tol`` relative
    to the leading pivot are dropped, giving the basic (minimum-column-support)
    least-squares solution rather than a division by a near-zero pivot.

    Parameters
    ----------
    Ui : ndarray of shape (n_observed, r)
        Rows of ``U`` gathered at the observed row indices.
    Vj : ndarray of shape (n_observed, r)
        Rows of ``V`` gathered at the observed column indices.
    tol : float, optional
        Relative pivot threshold for rank detection.

    Returns
    -------
    Q : ndarray of shape (n_observed, rank)
        Orthonormal basis of ``range(J_S)``.
    R : ndarray of shape (rank, rank)
        Upper triangular factor of the retained columns.
    perm : ndarray of shape (rank,)
        Indices of the retained columns of ``A``.
    """

    r = Ui.shape[1]

    # A[k, a * r + b] = Ui[k, a] * Vj[k, b]
    A = (Ui[:, :, None] * Vj[:, None, :]).reshape(-1, r**2)

    Q, R, perm = qr(A, mode="economic", pivoting=True)

    # Pivoting orders |diag(R)| non-increasingly, so the numerical rank is the
    # length of the leading run of pivots above the relative threshold.
    pivots = np.abs(np.diag(R))
    if pivots.size == 0 or pivots[0] == 0.0:
        rank = 0
    else:
        rank = int(np.count_nonzero(pivots > tol * pivots[0]))

    return Q[:, :rank], R[:rank, :rank], perm[:rank]


def _project_S_complement(Q, w):
    """Project w onto the orthogonal complement of the S tangent space.

    ``Q`` is an orthonormal basis of ``range(J_S)`` as returned by
    :func:`_obs_basis`, so the projector onto the complement is

    P_\\perp w = w - Q (Q^T w).

    This is exact and costs ``O(n_observed * r ** 2)``.
    """

    return w - Q @ (Q.T @ w)


def jacobian_S(Ui, Vj, S):
    """Compute J_S(dS).

    The Jacobian is

    J(dU, dV, dS) = J_S(dS) + J_U(dU) + J_V(dV).

    The S component of the Jacobian determines how changes in S (dS) contribute to
    the reconstruction error over the observed entries.

    J_S(dS) = P_\\Omega(U dS V^T).

    Evaluated only at the observed entries, this is

    J_S(dS)[k] = U[i] dS V[j]^T,

    which avoids forming the dense ``(m, n)`` product.
    """

    return (Ui @ S * Vj).sum(axis=1)


def jacobian_S_adj(Ui, Vj, w):
    """Compute J_S*(W).

    The Jacobian adjoint is defined with respect to the inner product by

    <J(dU, dV, dS), W> = <(dU, dV, dS), J*(W)>

    Thus,

    J_S*(W) = U^T P_\\Omega(W) V = sum_k w[k] outer(U[i], V[j]),

    where the sum runs over the observed entries ``k = (i, j)``.
    """

    return np.einsum("k,ka,kb->ab", w, Ui, Vj).ravel()


def _solve_S(Q, R, perm, b, r):
    """Compute optimal S given U and V.

    Solves the least squares problem to find the optimal S that
    minimizes the reconstruction difference on the observed entries:

    arg min_S ||P_\\Omega(U V S^T - M_observed)||_F^2

    where P_\\Omega is the projection onto the observed entries.
    This is the least-squares solution to the system

    J_S(dS) = -R

    Given the pivoted QR factorization ``A[:, perm] = Q R`` of the S-Jacobian
    (see :func:`_obs_basis`), the solution follows from a single triangular
    back-substitution, ``R s[perm] = Q^T b``, with the dropped (rank-deficient)
    coefficients left at zero. This is exact, unlike the iterative solve it
    replaces.
    """

    s = np.zeros(r**2, dtype=b.dtype)

    if perm.size:
        s[perm] = solve_triangular(R, Q.T @ b)

    return s.reshape(r, r)


def jacobian_UV(U, V, dU, dV, rows, cols, US, VST, Q):
    """Compute J_UV(dU, dV).

    The Jacobian is

    J(dU, dV, dS) = J_S(dS) + J_U(dU) + J_V(dV),

    where J_UV(dU, dV) = J_U(dU) + J_V(dV). The UV component of the Jacobian
    determines how changes in U and V (dU, dV) contribute to the reconstruction
    error over the observed entries.

    J_UV(dU, dV) = P_\\Omega(dU S V^T + U S dV^T)

    Evaluated only at the observed entries ``k = (i, j)``, this is

    J_UV(dU, dV)[k] = dU[i] . (V S^T)[j] + (U S)[i] . dV[j],

    so ``US = U @ S`` and ``VST = V @ S.T`` are precomputed once per outer
    iteration and no dense ``(m, n)`` array is ever formed.

    This is pre-composed with projection of the pair (dU, dV) onto the tangent
    space of (U, V), and post-composed with projection onto the orthogonal
    complement of the S tangent space.
    """

    # Project input onto (U, V) tangent space
    dU_t = dU - U @ (U.T @ dU)
    dV_t = dV - V @ (V.T @ dV)

    # Compute Jacobian over the observed entries only
    w = (dU_t[rows] * VST[cols]).sum(axis=1)
    w += (US[rows] * dV_t[cols]).sum(axis=1)

    # Project output onto complement of S tangent space
    return _project_S_complement(Q, w)


def jacobian_UV_adj(U, V, w, rows, cols, US, VST, Q):
    """Compute J*(W).

    The Jacobian adjoint is defined with respect to the inner product by

    <J(dU, dV, dS), W> = <(dU, dV, dS), J*(W)>

    Thus,

    J_UV*(W) = (P_\\Omega(W) V S^T, P_\\Omega(W)^T U S)

    which over the observed entries ``k = (i, j)`` is the scatter-accumulation

    dU[i] += w[k] * (V S^T)[j],   dV[j] += w[k] * (U S)[i].

    The input is first projected onto the orthogonal complement of the S tangent
    space and the output is projected back to the tangent space of (U, V). Both
    projections are self-adjoint, so this is the exact adjoint of
    :func:`jacobian_UV`.
    """

    m, r = U.shape
    n = V.shape[0]

    # Project input onto complement of S tangent space
    w = _project_S_complement(Q, w)

    # Compute Jacobian adjoint by scattering into the factor rows
    dU = np.empty((m, r), dtype=U.dtype)
    dV = np.empty((n, r), dtype=V.dtype)
    for q in range(r):
        dU[:, q] = np.bincount(rows, weights=w * VST[cols, q], minlength=m)
        dV[:, q] = np.bincount(cols, weights=w * US[rows, q], minlength=n)

    # Project output onto (U, V) tangent space
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


def solve_gauss_newton_step(U, V, rows, cols, US, VST, Q_S, residual, tol, damp):
    """Solve (J_UV* J_UV)dx = -J_UV* residual.

    The Gauss-Newton step is the vector dx = (dU, dV), where dU and dV are
    tangent vectors in their respective Grassmann manifolds. The step is the
    least-squares solution of the system J_UV dx = -residual, and it is
    computed using the LSMR algorithm.

    ``US``, ``VST`` and ``Q_S`` depend only on the current ``(U, V, S)``, which
    are fixed for the duration of this call, so they are computed once by the
    caller and closed over by the matrix-vector products rather than being
    rebuilt on every LSMR iteration.
    """

    nvars = U.size + V.size

    def matvec(x):
        dU, dV = unpack(x, U.shape, V.shape)
        return jacobian_UV(U, V, dU, dV, rows, cols, US, VST, Q_S)

    def rmatvec(y):
        dU, dV = jacobian_UV_adj(U, V, y, rows, cols, US, VST, Q_S)
        return pack(dU, dV)

    J_UV = LinearOperator(
        shape=(rows.size, nvars),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=U.dtype,
    )

    step = lsmr(J_UV, -residual.ravel(), atol=tol, btol=tol, damp=damp)[0]

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
    entries. It uses optimization on the Grassmann manifold to find the optimal
    low-rank approximation.

    Parameters
    ----------
    X : array_like of shape (n_samples, n_features)
        A matrix with observed values and NaN for missing entries.
    dimensions : int, optional
        The rank of the matrix to recover. Default is 3.
    max_iter : int, optional
        Maximum number of iterations. Default is 10000.

        .. versionchanged:: 0.7.4
            A non-positive or non-integer value now raises ``ValueError``
            instead of failing with an unrelated error partway through.
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
        If ``max_iter`` is not a positive integer.
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
       b. Update `U`, `V` with the manifold optimization step `dU`, `dV`.
       c. Project `U`, `V` back to Grassmann manifold.

    References
    ----------
    .. [1] Keshavan, R. H., Montanari, A., & Oh, S. (2010). Matrix completion from a
       few entries. IEEE transactions on information theory, 56(6), 2980-2998.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.ordination import optspace

    Create a low-rank matrix.

    >>> m, n, r = 100, 100, 5
    >>> rng = np.random.default_rng(42)
    >>> U_true = rng.normal(size=(m,r))
    >>> V_true = rng.normal(size=(n,r))
    >>> M_true = U_true @ V_true.T

    Mask some entries.

    >>> M_obs = M_true.copy()
    >>> p_obs = 0.4  # 40% observed
    >>> mask = rng.random((m, n)) < p_obs
    >>> M_obs[~mask] = np.nan

    Recover the matrix.

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

    elif not np.issubdtype(type(max_iter), np.integer) or max_iter < 1:
        raise ValueError("Max_iter must be a positive integer")

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
        # Gather the factor rows at the observed entries and factorize J_S
        Ui_curr, Vj_curr = U[rows], V[cols]
        Q_curr, R_curr, perm_curr = _obs_basis(Ui_curr, Vj_curr)

        # Compute optimal S given current U, V
        S_curr = _solve_S(Q_curr, R_curr, perm_curr, b, r)

        # Compute current error
        E_curr = jacobian_S(Ui_curr, Vj_curr, S_curr) - b

        # Current objective (Frobenius norm of error over observed entries)
        obj_curr = np.sum(E_curr**2)

        return obj_curr

    for _ in range(max_iter):
        # Gather the factor rows at the observed entries. The basis Q of
        # range(J_S) and its triangular factor depend only on (U, V), so they
        # are built once per iteration and reused by every Jacobian product.
        Ui, Vj = U[rows], V[cols]
        Q_S, R_S, perm_S = _obs_basis(Ui, Vj)

        # Compute optimal S given current U, V
        S = _solve_S(Q_S, R_S, perm_S, b, r)

        # Precompute the S-scaled factors used by the Jacobian products
        US = U @ S
        VST = V @ S.T

        # Compute current error
        R = jacobian_S(Ui, Vj, S) - b

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
            dU, dV = jacobian_UV_adj(U, V, -R, rows, cols, US, VST, Q_S)

            # Perform backtracking line search
            U, V, obj, alpha = line_search(
                U, V, dU, dV, obj, _compute_obj, alpha, tau, c
            )

        # Update via Gauss-Newton
        if method == "GN":
            # Compute Gauss-Newton step
            dU, dV = solve_gauss_newton_step(
                U, V, rows, cols, US, VST, Q_S, R, tol, damp
            )

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
