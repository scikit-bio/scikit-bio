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

from collections.abc import Callable
from typing import NamedTuple
from warnings import warn

import numpy as np
from scipy.linalg import svd, qr, solve_triangular
from scipy.sparse.linalg import svds, lsmr, LinearOperator


# Default tuning values for the S-Jacobian basis, shared by the ``chunk_size``
# parameter of :func:`_gram_from_chunks` and the ``kappa_max`` parameter of
# :func:`_s_jacobian_basis`. Both were fixed from the empirical sweeps reported
# on the pull request that introduced the Gram path: peak memory of the blocked
# accumulation stays flat from 1e5 to 1e7 observed entries at this block size,
# and the Gram solve tracks the pivoted QR solve to machine precision while
# ``kappa(A)`` stays below 1e3, reaching a relative error of 1e-12 by
# ``kappa(A) = 1e4``.
_GRAM_CHUNK_SIZE = 200_000
_KAPPA_MAX = 1e3


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


def _kron_rows(X, Y):
    """Build the row-wise Kronecker products of two sets of factor rows.

    Returns the ``(n_rows, r ** 2)`` array whose ``k``-th row is
    ``kron(X[k], Y[k])``,

    out[k, a * r + b] = X[k, a] * Y[k, b],

    which is the row-major ordering matching the ``s.reshape(r, r)`` convention
    used throughout this module. This is the matrix form of ``J_S`` restricted
    to the rows passed in, so it builds either all of ``A``
    (:func:`_obs_basis`) or one block of its rows (:func:`_gram_from_chunks`).

    Parameters
    ----------
    X, Y : ndarray of shape (n_rows, r)
        Factor rows gathered at the observed entries.

    Returns
    -------
    ndarray of shape (n_rows, r ** 2)
        The row-wise Kronecker products.
    """

    r = X.shape[1]

    return (X[:, :, None] * Y[:, None, :]).reshape(-1, r**2)


def _gram_from_chunks(Ui, Vj, chunk_size=_GRAM_CHUNK_SIZE):
    """Accumulate the Gram matrix of the S-Jacobian in blocks.

    The matrix form of ``J_S`` is ``A[k] = kron(Ui[k], Vj[k])`` with one row per
    observed entry (see :func:`_obs_basis`). Its Gram matrix

    G = A^T A

    is only ``(r ** 2, r ** 2)``, so it can be accumulated one block of rows at
    a time,

    G = sum_blocks A_block^T A_block,

    materializing at most ``chunk_size`` rows of ``A`` at once. Peak memory is
    therefore ``O(chunk_size * r ** 2)`` and independent of ``n_observed``,
    whereas forming all of ``A`` costs ``O(n_observed * r ** 2)`` and reaches
    gigabytes on large inputs. The arithmetic cost stays
    ``O(n_observed * r ** 4)``, the same as forming ``A`` in one shot and
    multiplying.

    Parameters
    ----------
    Ui : ndarray of shape (n_observed, r)
        Rows of ``U`` gathered at the observed row indices.
    Vj : ndarray of shape (n_observed, r)
        Rows of ``V`` gathered at the observed column indices.
    chunk_size : int, optional
        Number of rows of ``A`` materialized at a time. Defaults to
        ``_GRAM_CHUNK_SIZE``. The result does not depend on it, only the peak
        memory of the accumulation does.

    Returns
    -------
    ndarray of shape (r ** 2, r ** 2)
        The Gram matrix ``A^T A``, symmetric positive semidefinite.
    """

    n_observed, r = Ui.shape

    G = np.zeros((r**2, r**2), dtype=np.result_type(Ui, Vj))

    for start in range(0, n_observed, chunk_size):
        stop = start + chunk_size

        A_block = _kron_rows(Ui[start:stop], Vj[start:stop])

        G += A_block.T @ A_block

    return G


def _numerical_rank(magnitudes, tol):
    """Count the magnitudes that are significant relative to the largest one.

    Both factorizations of the S-Jacobian decide rank the same way, from a
    set of magnitudes in any order: the pivoted QR thresholds the absolute
    values of the diagonal of ``R`` (:func:`_obs_basis`), and the Gram path
    thresholds the eigenvalues of ``A^T A`` (:func:`_gram_is_well_conditioned`).
    In both cases a value counts as nonzero only while it exceeds ``tol`` times
    the largest value, and a largest value that is itself zero (or, for
    eigenvalues of a semidefinite matrix computed in floating point, negative)
    leaves no scale to measure against and makes the rank zero.

    Parameters
    ----------
    magnitudes : ndarray of shape (n,)
        Non-negative magnitudes, in any order.
    tol : float
        Relative threshold below which a magnitude counts as zero.

    Returns
    -------
    int
        Number of magnitudes exceeding ``tol * magnitudes.max()``, or 0 if the
        sequence is empty or its largest value is not positive.
    """

    if magnitudes.size == 0:
        return 0

    largest = magnitudes.max()
    if largest <= 0.0:
        return 0

    return int(np.count_nonzero(magnitudes > tol * largest))


def _gram_is_well_conditioned(eigvals, tol=1e-12, kappa_max=_KAPPA_MAX):
    """Test whether solving through the Gram matrix is numerically safe.

    The normal equations square the conditioning of the problem,

    kappa(G) = kappa(A^T A) = kappa(A) ** 2,

    so the Gram solve is only usable while ``A`` itself is well conditioned.
    ``G`` is symmetric positive semidefinite and at most ``(r ** 2, r ** 2)``
    (9 by 9 at the default rank of 3), so its extreme eigenvalues, and hence
    ``kappa(A) = sqrt(kappa(G))``, are cheap to obtain exactly.

    A smallest eigenvalue that is non-positive or negligible relative to the
    largest means ``A`` does not have full column rank, in which case
    ``kappa(A)`` is infinite and the Gram solve is rejected without dividing
    by a near-zero value.

    The eigenvalues are passed in rather than ``G`` itself because
    :func:`_s_jacobian_basis` needs the full eigendecomposition anyway to build
    the inverse of ``G``, so decomposing once and sharing the result here
    avoids a second decomposition of the same matrix.

    Parameters
    ----------
    eigvals : ndarray of shape (r ** 2,)
        Eigenvalues of the Gram matrix of the S-Jacobian in ascending order,
        as returned by ``np.linalg.eigh`` or ``np.linalg.eigvalsh``.
    tol : float, optional
        Relative threshold below which the smallest eigenvalue counts as zero.
    kappa_max : float, optional
        Largest ``kappa(A)`` for which the Gram solve is accepted. Defaults to
        ``_KAPPA_MAX``.

    Returns
    -------
    bool
        Whether ``kappa(A)`` is finite and below ``kappa_max``.
    """

    # np.linalg.eigh and np.linalg.eigvalsh both order eigenvalues ascending.
    largest, smallest = eigvals[-1], eigvals[0]

    # Anything short of full rank leaves kappa(A) infinite, which also covers a
    # non-positive largest eigenvalue, so the division below is safe.
    if _numerical_rank(eigvals, tol) < eigvals.size:
        return False

    return np.sqrt(largest / smallest) < kappa_max


def _obs_basis(Ui, Vj, tol=1e-12):
    """Factorize the S-Jacobian restricted to the observed entries.

    ``J_S`` is the linear map ``dS -> P_\\Omega(U dS V^T)``. Written as a matrix
    acting on ``dS.ravel()`` it has one row per observed entry ``k = (i, j)``,

    A[k, a * r + b] = U[i, a] * V[j, b],

    i.e. ``A[k] = kron(U[i], V[j])``, built by :func:`_kron_rows`. ``A`` is
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

    This factorization is backward stable for any conditioning of ``A``, but it
    materializes ``A``. It is therefore used as the fallback of
    :func:`_s_jacobian_basis`, which prefers the cheaper Gram matrix whenever
    ``A`` is well enough conditioned for it.

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

    A = _kron_rows(Ui, Vj)

    Q, R, perm = qr(A, mode="economic", pivoting=True)

    # Pivoting orders |diag(R)| non-increasingly, so the numerical rank is the
    # length of the leading run of pivots above the relative threshold.
    rank = _numerical_rank(np.abs(np.diag(R)), tol)

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


def _s_jacobian_basis(Ui, Vj, tol=1e-12, kappa_max=_KAPPA_MAX):
    """Build the reusable S-Jacobian operators for the current (U, V).

    Both downstream operations, the least-squares solve for ``S`` and the
    projection onto the orthogonal complement of ``range(J_S)``, need a
    factorization of the S-Jacobian ``A[k] = kron(Ui[k], Vj[k])``. That
    factorization depends only on ``(U, V)``, so it is built once and reused by
    every solve and Jacobian product computed against the current ``(U, V)``
    within this outer iteration. Line search trials do not reuse it: each
    trial objective in ``_compute_obj`` rebuilds its own factorization for its
    own candidate ``(U, V)``.

    Two factorizations are available and this function selects between them:

    - The Gram matrix ``G = A^T A`` (:func:`_gram_from_chunks`), accumulated in
      blocks so that ``A`` is never fully materialized. Memory is independent of
      ``n_observed`` and it is several times faster than the QR, but solving
      through ``G`` squares the conditioning of ``A``.
    - The pivoted QR factorization of ``A`` (:func:`_obs_basis`), which is
      backward stable at any conditioning and handles rank deficiency, but
      materializes all ``O(n_observed * r ** 2)`` entries of ``A``.

    The Gram path is used while ``kappa(A) < kappa_max``, which covers the
    well-conditioned case the algorithm spends nearly all of its iterations in,
    and the QR is the fallback above that. Rank-deficient problems have an
    infinite ``kappa(A)`` and so take the fallback automatically.

    Parameters
    ----------
    Ui : ndarray of shape (n_observed, r)
        Rows of ``U`` gathered at the observed row indices.
    Vj : ndarray of shape (n_observed, r)
        Rows of ``V`` gathered at the observed column indices.
    tol : float, optional
        Relative threshold for rank detection. The same numeric value is
        passed to both paths, but each applies it to a different quantity: the
        QR path thresholds pivot magnitudes of ``A``, which scale like the
        singular values of ``A``, while the Gram path thresholds eigenvalues of
        ``G``, which are those singular values squared. The two are therefore
        not calibrated identically in terms of ``kappa(A)``: the same ``tol``
        flags rank deficiency around ``kappa(A) > 1 / tol`` for the QR and
        around ``kappa(A) > 1 / sqrt(tol)`` for the Gram. Neither bound governs
        the choice of path in practice, because ``kappa_max`` is far smaller
        than both and is reached first.
    kappa_max : float, optional
        Largest ``kappa(A)`` for which the Gram path is used. Defaults to
        ``_KAPPA_MAX``. This is what actually decides which path is taken.

    Returns
    -------
    solve_fn : callable
        ``solve_fn(b)`` returns the ``(r, r)`` least-squares solution ``S`` of
        ``J_S(S) = b`` over the observed entries (see :func:`_solve_S`).
    project_fn : callable
        ``project_fn(w)`` returns ``w`` projected onto the orthogonal
        complement of ``range(J_S)`` (see :func:`_project_S_complement`).
    """

    r = Ui.shape[1]

    G = _gram_from_chunks(Ui, Vj)

    # G is symmetric positive semidefinite, so one eigendecomposition serves
    # both the conditioning test and the inverse built from it below.
    eigvals, eigvecs = np.linalg.eigh(G)

    if _gram_is_well_conditioned(eigvals, tol, kappa_max):
        # Normal equations: S = G^-1 A^T b, where A^T b is J_S*(b), which is
        # already computed without forming A. Every eigenvalue is positive and
        # above tol * max(eigvals) here, since that is what the conditioning
        # test just checked, so the reciprocals are safe to take directly.
        G_inv = (eigvecs / eigvals) @ eigvecs.T

        def solve_fn(b):
            return (G_inv @ jacobian_S_adj(Ui, Vj, b)).reshape(r, r)

        def project_fn(w):
            # P_perp w = w - A (A^T A)^-1 A^T w, applied A-free.
            return w - jacobian_S(Ui, Vj, solve_fn(w))

    else:
        Q, R, perm = _obs_basis(Ui, Vj, tol)

        def solve_fn(b):
            return _solve_S(Q, R, perm, b, r)

        def project_fn(w):
            return _project_S_complement(Q, w)

    return solve_fn, project_fn


class _GNContext(NamedTuple):
    """State of the UV-component Jacobian at a fixed ``(U, V, S)``.

    :func:`jacobian_UV` and :func:`jacobian_UV_adj` are applied many times per
    outer iteration, once for every matrix-vector product LSMR requests, but
    everything they need besides the direction being mapped depends only on the
    current ``(U, V, S)``. Bundling that state keeps it built once per outer
    iteration, in particular the gathers ``US[rows]`` and ``VST[cols]``, which
    are ``(n_observed, r)`` and would otherwise be repeated on every product.

    Build one with :func:`_build_gn_context`. A context belongs to the exact
    ``(U, V, S)`` it was built from, so a line search trial, which proposes a
    different ``(U, V)``, needs its own.

    Attributes
    ----------
    rows, cols : ndarray of shape (n_observed,)
        Row and column indices of the observed entries.
    US_rows : ndarray of shape (n_observed, r)
        Rows of ``U @ S`` gathered at the observed row indices.
    VST_cols : ndarray of shape (n_observed, r)
        Rows of ``V @ S.T`` gathered at the observed column indices.
    project_fn : callable
        Projection onto the orthogonal complement of ``range(J_S)``, as
        returned by :func:`_s_jacobian_basis`.
    """

    rows: np.ndarray
    cols: np.ndarray
    US_rows: np.ndarray
    VST_cols: np.ndarray
    project_fn: Callable


def _build_gn_context(rows, cols, US, VST, project_fn):
    """Gather the S-scaled factor rows once and bundle the Jacobian state.

    Parameters
    ----------
    rows, cols : ndarray of shape (n_observed,)
        Row and column indices of the observed entries.
    US : ndarray of shape (m, r)
        The product ``U @ S``.
    VST : ndarray of shape (n, r)
        The product ``V @ S.T``.
    project_fn : callable
        Projection onto the orthogonal complement of ``range(J_S)`` built for
        the same ``(U, V)`` by :func:`_s_jacobian_basis`.

    Returns
    -------
    _GNContext
        The bundled state, with ``US[rows]`` and ``VST[cols]`` already gathered.
    """

    return _GNContext(rows, cols, US[rows], VST[cols], project_fn)


def jacobian_UV(U, V, dU, dV, ctx):
    """Compute J_UV(dU, dV).

    The Jacobian is

    J(dU, dV, dS) = J_S(dS) + J_U(dU) + J_V(dV),

    where J_UV(dU, dV) = J_U(dU) + J_V(dV). The UV component of the Jacobian
    determines how changes in U and V (dU, dV) contribute to the reconstruction
    error over the observed entries.

    J_UV(dU, dV) = P_\\Omega(dU S V^T + U S dV^T)

    Evaluated only at the observed entries ``k = (i, j)``, this is

    J_UV(dU, dV)[k] = dU[i] . (V S^T)[j] + (U S)[i] . dV[j],

    so the rows of ``U @ S`` and ``V @ S.T`` are gathered at the observed
    entries once per outer iteration, into ``ctx`` (see :func:`_GNContext`),
    and no dense ``(m, n)`` array is ever formed.

    This is pre-composed with projection of the pair (dU, dV) onto the tangent
    space of (U, V), and post-composed with ``ctx.project_fn``, the projection
    onto the orthogonal complement of the S tangent space built for the current
    ``(U, V)`` by :func:`_s_jacobian_basis`.
    """

    # Project input onto (U, V) tangent space
    dU_t = dU - U @ (U.T @ dU)
    dV_t = dV - V @ (V.T @ dV)

    # Compute Jacobian over the observed entries only
    w = (dU_t[ctx.rows] * ctx.VST_cols).sum(axis=1)
    w += (ctx.US_rows * dV_t[ctx.cols]).sum(axis=1)

    # Project output onto complement of S tangent space
    return ctx.project_fn(w)


def jacobian_UV_adj(U, V, w, ctx):
    """Compute J*(W).

    The Jacobian adjoint is defined with respect to the inner product by

    <J(dU, dV, dS), W> = <(dU, dV, dS), J*(W)>

    Thus,

    J_UV*(W) = (P_\\Omega(W) V S^T, P_\\Omega(W)^T U S)

    which over the observed entries ``k = (i, j)`` is the scatter-accumulation

    dU[i] += w[k] * (V S^T)[j],   dV[j] += w[k] * (U S)[i].

    The input is first projected onto the orthogonal complement of the S tangent
    space by ``ctx.project_fn`` and the output is projected back to the tangent
    space of (U, V). Both projections are self-adjoint, so this is the exact
    adjoint of :func:`jacobian_UV`.
    """

    m, r = U.shape
    n = V.shape[0]

    # Project input onto complement of S tangent space
    w = ctx.project_fn(w)

    # Compute Jacobian adjoint by scattering into the factor rows
    dU = np.empty((m, r), dtype=U.dtype)
    dV = np.empty((n, r), dtype=V.dtype)
    for q in range(r):
        dU[:, q] = np.bincount(ctx.rows, weights=w * ctx.VST_cols[:, q], minlength=m)
        dV[:, q] = np.bincount(ctx.cols, weights=w * ctx.US_rows[:, q], minlength=n)

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


def solve_gauss_newton_step(U, V, ctx, residual, tol, damp):
    """Solve (J_UV* J_UV)dx = -J_UV* residual.

    The Gauss-Newton step is the vector dx = (dU, dV), where dU and dV are
    tangent vectors in their respective Grassmann manifolds. The step is the
    least-squares solution of the system J_UV dx = -residual, and it is
    computed using the LSMR algorithm.

    ``ctx`` holds the Jacobian state at the current ``(U, V, S)``, which are
    fixed for the duration of this call, so it is built once by the caller (see
    :func:`_build_gn_context`) and closed over by the matrix-vector products
    rather than being rebuilt on every LSMR iteration.
    """

    nvars = U.size + V.size

    def matvec(x):
        dU, dV = unpack(x, U.shape, V.shape)
        return jacobian_UV(U, V, dU, dV, ctx)

    def rmatvec(y):
        dU, dV = jacobian_UV_adj(U, V, y, ctx)
        return pack(dU, dV)

    J_UV = LinearOperator(
        shape=(ctx.rows.size, nvars),
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
        # Gather the factor rows at the observed entries and factorize J_S.
        # Both are specific to the (U, V) being tried, so every trial of the
        # line search rebuilds them rather than reusing anything gathered for
        # the (U, V) of the current outer iteration.
        Ui_curr, Vj_curr = U[rows], V[cols]
        solve_curr, _ = _s_jacobian_basis(Ui_curr, Vj_curr)

        # Compute optimal S given current U, V
        S_curr = solve_curr(b)

        # Compute current error
        E_curr = jacobian_S(Ui_curr, Vj_curr, S_curr) - b

        # Current objective (Frobenius norm of error over observed entries)
        obj_curr = np.sum(E_curr**2)

        return obj_curr

    for _ in range(max_iter):
        # Gather the factor rows at the observed entries. The factorization of
        # J_S depends only on (U, V), so it is built once per iteration and
        # reused by every solve and Jacobian product.
        Ui, Vj = U[rows], V[cols]
        solve_fn, project_fn = _s_jacobian_basis(Ui, Vj)

        # Compute optimal S given current U, V
        S = solve_fn(b)

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

        # Bundle the state the Jacobian products below reuse, gathering the
        # S-scaled factor rows at the observed entries once for this iteration
        # instead of once per product. The bundle is tied to this (U, V, S), so
        # the line search, which tries other (U, V) pairs, must not be given it.
        ctx = _build_gn_context(rows, cols, U @ S, V @ S.T, project_fn)

        # Update via gradient descent with line search
        if method == "GD":
            # Compute gradient directions
            dU, dV = jacobian_UV_adj(U, V, -R, ctx)

            # Perform backtracking line search
            U, V, obj, alpha = line_search(
                U, V, dU, dV, obj, _compute_obj, alpha, tau, c
            )

        # Update via Gauss-Newton
        if method == "GN":
            # Compute Gauss-Newton step
            dU, dV = solve_gauss_newton_step(U, V, ctx, R, tol, damp)

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
