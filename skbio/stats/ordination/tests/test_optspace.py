# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import warnings
from unittest.mock import patch

import numpy as np

from skbio.stats.ordination._optspace import (
    optspace,
    _KAPPA_MAX,
    _build_gn_context,
    _gram_from_chunks,
    _gram_is_well_conditioned,
    _numerical_rank,
    _obs_basis,
    _project_S_complement,
    _s_jacobian_basis,
    _solve_S,
    jacobian_S,
    jacobian_S_adj,
    jacobian_UV,
    jacobian_UV_adj,
    solve_gauss_newton_step,
)


class TestOptSpace(unittest.TestCase):
    """Tests for OptSpace matrix completion algorithm."""

    def setUp(self):
        """Set up test fixtures."""
        rng = np.random.default_rng(0)
        self.tol = 1e-5

        # Create a low-rank matrix
        self.m, self.n, self.r = 100, 100, 3
        self.U_true = rng.normal(size=(self.m, self.r))
        self.V_true = rng.normal(size=(self.n, self.r))
        self.M_true = self.U_true @ self.V_true.T

        # Introduce noise
        noise_level = 0.01
        Z = noise_level * rng.normal(size=self.M_true.shape)
        M_noisy = self.M_true + Z

        # Mask some entries
        p = 0.4
        self.M_obs = M_noisy.copy()
        self.mask = rng.random((self.m, self.n)) < p
        self.M_obs[~self.mask] = np.nan

    def compute_rmse(self, M_hat):
        return np.sqrt(np.sum((M_hat - self.M_true)**2) / (self.m * self.n))

    def test_basic_completion(self):
        """Test basic matrix completion on low-rank matrix."""
        M_hat = optspace(self.M_obs, self.r, tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.0041953807728951865
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_gauss_newton(self):
        """Test low-rank matrix completion with Gauss-Newton"""
        M_hat = optspace(self.M_obs, self.r, method='GN', tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.004193496955751318
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_gradient_descent(self):
        """Test low-rank matrix completion with gradient descent"""
        M_hat = optspace(self.M_obs, self.r, method='GD', tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.0041953807728951865
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_positive_dimension(self):
        """Test error for nonpositive dimension."""
        dim = 0 # Not positive
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("positive", str(context.exception))

    def test_exceed_dimension(self):
        """Test error for dimension exceeding array size."""
        dim = max(self.m, self.n) + 1 # Too large
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("less than", str(context.exception))

    def test_integer_dimension(self):
        """Test error for noninteger dimension."""
        dim = 1.1 # Not an integer
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("integer", str(context.exception))

    def test_non_2d_error(self):
        """Test error on non-2D input."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(np.random.randn(5, 5, 5))
        self.assertIn("2D", str(context.exception))

    def test_invalid_method(self):
        """Test invalid method parameter."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=self.r, method=1)
        self.assertIn("Method", str(context.exception))

    def test_svds_exception(self):
        """Test sparse SVDS failure and fallback to dense SVD."""
        with patch("skbio.stats.ordination._optspace.svds",
            side_effect=RuntimeError("PROPACK failed")
        ):
            with self.assertWarns(RuntimeWarning) as context:
                _ = optspace(self.M_obs, self.r)

        self.assertIn("Sparse SVD failed", str(context.warning))

    def test_max_iter_exception(self):
        """Test manifold optimization warning when convergence is not achieved."""
        with self.assertWarns(RuntimeWarning) as context:
            _ = optspace(self.M_obs, self.r, max_iter=1)

        self.assertIn("did not converge", str(context.warning))

    def test_zero_max_iter_error(self):
        """Test error for max_iter of zero."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, self.r, max_iter=0)
        self.assertIn("positive", str(context.exception))

    def test_negative_max_iter_error(self):
        """Test error for negative max_iter."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, self.r, max_iter=-1)
        self.assertIn("positive", str(context.exception))

    def test_numpy_integer_max_iter_accepted(self):
        """Test that a NumPy integer type is accepted for max_iter."""
        with self.assertWarns(RuntimeWarning):
            _ = optspace(self.M_obs, self.r, max_iter=np.int64(1))

    def test_non_integer_max_iter_error(self):
        """Test error for non-integer max_iter."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, self.r, max_iter=2.5)
        self.assertIn("integer", str(context.exception))

    def test_max_ls_exception(self):
        """Test line search warning when convergence is not achieved."""
        with patch(
            "skbio.stats.ordination._optspace.retract_grassmann",
            side_effect=lambda X, dX: X
        ):
            with self.assertWarns(RuntimeWarning) as context:
                _ = optspace(self.M_obs, self.r, method='GD')

        self.assertIn("depth limit was reached", str(context.warning))

    def test_unobserved_row_exception(self):
        """Test warning when a row is not observed"""
        M_obs_row = self.M_obs
        nan_row = 1
        M_obs_row[nan_row, :] = np.nan

        with self.assertWarns(RuntimeWarning) as context:
            M_hat = optspace(M_obs_row, self.r)

        self.assertIn("remain as NaN", str(context.warning))
        self.assertTrue(np.all(np.isnan(M_hat[nan_row, :])))
    
    def test_unobserved_col_exception(self):
        """Test warning when a column is not observed"""
        M_obs_col = self.M_obs
        nan_col = 1
        M_obs_col[:, nan_col] = np.nan

        with self.assertWarns(RuntimeWarning) as context:
            M_hat = optspace(M_obs_col, self.r)

        self.assertIn("remain as NaN", str(context.warning))
        self.assertTrue(np.all(np.isnan(M_hat[:, nan_col])))

    def test_fully_unobserved_error(self):
        """Test error when input is fully unobserved."""
        M_obs = self.M_obs * np.nan
        with self.assertRaises(ValueError) as context:
            _ = optspace(M_obs, dimensions=self.r)
        self.assertIn("requires at least", str(context.exception))


class TestJacobian(unittest.TestCase):
    """Tests for the entry-indexed Jacobian operators used by OptSpace.

    The Jacobians are evaluated only at the observed entries. These tests check
    them against the equivalent dense formulation, which materializes the full
    (m, n) reconstruction and is therefore obviously correct but far slower.
    """

    def setUp(self):
        """Set up a small, well-conditioned factorization."""
        rng = np.random.default_rng(3)
        self.m, self.n, self.r = 23, 17, 3

        self.U = np.linalg.qr(rng.normal(size=(self.m, self.r)))[0]
        self.V = np.linalg.qr(rng.normal(size=(self.n, self.r)))[0]
        self.S = rng.normal(size=(self.r, self.r))

        self.mask = rng.random((self.m, self.n)) < 0.5
        self.rows, self.cols = np.where(self.mask)
        self.n_observed = self.rows.size

        self.Ui = self.U[self.rows]
        self.Vj = self.V[self.cols]
        self.US = self.U @ self.S
        self.VST = self.V @ self.S.T
        self.Q, self.R, self.perm = _obs_basis(self.Ui, self.Vj)
        self.solve_fn, self.project_fn = _s_jacobian_basis(self.Ui, self.Vj)
        self.ctx = _build_gn_context(
            self.rows, self.cols, self.US, self.VST, self.project_fn
        )

        self.rng = rng

    def dense_basis(self):
        """Dense (n_observed, r ** 2) matrix form of J_S."""
        return np.einsum("ka,kb->kab", self.Ui, self.Vj).reshape(
            self.n_observed, self.r**2
        )

    def dense_project(self, w):
        """Exact projection onto the complement of range(J_S), densely."""
        A = self.dense_basis()
        s = np.linalg.lstsq(A, w, rcond=None)[0]
        return w - A @ s

    def dense_jacobian_UV(self, dU, dV):
        """Dense reference for jacobian_UV (forms the full (m, n) product)."""
        dU_t = dU - self.U @ (self.U.T @ dU)
        dV_t = dV - self.V @ (self.V.T @ dV)
        W = dU_t @ self.S @ self.V.T + self.U @ self.S @ dV_t.T
        return self.dense_project(W[self.rows, self.cols])

    def dense_jacobian_UV_adj(self, w):
        """Dense reference for jacobian_UV_adj."""
        W = np.zeros((self.m, self.n))
        W[self.mask] = self.dense_project(w)
        dU = W @ self.V @ self.S.T
        dV = W.T @ self.U @ self.S
        dU -= self.U @ (self.U.T @ dU)
        dV -= self.V @ (self.V.T @ dV)
        return dU, dV

    def test_jacobian_S_matches_dense(self):
        """Entry-indexed J_S agrees with the dense reconstruction."""
        expected = (self.U @ self.S @ self.V.T)[self.rows, self.cols]
        observed = jacobian_S(self.Ui, self.Vj, self.S)
        np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_jacobian_S_adj_matches_dense(self):
        """Entry-indexed J_S* agrees with the dense scatter form."""
        w = self.rng.normal(size=self.n_observed)
        W = np.zeros((self.m, self.n))
        W[self.mask] = w
        expected = (self.U.T @ W @ self.V).ravel()
        observed = jacobian_S_adj(self.Ui, self.Vj, w)
        np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_obs_basis_reshape_convention(self):
        """The QR basis is built from the same s.reshape(r, r) ordering."""
        A = self.dense_basis()
        np.testing.assert_allclose(
            A @ self.S.ravel(), jacobian_S(self.Ui, self.Vj, self.S), atol=1e-12
        )

        # A[:, perm] == Q @ R for the retained columns
        np.testing.assert_allclose(A[:, self.perm], self.Q @ self.R, atol=1e-12)

    def test_solve_S_is_exact(self):
        """The triangular solve for S attains the least-squares minimum."""
        b = self.rng.normal(size=self.n_observed)
        A = self.dense_basis()

        S_solved = _solve_S(self.Q, self.R, self.perm, b, self.r)
        S_lstsq = np.linalg.lstsq(A, b, rcond=None)[0].reshape(self.r, self.r)

        np.testing.assert_allclose(S_solved, S_lstsq, atol=1e-10)

        # The residual must not exceed the reference least-squares residual.
        residual = np.linalg.norm(A @ S_solved.ravel() - b)
        reference = np.linalg.norm(A @ S_lstsq.ravel() - b)
        self.assertLessEqual(residual, reference * (1 + 1e-10))

    def test_project_S_complement(self):
        """P_perp is idempotent and annihilates the range of J_S."""
        w = self.rng.normal(size=self.n_observed)
        Pw = _project_S_complement(self.Q, w)

        # Idempotence
        np.testing.assert_allclose(
            _project_S_complement(self.Q, Pw), Pw, atol=1e-12
        )

        # Orthogonality: the projected vector has no J_S component left
        np.testing.assert_allclose(
            jacobian_S_adj(self.Ui, self.Vj, Pw),
            np.zeros(self.r**2),
            atol=1e-12,
        )

        # Agrees with the dense least-squares projection
        np.testing.assert_allclose(Pw, self.dense_project(w), atol=1e-9)

    def test_jacobian_UV_matches_dense(self):
        """Entry-indexed J_UV agrees with the dense formulation."""
        dU = self.rng.normal(size=(self.m, self.r))
        dV = self.rng.normal(size=(self.n, self.r))

        observed = jacobian_UV(self.U, self.V, dU, dV, self.ctx)
        np.testing.assert_allclose(
            observed, self.dense_jacobian_UV(dU, dV), atol=1e-9
        )

    def test_jacobian_UV_adj_matches_dense(self):
        """Entry-indexed J_UV* agrees with the dense formulation."""
        w = self.rng.normal(size=self.n_observed)

        dU, dV = jacobian_UV_adj(self.U, self.V, w, self.ctx)
        dU_dense, dV_dense = self.dense_jacobian_UV_adj(w)

        np.testing.assert_allclose(dU, dU_dense, atol=1e-9)
        np.testing.assert_allclose(dV, dV_dense, atol=1e-9)

    def test_jacobian_UV_adjoint_identity(self):
        """<J x, y> == <x, J* y> for random x and y.

        LSMR assumes matvec and rmatvec form an exact adjoint pair, so this is
        the defining property of the operator used by the Gauss-Newton step.
        """
        for seed in range(5):
            rng = np.random.default_rng(100 + seed)
            dU = rng.normal(size=(self.m, self.r))
            dV = rng.normal(size=(self.n, self.r))
            y = rng.normal(size=self.n_observed)

            Jx = jacobian_UV(self.U, self.V, dU, dV, self.ctx)
            JtY_U, JtY_V = jacobian_UV_adj(self.U, self.V, y, self.ctx)

            lhs = Jx @ y
            rhs = np.sum(dU * JtY_U) + np.sum(dV * JtY_V)
            self.assertAlmostEqual(lhs, rhs, delta=1e-12 * max(1.0, abs(lhs)))


class TestNumericalRank(unittest.TestCase):
    """Tests for the shared relative-threshold rank count.

    Both factorizations of the S-Jacobian decide rank from a non-increasing
    sequence of magnitudes, the pivoted QR from |diag(R)| and the Gram path
    from the eigenvalues of A^T A, so they share this helper.
    """

    def test_well_separated(self):
        """Magnitudes above and below the threshold are counted apart."""
        magnitudes = np.array([1.0, 1e-2, 1e-4, 1e-14, 1e-18])
        self.assertEqual(_numerical_rank(magnitudes, 1e-12), 3)

    def test_threshold_is_strict(self):
        """A magnitude exactly at the threshold does not count."""
        magnitudes = np.array([1.0, 1e-12])
        self.assertEqual(_numerical_rank(magnitudes, 1e-12), 1)

    def test_all_below_threshold(self):
        """Everything trailing the leading value can fall below it."""
        magnitudes = np.array([1.0, 1e-13, 1e-15])
        self.assertEqual(_numerical_rank(magnitudes, 1e-12), 1)

    def test_none_above_threshold_with_loose_tol(self):
        """A tolerance of one leaves nothing above the leading value."""
        magnitudes = np.array([1.0, 0.5, 0.25])
        self.assertEqual(_numerical_rank(magnitudes, 1.0), 0)

    def test_zero_leading_value(self):
        """An all-zero sequence has no scale to measure against."""
        self.assertEqual(_numerical_rank(np.zeros(5), 1e-12), 0)

    def test_negative_leading_value(self):
        """A non-positive leading value is rejected, not divided by.

        Eigenvalues of a semidefinite matrix can come back slightly negative in
        floating point, which the Gram path must treat as rank zero rather than
        thresholding against a negative scale.
        """
        magnitudes = np.array([-1e-18, -1.0])
        self.assertEqual(_numerical_rank(magnitudes, 1e-12), 0)

    def test_empty(self):
        """An empty sequence has rank zero."""
        self.assertEqual(_numerical_rank(np.array([]), 1e-12), 0)

    def test_returns_python_int(self):
        """The count is a plain int, usable directly as a slice bound."""
        rank = _numerical_rank(np.array([1.0, 0.5]), 1e-12)
        self.assertIsInstance(rank, int)

    def test_matches_obs_basis_rank(self):
        """The helper agrees with the rank the pivoted QR reports."""
        rng = np.random.default_rng(7)

        # A rank-r (not r ** 2) basis, from degenerate factor rows.
        r = 3
        Ui = np.tile(np.linalg.qr(rng.normal(size=(8, r)))[0][:1], (8, 1))
        Vj = np.linalg.qr(rng.normal(size=(8, r)))[0]

        Q, _, _ = _obs_basis(Ui, Vj)
        A = (Ui[:, :, None] * Vj[:, None, :]).reshape(-1, r**2)

        # Singular values come back non-increasing, like pivoted |diag(R)|.
        self.assertEqual(Q.shape[1], r)
        self.assertEqual(_numerical_rank(np.linalg.svd(A)[1], 1e-12), r)

    def test_matches_gram_conditioning_test(self):
        """Full rank of the ascending eigenvalues is what the Gram test needs."""
        # Well conditioned: full rank, so the conditioning test can proceed.
        eigvals = np.array([1.0, 2.0, 4.0])
        self.assertEqual(_numerical_rank(eigvals, 1e-12), eigvals.size)
        self.assertTrue(_gram_is_well_conditioned(eigvals))

        # Singular: rank deficient, so the Gram path is rejected.
        eigvals = np.array([0.0, 2.0, 4.0])
        self.assertLess(_numerical_rank(eigvals, 1e-12), eigvals.size)
        self.assertFalse(_gram_is_well_conditioned(eigvals))

    def test_order_independent(self):
        """The rank does not depend on the input order."""
        magnitudes = np.array([1e-4, 1.0, 1e-14, 1e-2, 1e-18])
        rank = _numerical_rank(magnitudes, 1e-12)
        self.assertEqual(rank, _numerical_rank(np.sort(magnitudes), 1e-12))
        self.assertEqual(rank, _numerical_rank(magnitudes[::-1], 1e-12))


class TestGNContext(unittest.TestCase):
    """Tests for the per-iteration bundle of UV-Jacobian state.

    The point of the bundle is that the (n_observed, r) gathers of ``U @ S``
    and ``V @ S.T`` happen once per outer iteration instead of once per LSMR
    matrix-vector product, and that each candidate ``(U, V)`` gets its own.
    """

    @staticmethod
    def problem(m=40, n=30, r=3, p=0.5, seed=5):
        """A partially observed low-rank matrix."""
        rng = np.random.default_rng(seed)
        M = rng.normal(size=(m, r)) @ rng.normal(size=(n, r)).T
        M[~(rng.random((m, n)) < p)] = np.nan
        return M

    def test_context_holds_the_gathered_rows(self):
        """The bundle carries the gathers, not the ungathered factors."""
        rng = np.random.default_rng(11)
        m, n, r = 9, 7, 3

        rows, cols = np.where(rng.random((m, n)) < 0.5)
        US = rng.normal(size=(m, r))
        VST = rng.normal(size=(n, r))

        def project_fn(w):
            return w

        ctx = _build_gn_context(rows, cols, US, VST, project_fn)

        np.testing.assert_array_equal(ctx.rows, rows)
        np.testing.assert_array_equal(ctx.cols, cols)
        np.testing.assert_array_equal(ctx.US_rows, US[rows])
        np.testing.assert_array_equal(ctx.VST_cols, VST[cols])
        self.assertIs(ctx.project_fn, project_fn)

    def test_gather_is_hoisted_out_of_the_lsmr_products(self):
        """Gauss-Newton builds one context per step, not one per product.

        Each Gauss-Newton step hands its context to LSMR, which asks for many
        matrix-vector products against it, so the context count must track the
        number of steps and stay far below the number of products.
        """
        M = self.problem()

        with patch(
            "skbio.stats.ordination._optspace._build_gn_context",
            wraps=_build_gn_context,
        ) as mock_ctx:
            with patch(
                "skbio.stats.ordination._optspace.solve_gauss_newton_step",
                wraps=solve_gauss_newton_step,
            ) as mock_step:
                with patch(
                    "skbio.stats.ordination._optspace.jacobian_UV",
                    wraps=jacobian_UV,
                ) as mock_matvec:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        optspace(M, dimensions=3, max_iter=5, method="GN")

        self.assertGreater(mock_ctx.call_count, 0)

        # One context per Gauss-Newton step.
        self.assertEqual(mock_ctx.call_count, mock_step.call_count)

        # Many matrix-vector products per context.
        self.assertGreater(mock_matvec.call_count, mock_ctx.call_count)

    def test_gradient_descent_builds_one_context_per_iteration(self):
        """Gradient descent needs one context per iteration and no more.

        It applies the adjoint once per iteration to get the gradient, so the
        context count must match that, and the line search must not pull a
        context of its own into being.
        """
        M = self.problem()

        with patch(
            "skbio.stats.ordination._optspace._build_gn_context",
            wraps=_build_gn_context,
        ) as mock_ctx:
            with patch(
                "skbio.stats.ordination._optspace.jacobian_UV_adj",
                wraps=jacobian_UV_adj,
            ) as mock_adj:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    optspace(M, dimensions=3, max_iter=5, method="GD")

        self.assertGreater(mock_ctx.call_count, 0)
        self.assertEqual(mock_ctx.call_count, mock_adj.call_count)

    def test_line_search_trials_refactorize_for_their_own_factors(self):
        """Every line search trial gets its own J_S factorization.

        The trial objective is evaluated at a retracted ``(U, V)`` that differs
        from the one the outer iteration gathered for, so reusing the outer
        iteration's factorization would evaluate the wrong objective. The
        trials therefore drive ``_s_jacobian_basis`` past the one call per
        iteration that builds the context.
        """
        M = self.problem()

        with patch(
            "skbio.stats.ordination._optspace._build_gn_context",
            wraps=_build_gn_context,
        ) as mock_ctx:
            with patch(
                "skbio.stats.ordination._optspace._s_jacobian_basis",
                wraps=_s_jacobian_basis,
            ) as mock_basis:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    optspace(M, dimensions=3, max_iter=5, method="GD")

        self.assertGreater(mock_basis.call_count, mock_ctx.call_count)


class TestGramBasis(unittest.TestCase):
    """Tests for the blocked Gram accumulation of the S-Jacobian."""

    @staticmethod
    def dense_A(Ui, Vj):
        """Dense (n_observed, r ** 2) matrix form of J_S."""
        r = Ui.shape[1]
        return (Ui[:, :, None] * Vj[:, None, :]).reshape(-1, r**2)

    @staticmethod
    def random_factors(n_observed, r, seed):
        """Orthonormal, well-conditioned factor rows."""
        rng = np.random.default_rng(seed)
        Ui = np.linalg.qr(rng.normal(size=(n_observed, r)))[0]
        Vj = np.linalg.qr(rng.normal(size=(n_observed, r)))[0]
        return Ui, Vj, rng

    def test_matches_direct_gram(self):
        """The blocked accumulation equals A^T A formed in one shot."""
        for n_observed, r, seed in [(200, 3, 0), (1000, 3, 1), (500, 5, 2)]:
            with self.subTest(n_observed=n_observed, r=r):
                Ui, Vj, _ = self.random_factors(n_observed, r, seed)
                A = self.dense_A(Ui, Vj)
                np.testing.assert_allclose(
                    _gram_from_chunks(Ui, Vj), A.T @ A, atol=1e-12
                )

    def test_matches_direct_gram_across_block_sizes(self):
        """Any block size, down to a single row, gives the same Gram matrix."""
        Ui, Vj, _ = self.random_factors(200, 3, 3)
        expected = self.dense_A(Ui, Vj).T @ self.dense_A(Ui, Vj)

        # 1 and 7 split the rows into many blocks, 1000 into a single block.
        for block in [1, 7, 137, 1000]:
            with self.subTest(block=block):
                observed = _gram_from_chunks(Ui, Vj, chunk_size=block)
                np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_gram_is_symmetric_psd(self):
        """The Gram matrix is symmetric with non-negative eigenvalues."""
        Ui, Vj, _ = self.random_factors(300, 3, 4)
        G = _gram_from_chunks(Ui, Vj)

        self.assertEqual(G.shape, (9, 9))
        np.testing.assert_allclose(G, G.T, atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(G).min(), -1e-12)


class TestSJacobianBasis(unittest.TestCase):
    """Tests for the condition-number-adaptive S-Jacobian dispatcher.

    The Gram path solves the normal equations, which squares the conditioning
    of the problem, so it is only selected while ``kappa(A)`` stays below
    ``_KAPPA_MAX``. These tests check that both paths produce the same
    least-squares solution and projection when they are both usable, and that
    the dispatcher switches to the pivoted QR above the threshold.
    """

    dense_A = staticmethod(TestGramBasis.dense_A)

    @staticmethod
    def factors_with_kappa(n_observed, r, kappa, seed):
        """Factor rows whose Jacobian has a condition number near ``kappa``.

        Scaling the columns of the orthonormal ``Ui`` built by
        :meth:`TestGramBasis.random_factors` by a geometric spread of singular
        values sets ``kappa(Ui)``, and hence ``kappa(A)``, directly.
        """
        Ui, Vj, rng = TestGramBasis.random_factors(n_observed, r, seed)
        return Ui * np.geomspace(kappa, 1.0, r), Vj, rng

    def assert_path(self, Ui, Vj, expected):
        """Run the dispatcher and assert which factorization it selected.

        ``_obs_basis`` is reached only by the QR fallback, so whether it was
        called identifies the path taken.
        """
        with patch(
            "skbio.stats.ordination._optspace._obs_basis", wraps=_obs_basis
        ) as mock_obs_basis:
            solve_fn, project_fn = _s_jacobian_basis(Ui, Vj)

        self.assertEqual("qr" if mock_obs_basis.called else "gram", expected)

        return solve_fn, project_fn

    def test_gram_path_matches_lstsq_and_qr(self):
        """On a well-conditioned basis, solve_fn attains the lstsq minimum."""
        Ui, Vj, rng = TestGramBasis.random_factors(400, 3, 5)
        A = self.dense_A(Ui, Vj)
        b = rng.normal(size=Ui.shape[0])

        solve_fn, _ = self.assert_path(Ui, Vj, "gram")

        S_gram = solve_fn(b)
        S_lstsq = np.linalg.lstsq(A, b, rcond=None)[0].reshape(3, 3)
        S_qr = _solve_S(*_obs_basis(Ui, Vj), b, 3)

        np.testing.assert_allclose(S_gram, S_lstsq, atol=1e-12)
        np.testing.assert_allclose(S_gram, S_qr, atol=1e-12)

        # The residual must not exceed the reference least-squares residual.
        residual = np.linalg.norm(A @ S_gram.ravel() - b)
        reference = np.linalg.norm(A @ S_lstsq.ravel() - b)
        self.assertLessEqual(residual, reference * (1 + 1e-10))

    def test_gram_path_projection_matches_dense_and_qr(self):
        """On a well-conditioned basis, project_fn is the exact projector."""
        Ui, Vj, rng = TestGramBasis.random_factors(400, 3, 6)
        A = self.dense_A(Ui, Vj)
        w = rng.normal(size=Ui.shape[0])

        _, project_fn = self.assert_path(Ui, Vj, "gram")

        Pw = project_fn(w)

        # Agrees with the dense least-squares projection and with the QR
        np.testing.assert_allclose(
            Pw, w - A @ np.linalg.lstsq(A, w, rcond=None)[0], atol=1e-12
        )
        np.testing.assert_allclose(
            Pw, _project_S_complement(_obs_basis(Ui, Vj)[0], w), atol=1e-12
        )

        # Idempotence, and no component of range(J_S) left
        np.testing.assert_allclose(project_fn(Pw), Pw, atol=1e-12)
        np.testing.assert_allclose(jacobian_S_adj(Ui, Vj, Pw), np.zeros(9), atol=1e-12)

    def test_well_conditioned_predicate_tracks_kappa(self):
        """The dispatch predicate flips at kappa(A) = _KAPPA_MAX."""
        for kappa, expected in [(1.0, True), (1e2, True), (1e5, False)]:
            with self.subTest(kappa=kappa):
                Ui, Vj, _ = self.factors_with_kappa(300, 3, kappa, 7)
                G = _gram_from_chunks(Ui, Vj)
                eigvals = np.linalg.eigvalsh(G)

                # The predicate reports kappa(A) = sqrt(kappa(G)) < _KAPPA_MAX
                self.assertEqual(_gram_is_well_conditioned(eigvals), expected)
                self.assertEqual(np.sqrt(np.linalg.cond(G)) < _KAPPA_MAX, expected)

    def test_dispatch_honors_kappa_max(self):
        """A tighter kappa_max sends an otherwise Gram-eligible basis to QR.

        ``kappa(A)`` is about 1e2 for this basis, below the default threshold
        of 1e3 but above the 1e1 passed here.
        """
        Ui, Vj, _ = self.factors_with_kappa(300, 3, 1e2, 8)

        with patch(
            "skbio.stats.ordination._optspace._obs_basis", wraps=_obs_basis
        ) as mock_obs_basis:
            _s_jacobian_basis(Ui, Vj, kappa_max=1e1)

        self.assertTrue(mock_obs_basis.called)

    def test_dispatch_below_threshold_uses_gram(self):
        """A basis with kappa(A) below the threshold takes the Gram path."""
        Ui, Vj, _ = self.factors_with_kappa(300, 3, 1e2, 8)
        self.assertLess(np.linalg.cond(self.dense_A(Ui, Vj)), _KAPPA_MAX)

        self.assert_path(Ui, Vj, "gram")

    def test_dispatch_above_threshold_uses_qr(self):
        """A basis with kappa(A) above the threshold takes the QR path."""
        Ui, Vj, rng = self.factors_with_kappa(300, 3, 1e6, 9)
        A = self.dense_A(Ui, Vj)
        self.assertGreater(np.linalg.cond(A), _KAPPA_MAX)

        solve_fn, _ = self.assert_path(Ui, Vj, "qr")

        # The QR stays accurate where the Gram solve would have lost roughly
        # kappa(A) ** 2 * eps of precision.
        b = rng.normal(size=Ui.shape[0])
        S_lstsq = np.linalg.lstsq(A, b, rcond=None)[0].reshape(3, 3)
        np.testing.assert_allclose(solve_fn(b), S_lstsq, rtol=1e-8)


class TestRankDeficient(unittest.TestCase):
    """Tests for the rank-deficient fallback in the S-Jacobian factorization."""

    @staticmethod
    def fewer_observations_than_unknowns(r=3, seed=1):
        """A 4 by 4 problem with 6 observed entries but r ** 2 unknowns in S.

        Returns the partially observed matrix along with the factor rows
        gathered at its observed entries, which form a rank-deficient ``A``.
        The observed values of ``M`` are ``M[~np.isnan(M)]``, in the same
        row-major order as the returned factor rows.
        """
        rng = np.random.default_rng(seed)

        M = np.full((4, 4), np.nan)
        for i, j in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 3)]:
            M[i, j] = rng.normal()

        rows, cols = np.where(~np.isnan(M))
        U = np.linalg.qr(rng.normal(size=(4, r)))[0]
        V = np.linalg.qr(rng.normal(size=(4, r)))[0]

        return M, U[rows], V[cols], rng

    @staticmethod
    def degenerate_factors(r=3, seed=1):
        """Factor rows of a fully observed 6 by 6 problem, with U degenerate.

        Repeating a single row of ``U`` across all of its rows collapses ``A``
        to rank ``r`` instead of ``r ** 2``.
        """
        rng = np.random.default_rng(seed)

        rows, cols = np.where(np.ones((6, 6), dtype=bool))
        U = np.tile(np.linalg.qr(rng.normal(size=(6, r)))[0][:1], (6, 1))
        V = np.linalg.qr(rng.normal(size=(6, r)))[0]

        return U[rows], V[cols], rng

    def test_fewer_observations_than_unknowns(self):
        """A is rank deficient when n_observed < r ** 2."""
        r = 3
        M, Ui, Vj, _ = self.fewer_observations_than_unknowns(r)

        Q, R, perm = _obs_basis(Ui, Vj)
        self.assertLess(Q.shape[1], r**2)
        self.assertEqual(Q.shape[1], R.shape[0])
        self.assertEqual(Q.shape[1], perm.size)

        S = _solve_S(Q, R, perm, M[~np.isnan(M)], r)
        self.assertTrue(np.all(np.isfinite(S)))

        # The dropped coefficients are exactly zero, not garbage.
        dropped = np.setdiff1d(np.arange(r**2), perm)
        np.testing.assert_array_equal(S.ravel()[dropped], 0.0)

    def test_degenerate_factors(self):
        """A is rank deficient when the factor rows are degenerate."""
        r = 3
        Ui, Vj, rng = self.degenerate_factors(r)

        Q, R, perm = _obs_basis(Ui, Vj)
        self.assertEqual(Q.shape[1], r)

        S = _solve_S(Q, R, perm, rng.normal(size=Ui.shape[0]), r)
        self.assertTrue(np.all(np.isfinite(S)))

    def test_completion_of_degenerate_input(self):
        """Matrix completion of a rank-deficient problem does not blow up."""
        M = self.fewer_observations_than_unknowns()[0]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            M_hat = optspace(M, dimensions=3, max_iter=50)

        self.assertTrue(np.all(np.isfinite(M_hat)))

    def test_rank_deficient_bases_dispatch_to_qr(self):
        """A rank-deficient A always routes to the QR fallback.

        A rank-deficient ``A`` has a singular Gram matrix, so ``kappa(A)`` is
        infinite and exceeds ``_KAPPA_MAX`` by construction. This is what makes
        the Gram path safe to prefer: it is never reached by a basis whose
        dropped columns the pivoted QR has to handle.
        """
        r = 3
        _, Ui_few, Vj_few, rng = self.fewer_observations_than_unknowns(r)
        Ui_deg, Vj_deg, _ = self.degenerate_factors(r)

        bases = [
            ("fewer observations than unknowns", Ui_few, Vj_few),
            ("degenerate factors", Ui_deg, Vj_deg),
        ]

        for name, Ui, Vj in bases:
            with self.subTest(basis=name):
                G = _gram_from_chunks(Ui, Vj)
                self.assertFalse(_gram_is_well_conditioned(np.linalg.eigvalsh(G)))

                with patch(
                    "skbio.stats.ordination._optspace._obs_basis", wraps=_obs_basis
                ) as mock_obs_basis:
                    solve_fn, _ = _s_jacobian_basis(Ui, Vj)
                self.assertTrue(mock_obs_basis.called)

                # The fallback keeps the basic solution of the QR path.
                Q, R, perm = _obs_basis(Ui, Vj)
                b = rng.normal(size=Ui.shape[0])
                np.testing.assert_allclose(
                    solve_fn(b), _solve_S(Q, R, perm, b, r), atol=1e-12
                )

    def test_degenerate_completion_dispatches_to_qr(self):
        """Every iteration of a degenerate completion takes the QR path."""
        M = self.fewer_observations_than_unknowns()[0]

        # The dispatcher builds the Gram matrix on every call but reaches
        # _obs_basis only when it falls back, so equal call counts mean every
        # iteration took the QR path.
        with patch(
            "skbio.stats.ordination._optspace._gram_from_chunks",
            wraps=_gram_from_chunks,
        ) as mock_gram:
            with patch(
                "skbio.stats.ordination._optspace._obs_basis", wraps=_obs_basis
            ) as mock_obs_basis:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    optspace(M, dimensions=3, max_iter=50)

        self.assertGreater(mock_gram.call_count, 0)
        self.assertEqual(mock_obs_basis.call_count, mock_gram.call_count)


class TestSolverAgreement(unittest.TestCase):
    """Gradient descent and Gauss-Newton share the exact S solve."""

    def test_gd_and_gn_agree(self):
        """Both methods converge to the same completed matrix."""
        rng = np.random.default_rng(0)
        m, n, r = 100, 100, 3
        M_true = rng.normal(size=(m, r)) @ rng.normal(size=(n, r)).T
        M_obs = M_true + 0.01 * rng.normal(size=M_true.shape)
        M_obs[~(rng.random((m, n)) < 0.4)] = np.nan

        M_gd = optspace(M_obs, r, method="GD", tol=1e-5)
        M_gn = optspace(M_obs, r, method="GN", tol=1e-5)

        rel_diff = np.linalg.norm(M_gd - M_gn) / np.linalg.norm(M_gn)
        self.assertLess(rel_diff, 1e-3)


# TODO: More varied sizes, shapes, and missingness among matrices
# Vary total size (m*n), ratio (m/n), portion p of observed entries,
# and different singular value structures.


if __name__ == '__main__':
    unittest.main()