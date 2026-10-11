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
    _gram_from_chunks,
    _gram_is_well_conditioned,
    _obs_basis,
    _project_S_complement,
    _s_jacobian_basis,
    _solve_S,
    jacobian_S,
    jacobian_S_adj,
    jacobian_UV,
    jacobian_UV_adj,
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

        observed = jacobian_UV(
            self.U,
            self.V,
            dU,
            dV,
            self.rows,
            self.cols,
            self.US,
            self.VST,
            self.project_fn,
        )
        np.testing.assert_allclose(
            observed, self.dense_jacobian_UV(dU, dV), atol=1e-9
        )

    def test_jacobian_UV_adj_matches_dense(self):
        """Entry-indexed J_UV* agrees with the dense formulation."""
        w = self.rng.normal(size=self.n_observed)

        dU, dV = jacobian_UV_adj(
            self.U,
            self.V,
            w,
            self.rows,
            self.cols,
            self.US,
            self.VST,
            self.project_fn,
        )
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

            Jx = jacobian_UV(
                self.U,
                self.V,
                dU,
                dV,
                self.rows,
                self.cols,
                self.US,
                self.VST,
                self.project_fn,
            )
            JtY_U, JtY_V = jacobian_UV_adj(
                self.U,
                self.V,
                y,
                self.rows,
                self.cols,
                self.US,
                self.VST,
                self.project_fn,
            )

            lhs = Jx @ y
            rhs = np.sum(dU * JtY_U) + np.sum(dV * JtY_V)
            self.assertAlmostEqual(lhs, rhs, delta=1e-12 * max(1.0, abs(lhs)))


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
                with patch("skbio.stats.ordination._optspace._GRAM_CHUNK_SIZE", block):
                    observed = _gram_from_chunks(Ui, Vj)
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

        Scaling the columns of an orthonormal factor by a geometric spread of
        singular values sets ``kappa(Ui)``, and hence ``kappa(A)``, directly.
        """
        rng = np.random.default_rng(seed)
        Ui = np.linalg.qr(rng.normal(size=(n_observed, r)))[0]
        Ui = Ui * np.geomspace(kappa, 1.0, r)
        Vj = np.linalg.qr(rng.normal(size=(n_observed, r)))[0]
        return Ui, Vj, rng

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

                # The predicate reports kappa(A) = sqrt(kappa(G)) < _KAPPA_MAX
                self.assertEqual(_gram_is_well_conditioned(G), expected)
                self.assertEqual(np.sqrt(np.linalg.cond(G)) < _KAPPA_MAX, expected)

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

    def test_fewer_observations_than_unknowns(self):
        """A is rank deficient when n_observed < r ** 2."""
        rng = np.random.default_rng(1)
        r = 3

        # 6 observed entries but 9 unknowns in S
        M = np.full((4, 4), np.nan)
        for i, j in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 3)]:
            M[i, j] = rng.normal()

        rows, cols = np.where(~np.isnan(M))
        U = np.linalg.qr(rng.normal(size=(4, r)))[0]
        V = np.linalg.qr(rng.normal(size=(4, r)))[0]

        Q, R, perm = _obs_basis(U[rows], V[cols])
        self.assertLess(Q.shape[1], r**2)
        self.assertEqual(Q.shape[1], R.shape[0])
        self.assertEqual(Q.shape[1], perm.size)

        S = _solve_S(Q, R, perm, M[rows, cols], r)
        self.assertTrue(np.all(np.isfinite(S)))

        # The dropped coefficients are exactly zero, not garbage.
        dropped = np.setdiff1d(np.arange(r**2), perm)
        np.testing.assert_array_equal(S.ravel()[dropped], 0.0)

    def test_degenerate_factors(self):
        """A is rank deficient when the factor rows are degenerate."""
        rng = np.random.default_rng(1)
        r = 3

        rows, cols = np.where(np.ones((6, 6), dtype=bool))
        U = np.tile(np.linalg.qr(rng.normal(size=(6, r)))[0][:1], (6, 1))
        V = np.linalg.qr(rng.normal(size=(6, r)))[0]

        Q, R, perm = _obs_basis(U[rows], V[cols])
        self.assertEqual(Q.shape[1], r)

        S = _solve_S(Q, R, perm, rng.normal(size=rows.size), r)
        self.assertTrue(np.all(np.isfinite(S)))

    def test_completion_of_degenerate_input(self):
        """Matrix completion of a rank-deficient problem does not blow up."""
        rng = np.random.default_rng(1)

        M = np.full((4, 4), np.nan)
        for i, j in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 3)]:
            M[i, j] = rng.normal()

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
        rng = np.random.default_rng(1)
        r = 3

        # Fewer observed entries than unknowns in S
        M = np.full((4, 4), np.nan)
        for i, j in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 3)]:
            M[i, j] = rng.normal()
        rows, cols = np.where(~np.isnan(M))
        U = np.linalg.qr(rng.normal(size=(4, r)))[0]
        V = np.linalg.qr(rng.normal(size=(4, r)))[0]
        bases = [("fewer observations than unknowns", U[rows], V[cols])]

        # Degenerate factor rows
        rows, cols = np.where(np.ones((6, 6), dtype=bool))
        U = np.tile(np.linalg.qr(rng.normal(size=(6, r)))[0][:1], (6, 1))
        V = np.linalg.qr(rng.normal(size=(6, r)))[0]
        bases.append(("degenerate factors", U[rows], V[cols]))

        for name, Ui, Vj in bases:
            with self.subTest(basis=name):
                G = _gram_from_chunks(Ui, Vj)
                self.assertFalse(_gram_is_well_conditioned(G))

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
        rng = np.random.default_rng(1)

        M = np.full((4, 4), np.nan)
        for i, j in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (2, 3)]:
            M[i, j] = rng.normal()

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