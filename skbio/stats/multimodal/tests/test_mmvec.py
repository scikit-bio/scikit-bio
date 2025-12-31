# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Tests for MMvec implementation."""

import io
import sys
import unittest

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from skbio.stats.composition import clr_inv as softmax
from skbio.stats.multimodal import mmvec, MMvecResults
from skbio.stats.multimodal.tests._simulation import random_multimodal


class TestMMvecRecovery(unittest.TestCase):
    """Test that MMvec recovers true embeddings from synthetic data."""

    def setUp(self):
        """Build small simulation matching original mmvec test."""
        np.random.seed(1)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=8,
            num_samples=150,
            latent_dim=2,
            sigmaQ=2,
            microbe_total=1000,
            metabolite_total=10000,
            seed=1,
        )
        (
            self.microbes,
            self.metabolites,
            self.X,
            self.B,
            self.U,
            self.Ubias,
            self.V,
            self.Vbias,
        ) = res

        num_train = 10
        self.trainX = self.microbes.iloc[:-num_train]
        self.testX = self.microbes.iloc[-num_train:]
        self.trainY = self.metabolites.iloc[:-num_train]
        self.testY = self.metabolites.iloc[-num_train:]

    def test_recovers_embeddings(self):
        """Verify model recovers true embeddings from synthetic data."""
        np.random.seed(1)

        result = mmvec(
            self.trainX,
            self.trainY,
            test_microbes=self.testX,
            test_metabolites=self.testY,
            n_components=2,
            optimizer="adam",
            max_iter=1000,
            batch_size=50,
            learning_rate=1e-3,
            beta_1=0.8,
            beta_2=0.9,
            random_state=0,
        )

        # Check result type
        self.assertIsInstance(result, MMvecResults)

        # Verify U embeddings are recovered
        # Compare pairwise distances using Spearman correlation
        u_r, u_p = spearmanr(
            pdist(result.microbe_embeddings.values[:, :-1]),  # exclude bias
            pdist(self.U),
        )
        self.assertGreater(u_r, 0.5, f"U correlation too low: {u_r}")
        self.assertLess(u_p, 0.05, f"U p-value too high: {u_p}")

        # Verify V embeddings are recovered
        # Compare pairwise distances between metabolites
        v_r, v_p = spearmanr(
            pdist(result.metabolite_embeddings.values[1:, :-1]),  # exclude ref & bias
            pdist(self.V.T),  # V is (latent_dim, n_metabolites-1), so V.T is (n_met-1, dim)
        )
        self.assertGreater(v_r, 0.5, f"V correlation too low: {v_r}")
        self.assertLess(v_p, 0.05, f"V p-value too high: {v_p}")

        # Verify conditional probabilities match
        d1 = self.U.shape[0]

        # Compute expected probabilities from true parameters
        U_ = np.hstack((np.ones((d1, 1)), self.Ubias, self.U))
        V_ = np.vstack(
            (self.Vbias, np.ones((1, self.V.shape[1])), self.V)
        )
        exp = softmax(np.hstack((np.zeros((d1, 1)), U_ @ V_)))

        # Compute result probabilities
        res = softmax(result.ranks.values)

        s_r, s_p = spearmanr(np.ravel(res), np.ravel(exp))
        self.assertGreater(s_r, 0.5, f"Probability correlation too low: {s_r}")
        self.assertLess(s_p, 0.05, f"Probability p-value too high: {s_p}")

    def test_cv_error_reasonable(self):
        """Cross-validation error should be reasonable."""
        np.random.seed(1)

        result = mmvec(
            self.trainX,
            self.trainY,
            test_microbes=self.testX,
            test_metabolites=self.testY,
            n_components=2,
            optimizer="adam",
            max_iter=1000,
            batch_size=50,
            learning_rate=1e-3,
            beta_1=0.8,
            beta_2=0.9,
            random_state=0,
        )

        # Check that CV error is reasonable (< 500 as in original test)
        final_cv = result.convergence["cv_error"].dropna().iloc[-1]
        self.assertLess(final_cv, 500, f"CV error too high: {final_cv}")


class TestMMvecGradients(unittest.TestCase):
    """Test that analytical gradients match numerical gradients."""

    def test_multinomial_gradient(self):
        """Verify multinomial softmax gradient is correct."""
        from skbio.stats.multimodal._mmvec import _multinomial_loglik_and_grad

        # Small test case
        logits = np.array([[0.0, 1.0, 2.0], [0.0, -1.0, 1.0]])
        y = np.array([[5.0, 10.0, 15.0], [20.0, 5.0, 5.0]])

        # Get analytical gradient
        _, grad_analytical = _multinomial_loglik_and_grad(logits, y)

        # Compute numerical gradient
        eps = 1e-5
        grad_numerical = np.zeros_like(logits)
        for i in range(logits.shape[0]):
            for j in range(logits.shape[1]):
                logits_plus = logits.copy()
                logits_plus[i, j] += eps
                logits_minus = logits.copy()
                logits_minus[i, j] -= eps

                loss_plus, _ = _multinomial_loglik_and_grad(logits_plus, y)
                loss_minus, _ = _multinomial_loglik_and_grad(logits_minus, y)

                grad_numerical[i, j] = (loss_plus - loss_minus) / (2 * eps)

        np.testing.assert_allclose(
            grad_analytical, grad_numerical, rtol=1e-4, atol=1e-6
        )

    def test_full_model_gradients(self):
        """Verify all model gradients against numerical differentiation."""
        from skbio.stats.multimodal._mmvec import _MMvecModel

        np.random.seed(42)
        n_microbes = 5
        n_metabolites = 8
        n_components = 2
        n_samples = 20

        # Create model
        rng = np.random.default_rng(42)
        model = _MMvecModel(
            n_microbes=n_microbes,
            n_metabolites=n_metabolites,
            n_components=n_components,
            rng=rng,
        )

        # Generate small test data
        X = np.random.randint(0, 100, size=(n_samples, n_microbes)).astype(
            np.float64
        )
        Y = np.random.randint(0, 100, size=(n_samples, n_metabolites)).astype(
            np.float64
        )

        # Compute analytical gradients with fixed seed for reproducibility
        batch_seed = 123
        _, grads = model.loss_and_gradients(
            X, Y, batch_size=10, rng=np.random.default_rng(batch_seed)
        )

        # Verify each parameter numerically
        # Use same seed for each call to get consistent batch
        eps = 1e-5
        for param_name in ["U", "b_U", "V", "b_V"]:
            param = getattr(model, param_name)
            numerical_grad = np.zeros_like(param)

            for idx in np.ndindex(param.shape):
                original = param[idx]
                param[idx] = original + eps
                loss_plus, _ = model.loss_and_gradients(
                    X, Y, batch_size=10, rng=np.random.default_rng(batch_seed)
                )
                param[idx] = original - eps
                loss_minus, _ = model.loss_and_gradients(
                    X, Y, batch_size=10, rng=np.random.default_rng(batch_seed)
                )
                param[idx] = original  # restore

                numerical_grad[idx] = (loss_plus - loss_minus) / (2 * eps)

            # Check relative error
            rel_error = np.abs(grads[param_name] - numerical_grad) / (
                np.abs(numerical_grad) + 1e-8
            )
            max_error = rel_error.max()
            self.assertLess(
                max_error,
                1e-3,
                f"{param_name}: max relative error {max_error}",
            )


class TestMMvecBasic(unittest.TestCase):
    """Basic functionality tests for MMvec."""

    def test_output_shapes(self):
        """Check that output shapes are correct."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=10,
            num_metabolites=15,
            num_samples=50,
            latent_dim=3,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result = mmvec(
            microbes,
            metabolites,
            n_components=3,
            max_iter=10,
            random_state=42,
        )

        # Check shapes
        n_microbes = microbes.shape[1]
        n_metabolites = metabolites.shape[1]
        n_components = 3

        self.assertEqual(
            result.microbe_embeddings.shape,
            (n_microbes, n_components + 1),  # +1 for bias
        )
        self.assertEqual(
            result.metabolite_embeddings.shape,
            (n_metabolites, n_components + 1),  # +1 for bias
        )
        self.assertEqual(result.ranks.shape, (n_microbes, n_metabolites))

    def test_ranks_row_centered(self):
        """Check that ranks are row-centered."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            random_state=42,
        )

        # Ranks should be approximately row-centered
        row_means = result.ranks.values.mean(axis=1)
        np.testing.assert_allclose(row_means, 0, atol=1e-10)

    def test_reproducibility(self):
        """Check that results are reproducible with same random state."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result1 = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            random_state=123,
        )

        result2 = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            random_state=123,
        )

        np.testing.assert_allclose(
            result1.ranks.values, result2.ranks.values, rtol=1e-10
        )

    def test_batch_normalization_modes(self):
        """Test that both batch normalization modes work."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        # Legacy mode (batch_normalization only applies to Adam)
        result_legacy = mmvec(
            microbes,
            metabolites,
            n_components=2,
            optimizer="adam",
            max_iter=50,
            batch_normalization="legacy",
            random_state=42,
        )
        self.assertIsInstance(result_legacy, MMvecResults)

        # Unbiased mode
        result_unbiased = mmvec(
            microbes,
            metabolites,
            n_components=2,
            optimizer="adam",
            max_iter=50,
            batch_normalization="unbiased",
            random_state=42,
        )
        self.assertIsInstance(result_unbiased, MMvecResults)

        # Results should differ due to different normalization
        self.assertFalse(
            np.allclose(result_legacy.ranks.values, result_unbiased.ranks.values)
        )


class TestMMvecResults(unittest.TestCase):
    """Test MMvecResults class."""

    def test_probabilities(self):
        """Test probability computation from ranks."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            random_state=42,
        )

        probs = result.probabilities()

        # Probabilities should sum to 1 per row
        np.testing.assert_allclose(probs.sum(axis=1), 1.0, rtol=1e-6)

        # Probabilities should be positive
        self.assertTrue((probs.values >= 0).all())


class TestMMvecValidation(unittest.TestCase):
    """Test input validation for MMvec."""

    def setUp(self):
        """Create valid test data."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]

    def test_mismatched_samples_raises(self):
        """Mismatched sample counts should raise ValueError."""
        microbes = self.microbes.iloc[:40]  # 40 samples
        metabolites = self.metabolites.iloc[:50]  # 50 samples

        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, metabolites, max_iter=1)

        self.assertIn("same number of samples", str(ctx.exception))

    def test_zero_microbe_column_raises(self):
        """All-zero microbe column should raise ValueError."""
        microbes = self.microbes.copy()
        microbes.iloc[:, 0] = 0  # Set first column to all zeros

        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, self.metabolites, max_iter=1)

        self.assertIn("all-zero columns", str(ctx.exception))
        self.assertIn("OTU_0", str(ctx.exception))

    def test_zero_metabolite_column_raises(self):
        """All-zero metabolite column should raise ValueError."""
        metabolites = self.metabolites.copy()
        metabolites.iloc[:, 0] = 0  # Set first column to all zeros

        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, metabolites, max_iter=1)

        self.assertIn("all-zero columns", str(ctx.exception))
        self.assertIn("metabolite_0", str(ctx.exception))

    def test_zero_sample_row_raises(self):
        """All-zero sample row should raise ValueError."""
        microbes = self.microbes.copy()
        microbes.iloc[0, :] = 0  # Set first row to all zeros

        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, self.metabolites, max_iter=1)

        self.assertIn("all-zero counts", str(ctx.exception))

    def test_zero_u_prior_scale_raises(self):
        """Zero u_prior_scale should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            mmvec(
                self.microbes,
                self.metabolites,
                u_prior_scale=0,
                max_iter=1,
            )

        self.assertIn("u_prior_scale must be positive", str(ctx.exception))

    def test_negative_u_prior_scale_raises(self):
        """Negative u_prior_scale should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            mmvec(
                self.microbes,
                self.metabolites,
                u_prior_scale=-1.0,
                max_iter=1,
            )

        self.assertIn("u_prior_scale must be positive", str(ctx.exception))

    def test_zero_v_prior_scale_raises(self):
        """Zero v_prior_scale should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            mmvec(
                self.microbes,
                self.metabolites,
                v_prior_scale=0,
                max_iter=1,
            )

        self.assertIn("v_prior_scale must be positive", str(ctx.exception))

    def test_negative_v_prior_scale_raises(self):
        """Negative v_prior_scale should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            mmvec(
                self.microbes,
                self.metabolites,
                v_prior_scale=-0.5,
                max_iter=1,
            )

        self.assertIn("v_prior_scale must be positive", str(ctx.exception))


class TestMMvecInputTypes(unittest.TestCase):
    """Test different input types for MMvec."""

    def test_sparse_pandas_input(self):
        """Test with sparse pandas DataFrame."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        microbes_dense, metabolites_dense = res[0], res[1]

        # Convert to sparse pandas DataFrame
        microbes_sparse = microbes_dense.astype(pd.SparseDtype("float64", 0))
        metabolites_sparse = metabolites_dense.astype(
            pd.SparseDtype("float64", 0)
        )

        result = mmvec(
            microbes_sparse,
            metabolites_sparse,
            n_components=2,
            max_iter=10,
            random_state=42,
        )

        # Should produce valid results
        self.assertEqual(result.ranks.shape, (8, 10))
        np.testing.assert_allclose(
            result.ranks.values.mean(axis=1), 0, atol=1e-10
        )


class TestMMvecParameterBehavior(unittest.TestCase):
    """Test parameter behavior for MMvec."""

    def setUp(self):
        """Create test data."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]

    def test_nonzero_prior_mean_affects_results(self):
        """Non-zero prior means should affect results."""
        result_default = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            u_prior_mean=0.0,
            v_prior_mean=0.0,
            random_state=42,
        )

        result_nonzero = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            u_prior_mean=5.0,
            v_prior_mean=5.0,
            random_state=42,
        )

        # Results should differ
        self.assertFalse(
            np.allclose(
                result_default.microbe_embeddings.values,
                result_nonzero.microbe_embeddings.values,
            )
        )

    def test_verbose_output(self):
        """Verbose mode should produce output."""
        # Capture stdout
        captured = io.StringIO()
        sys.stdout = captured

        try:
            mmvec(
                self.microbes,
                self.metabolites,
                n_components=2,
                optimizer="adam",
                max_iter=5,
                verbose=True,
                random_state=42,
            )
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        self.assertIn("Epoch", output)
        self.assertIn("Loss", output)


class TestMMvecOutputVerification(unittest.TestCase):
    """Test output structure and properties."""

    def setUp(self):
        """Create test data with specific IDs."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]

    def test_index_columns_preserved(self):
        """Index and column names should be preserved from input."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=10,
            random_state=42,
        )

        # Microbe IDs preserved
        self.assertEqual(
            list(result.microbe_embeddings.index),
            list(self.microbes.columns),
        )
        self.assertEqual(
            list(result.ranks.index),
            list(self.microbes.columns),
        )

        # Metabolite IDs preserved
        self.assertEqual(
            list(result.metabolite_embeddings.index),
            list(self.metabolites.columns),
        )
        self.assertEqual(
            list(result.ranks.columns),
            list(self.metabolites.columns),
        )

    def test_convergence_structure(self):
        """Convergence DataFrame should have correct structure."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=10,
            random_state=42,
        )

        # Check columns exist
        self.assertIn("iteration", result.convergence.columns)
        self.assertIn("loss", result.convergence.columns)
        self.assertIn("cv_error", result.convergence.columns)

        # Iterations should be sequential positive integers
        iterations = result.convergence["iteration"].values
        self.assertEqual(iterations[0], 1)
        self.assertTrue(np.all(np.diff(iterations) == 1))

    def test_loss_decreases(self):
        """Loss should generally decrease during training."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="adam",
            max_iter=100,
            random_state=42,
        )

        losses = result.convergence["loss"].values

        # Compare first 10% average to last 10% average
        n = len(losses)
        early_loss = losses[: n // 10].mean()
        late_loss = losses[-n // 10 :].mean()

        self.assertLess(late_loss, early_loss)


class TestMMvecLBFGS(unittest.TestCase):
    """Test L-BFGS optimizer for MMvec."""

    def setUp(self):
        """Create test data."""
        np.random.seed(42)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=10,
            num_samples=50,
            latent_dim=2,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]
        self.U = res[4]
        self.V = res[6]

    def test_lbfgs_produces_valid_results(self):
        """L-BFGS optimizer should produce valid results."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=100,
            random_state=42,
        )

        # Check result type
        self.assertIsInstance(result, MMvecResults)

        # Check shapes
        self.assertEqual(result.microbe_embeddings.shape, (8, 3))
        self.assertEqual(result.metabolite_embeddings.shape, (10, 3))
        self.assertEqual(result.ranks.shape, (8, 10))

        # Ranks should be row-centered
        row_means = result.ranks.values.mean(axis=1)
        np.testing.assert_allclose(row_means, 0, atol=1e-10)

    def test_lbfgs_reproducibility(self):
        """L-BFGS should be deterministic with same seed."""
        result1 = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            random_state=123,
        )

        result2 = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            random_state=123,
        )

        np.testing.assert_allclose(
            result1.ranks.values, result2.ranks.values, rtol=1e-10
        )

    def test_lbfgs_recovers_structure(self):
        """L-BFGS should recover embedding structure from synthetic data."""
        # Use more data for better recovery
        np.random.seed(1)
        res = random_multimodal(
            num_microbes=8,
            num_metabolites=8,
            num_samples=150,
            latent_dim=2,
            sigmaQ=2,
            microbe_total=1000,
            metabolite_total=10000,
            seed=1,
        )
        microbes, metabolites = res[0], res[1]
        U_true = res[4]

        result = mmvec(
            microbes,
            metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=500,
            random_state=42,
        )

        # Check correlation of pairwise distances
        u_r, u_p = spearmanr(
            pdist(result.microbe_embeddings.values[:, :-1]),
            pdist(U_true),
        )
        self.assertGreater(u_r, 0.3, f"U correlation too low: {u_r}")

    def test_lbfgs_with_test_data(self):
        """L-BFGS should work with test data for CV error."""
        # Split data
        train_microbes = self.microbes.iloc[:40]
        test_microbes = self.microbes.iloc[40:]
        train_metabolites = self.metabolites.iloc[:40]
        test_metabolites = self.metabolites.iloc[40:]

        result = mmvec(
            train_microbes,
            train_metabolites,
            test_microbes=test_microbes,
            test_metabolites=test_metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=100,
            random_state=42,
        )

        # CV error should be tracked
        cv_errors = result.convergence["cv_error"].dropna()
        self.assertGreater(len(cv_errors), 0)

    def test_invalid_optimizer_raises(self):
        """Invalid optimizer should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            mmvec(
                self.microbes,
                self.metabolites,
                optimizer="invalid",
                max_iter=1,
            )

        self.assertIn("optimizer must be", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
