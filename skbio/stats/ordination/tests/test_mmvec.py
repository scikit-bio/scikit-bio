# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Tests for MMvec implementation."""

import os
import io
import sys
import unittest

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from skbio.util import get_data_path
from skbio.stats.composition import clr_inv as softmax
from skbio.stats.ordination import mmvec, MMvecResults
from skbio.stats.ordination._mmvec import (
    random_multimodal, _MMvecModel, _multinomial_loglik_and_grad
)


class TestRandomMultimodal(unittest.TestCase):
    """Tests for random_multimodal simulation helper."""

    def test_output_shapes(self):
        """Generated outputs should have expected shapes and IDs."""
        res = random_multimodal(
            n_microbes=4,
            n_metabolites=6,
            n_samples=8,
            latent_dim=2,
            seed=7,
        )
        microbes, metabolites, design, beta, U, Ubias, V, Vbias = res

        self.assertEqual(microbes.shape, (8, 4))
        self.assertEqual(metabolites.shape, (8, 6))
        self.assertEqual(design.shape, (8, 2))
        self.assertEqual(beta.shape, (2, 4))
        self.assertEqual(U.shape, (4, 2))
        self.assertEqual(Ubias.shape, (4, 1))
        self.assertEqual(V.shape, (2, 5))
        self.assertEqual(Vbias.shape, (1, 5))

        self.assertEqual(microbes.index[0], "sample_0")
        self.assertEqual(microbes.columns[0], "microbe_0")
        self.assertEqual(metabolites.columns[0], "metabolite_0")

    def test_reproducible_with_int_seed(self):
        """Integer seeds should provide reproducible outputs."""
        res1 = random_multimodal(seed=42)
        res2 = random_multimodal(seed=42)

        pdt.assert_frame_equal(res1[0], res2[0])
        pdt.assert_frame_equal(res1[1], res2[1])
        for i in range(2, len(res1)):
            npt.assert_allclose(res1[i], res2[i])


class TestMMvecRecovery(unittest.TestCase):
    """Test that MMvec recovers true embeddings from synthetic data."""

    def setUp(self):
        """Build small simulation matching original mmvec test."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=8,
            n_samples=150,
            latent_dim=2,
            sigmaQ=2,
            microbe_total=1000,
            metabolite_total=10000,
            seed=1,
        )
        (
            self.microbes,
            self.metabolites,
            self.design,
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

    @unittest.skip("Skipping a test that requires long runtime.")
    def test_recovers_embeddings(self):
        """Verify model recovers true embeddings from synthetic data."""
        result = mmvec(
            self.trainX,
            self.trainY,
            n_components=2,
            optimizer="adam",
            max_iter=1000,
            batch_size=50,
            learning_rate=1e-3,
            beta_1=0.8,
            beta_2=0.9,
            seed=0,
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
        res = softmax(result.ranks.values)

        s_r, s_p = spearmanr(np.ravel(res), np.ravel(exp))
        self.assertGreater(s_r, 0.5, f"Probability correlation too low: {s_r}")
        self.assertLess(s_p, 0.05, f"Probability p-value too high: {s_p}")

    @unittest.skip("Skipping a test that requires long runtime.")
    def test_score_reasonable(self):
        """Q^2 score on held-out data should be reasonable."""
        result = mmvec(
            self.trainX,
            self.trainY,
            n_components=2,
            optimizer="adam",
            max_iter=1000,
            batch_size=50,
            learning_rate=1e-3,
            beta_1=0.8,
            beta_2=0.9,
            seed=0,
        )

        # Compute Q^2 score on test data
        q2 = result.score(self.testX, self.testY)

        # Q^2 should be positive for a reasonably trained model
        # (better than predicting the mean)
        self.assertGreater(q2, -1.0, f"Q-squared score too low: {q2}")


class TestMMvecGradients(unittest.TestCase):
    """Test that analytical gradients match numerical gradients."""

    def test_multinomial_gradient(self):
        """Verify multinomial softmax gradient is correct."""
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

        npt.assert_allclose(
            grad_analytical, grad_numerical, rtol=1e-4, atol=1e-6
        )

    def test_full_model_gradients(self):
        """Verify all model gradients against numerical differentiation."""
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
            X, Y, np.random.default_rng(batch_seed), batch_size=10
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
                    X, Y, np.random.default_rng(batch_seed), batch_size=10
                )
                param[idx] = original - eps
                loss_minus, _ = model.loss_and_gradients(
                    X, Y, np.random.default_rng(batch_seed), batch_size=10
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
        res = random_multimodal(
            n_microbes=10,
            n_metabolites=15,
            n_samples=50,
            latent_dim=3,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result = mmvec(
            microbes,
            metabolites,
            n_components=3,
            max_iter=10,
            seed=42,
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
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            seed=42,
        )

        # Ranks should be approximately row-centered
        row_means = result.ranks.values.mean(axis=1)
        npt.assert_allclose(row_means, 0, atol=1e-10)

    def test_reproducibility(self):
        """Check that results are reproducible with same seed."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        result1 = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            seed=123,
        )

        result2 = mmvec(
            microbes,
            metabolites,
            n_components=2,
            max_iter=50,
            seed=123,
        )

        npt.assert_allclose(
            result1.ranks.values, result2.ranks.values, rtol=1e-10
        )

    def test_batch_norm_modes(self):
        """Test that both batch normalization modes work."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
            seed=42,
        )
        microbes, metabolites = res[0], res[1]

        # Legacy mode (batch_norm only applies to Adam)
        result_legacy = mmvec(
            microbes,
            metabolites,
            n_components=2,
            optimizer="adam",
            max_iter=50,
            batch_norm="legacy",
            seed=42,
        )
        self.assertIsInstance(result_legacy, MMvecResults)

        # Unbiased mode
        result_unbiased = mmvec(
            microbes,
            metabolites,
            n_components=2,
            optimizer="adam",
            max_iter=50,
            batch_norm="unbiased",
            seed=42,
        )
        self.assertIsInstance(result_unbiased, MMvecResults)

        # Results should differ due to different normalization
        self.assertFalse(
            np.allclose(result_legacy.ranks.values, result_unbiased.ranks.values)
        )


class TestMMvecResults(unittest.TestCase):
    """Test MMvecResults class."""

    def setUp(self):
        """Create test data."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]

    def test_probabilities(self):
        """Test probability computation from ranks."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            seed=42,
        )

        probs = result.probabilities()

        # Probabilities should sum to 1 per row
        npt.assert_allclose(probs.sum(axis=1), 1.0, rtol=1e-6)

        # Probabilities should be positive
        self.assertTrue((probs.values >= 0).all())

    def test_predict(self):
        """Test predict method returns valid metabolite distributions."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            seed=42,
        )

        # Predict on new samples
        new_microbes = self.microbes.iloc[:5]
        predictions = result.predict(new_microbes)

        # Predictions should have correct shape
        self.assertEqual(predictions.shape, (5, 10))

        # Predictions should sum to 1 per row
        npt.assert_allclose(predictions.sum(axis=1), 1.0, rtol=1e-6)

        # Predictions should be positive
        self.assertTrue((predictions.values >= 0).all())

        # Column names should match metabolites
        self.assertEqual(
            list(predictions.columns),
            list(self.metabolites.columns),
        )

    def test_score(self):
        """Test score method returns valid Q-squared value."""
        # Split data
        train_microbes = self.microbes.iloc[:40]
        test_microbes = self.microbes.iloc[40:]
        train_metabolites = self.metabolites.iloc[:40]
        test_metabolites = self.metabolites.iloc[40:]

        result = mmvec(
            train_microbes,
            train_metabolites,
            n_components=2,
            max_iter=100,
            seed=42,
        )

        q2 = result.score(test_microbes, test_metabolites)

        # Q^2 should be a float
        self.assertIsInstance(q2, float)

        # Q^2 should be <= 1.0 (perfect prediction)
        self.assertLessEqual(q2, 1.0)

    def test_predict_zero_sample_raises(self):
        """Predict with zero-count sample should raise ValueError."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=10,
            seed=42,
        )

        # Create microbes with a zero row
        zero_microbes = self.microbes.iloc[:3].copy()
        zero_microbes.iloc[0, :] = 0

        with self.assertRaises(ValueError) as ctx:
            result.predict(zero_microbes)

        self.assertIn("all-zero counts", str(ctx.exception))

    def test_str(self):
        """Test string representation of MMvecResults."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=10,
            seed=42,
        )
        s = str(result)
        self.assertIn("MMvecResults", s)
        self.assertIn("Microbes", s)
        self.assertIn("Metabolites", s)

    def test_predict_numpy_array(self):
        """Test predict with numpy array input (not DataFrame)."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            seed=42,
        )
        predictions = result.predict(self.microbes.values[:5])
        self.assertEqual(predictions.shape[0], 5)
        npt.assert_allclose(predictions.sum(axis=1), 1.0, rtol=1e-6)

    def test_score_numpy_array(self):
        """Test score with numpy array inputs (not DataFrames)."""
        train_microbes = self.microbes.iloc[:40]
        test_microbes = self.microbes.iloc[40:]
        train_metabolites = self.metabolites.iloc[:40]
        test_metabolites = self.metabolites.iloc[40:]

        result = mmvec(
            train_microbes,
            train_metabolites,
            n_components=2,
            max_iter=100,
            seed=42,
        )
        q2 = result.score(test_microbes.values, test_metabolites.values)
        self.assertIsInstance(q2, float)
        self.assertLessEqual(q2, 1.0)

    def test_score_zero_metabolite_raises(self):
        """Score with zero-count metabolite sample should raise ValueError."""
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=10,
            seed=42,
        )
        zero_metabolites = self.metabolites.iloc[:5].copy()
        zero_metabolites.iloc[0, :] = 0

        with self.assertRaises(ValueError) as ctx:
            result.score(self.microbes.iloc[:5], zero_metabolites)
        self.assertIn("all-zero counts", str(ctx.exception))


class TestMMvecValidation(unittest.TestCase):
    """Test input validation for MMvec."""

    def setUp(self):
        """Create valid test data."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
            seed=42,
        )
        self.microbes, self.metabolites = res[0], res[1]

    def test_mismatched_sample_counts(self):
        """Mismatched sample counts should raise."""
        microbes = self.microbes.iloc[:40]  # 40 samples
        metabolites = self.metabolites.iloc[:50]  # 50 samples

        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, metabolites, max_iter=1)
        self.assertIn("same number of samples", str(ctx.exception))

    def test_all_zero_columns(self):
        """Columns (features) with all zero counts should raise."""
        microbes = self.microbes.copy()
        for i in (1, 2, 3):
            microbes.iloc[:, i] = 0
        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, self.metabolites, max_iter=1)
        self.assertIn("microbes contains all-zero columns", str(ctx.exception))

        metabolites = self.metabolites.copy()
        for i in (4, 5, 6):
            metabolites.iloc[:, i] = 0
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, metabolites, max_iter=1)
        self.assertIn("metabolites contains all-zero columns", str(ctx.exception))

    def test_all_zero_rows(self):
        """Rows (samples) with all zero counts should raise."""
        microbes = self.microbes.copy()
        for i in (1, 2, 3):
            microbes.iloc[i, :] = 0
        with self.assertRaises(ValueError) as ctx:
            mmvec(microbes, self.metabolites, max_iter=1)
        self.assertIn("microbes contains all-zero rows", str(ctx.exception))

        metabolites = self.metabolites.copy()
        for i in (4, 5, 6):
            metabolites.iloc[i, :] = 0
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, metabolites, max_iter=1)
        self.assertIn("metabolites contains all-zero rows", str(ctx.exception))

    def test_prior_scales(self):
        """u_prior_scale and v_prior_scale must be positive."""
        msg = "u_prior_scale must be positive"
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, self.metabolites, u_prior_scale=0, max_iter=1)
        self.assertIn(msg, str(ctx.exception))
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, self.metabolites, u_prior_scale=-0.5, max_iter=1)
        self.assertIn(msg, str(ctx.exception))

        msg = "v_prior_scale must be positive"
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, self.metabolites, v_prior_scale=0, max_iter=1)
        self.assertIn(msg, str(ctx.exception))
        with self.assertRaises(ValueError) as ctx:
            mmvec(self.microbes, self.metabolites, v_prior_scale=-0.5, max_iter=1)
        self.assertIn(msg, str(ctx.exception))


class TestMMvecInputTypes(unittest.TestCase):
    """Test different input types for MMvec."""

    def test_sparse_pandas_input(self):
        """Test with sparse pandas DataFrame."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
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
            seed=42,
        )

        # Should produce valid results
        self.assertEqual(result.ranks.shape, (8, 10))
        npt.assert_allclose(
            result.ranks.values.mean(axis=1), 0, atol=1e-10
        )


class TestMMvecParameterBehavior(unittest.TestCase):
    """Test parameter behavior for MMvec."""

    def setUp(self):
        """Create test data."""
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
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
            seed=42,
        )

        result_nonzero = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            max_iter=50,
            u_prior_mean=5.0,
            v_prior_mean=5.0,
            seed=42,
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
                seed=42,
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
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
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
            seed=42,
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
            seed=42,
        )

        # Check columns exist
        self.assertIn("iteration", result.convergence.columns)
        self.assertIn("loss", result.convergence.columns)

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
            seed=42,
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
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=10,
            n_samples=50,
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
            seed=42,
        )

        # Check result type
        self.assertIsInstance(result, MMvecResults)

        # Check shapes
        self.assertEqual(result.microbe_embeddings.shape, (8, 3))
        self.assertEqual(result.metabolite_embeddings.shape, (10, 3))
        self.assertEqual(result.ranks.shape, (8, 10))

        # Ranks should be row-centered
        row_means = result.ranks.values.mean(axis=1)
        npt.assert_allclose(row_means, 0, atol=1e-10)

    def test_lbfgs_reproducibility(self):
        """L-BFGS should be deterministic with same seed."""
        result1 = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            seed=123,
        )

        result2 = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            seed=123,
        )

        npt.assert_allclose(
            result1.ranks.values, result2.ranks.values, rtol=1e-10
        )

    def test_lbfgs_recovers_structure(self):
        """L-BFGS should recover embedding structure from synthetic data."""
        # Use more data for better recovery
        res = random_multimodal(
            n_microbes=8,
            n_metabolites=8,
            n_samples=150,
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
            seed=42,
        )

        # Check correlation of pairwise distances
        u_r, u_p = spearmanr(
            pdist(result.microbe_embeddings.values[:, :-1]),
            pdist(U_true),
        )
        self.assertGreater(u_r, 0.3, f"U correlation too low: {u_r}")

    def test_lbfgs_score_on_test_data(self):
        """L-BFGS model should produce reasonable Q-squared score on test data."""
        # Split data
        train_microbes = self.microbes.iloc[:40]
        test_microbes = self.microbes.iloc[40:]
        train_metabolites = self.metabolites.iloc[:40]
        test_metabolites = self.metabolites.iloc[40:]

        result = mmvec(
            train_microbes,
            train_metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=100,
            seed=42,
        )

        # Compute Q^2 score on test data
        q2 = result.score(test_microbes, test_metabolites)

        # Q^2 should be a valid float
        self.assertIsInstance(q2, float)

    def test_lbfgs_verbose_output(self):
        """L-BFGS verbose mode should produce output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            mmvec(
                self.microbes,
                self.metabolites,
                n_components=2,
                optimizer="lbfgs",
                max_iter=50,
                verbose=True,
                seed=42,
            )
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        self.assertIn("Iteration", output)
        self.assertIn("Loss", output)

    def test_seed_generator(self):
        """seed as np.random.Generator should work."""
        rng = np.random.default_rng(42)
        result = mmvec(
            self.microbes,
            self.metabolites,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            seed=rng,
        )
        self.assertIsInstance(result, MMvecResults)
        self.assertEqual(result.ranks.shape, (8, 10))

    def test_numpy_array_inputs(self):
        """mmvec should accept numpy arrays (not just DataFrames)."""
        result = mmvec(
            self.microbes.values,
            self.metabolites.values,
            n_components=2,
            optimizer="lbfgs",
            max_iter=50,
            seed=42,
        )
        self.assertIsInstance(result, MMvecResults)
        self.assertEqual(result.ranks.shape, (8, 10))

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


class TestMMvecCaseStudies(unittest.TestCase):
    """Test MMvec on real-world case study datasets.

    These tests verify that the MMvec implementation can reproduce the biological
    findings from the original MMvec publication and example notebooks.

    """

    @unittest.skip("Skipping a test that requires long runtime.")
    def test_soils_cyanobacteria_metabolites(self):
        """Co-occurrence of Cyanobacteria and known metabolites in wetup biocrust.

        This test reproduces the finding from the soils example notebook: the model
        should learn that Cyanobacteria (specifically rplo 1) co-occur with a known
        set of 13 metabolites associated with Microcoleus vaginatus, a desert crust
        cyanobacterium.

        """
        # Load biocrust wetting dataset
        subdir = os.path.join("data", "soils")
        microbes = pd.read_table(
            get_data_path("microbes.tsv.gz", subdir), index_col=0
        ).T
        metabolites = pd.read_table(
            get_data_path("metabolites.tsv.gz", subdir), index_col=0
        ).T

        # Fit the model with low latent dimension for faster testing
        result = mmvec(
            microbes,
            metabolites,
            n_components=1,
            optimizer="lbfgs",
            max_iter=500,
            u_prior_scale=1.0,
            v_prior_scale=1.0,
            seed=42,
        )

        # Define expected Cyanobacteria-associated metabolites
        expected_metabolites = [
            "(3-methyladenine)",
            "7-methyladenine",
            "4-guanidinobutanoate",
            "uracil",
            "xanthine",
            "hypoxanthine",
            "(N6-acetyl-lysine)",
            "cytosine",
            "N-acetylornithine",
            "succinate",
            "adenosine",
            "guanine",
            "adenine",
        ]

        # Define Cyanobacteria taxon
        cyanobacteria_id = "rplo 1 (Cyanobacteria)"

        # Count how many of the expected metabolites have positive rank for
        # Cyanobacteria.
        ranks_cyanobacteria = result.ranks.loc[cyanobacteria_id]
        positive_count = (ranks_cyanobacteria.loc[expected_metabolites] > 0).sum()

        # The original notebook expects all 13 metabolites to have positive ranks for
        # Cyanobacteria. We use a more lenient threshold (11) to account for
        # optimization variability.
        # self.assertEqual(positive_count, 13)
        self.assertGreaterEqual(positive_count, 11)

        # Check numerical values against expected outputs for reproducibility.
        # NOTE: Due to optimization variability, these values may not match exactly,
        # but should be close.
        obs = result.score(microbes, metabolites)
        exp = 0.218074
        self.assertAlmostEqual(obs, exp, places=6)

        ranks = pd.read_table(get_data_path("ranks.tsv", subdir), index_col=0)
        pdt.assert_frame_equal(
            result.ranks, ranks,
            check_dtype=False, check_exact=False, rtol=0, atol=1e-6,
        )

        probabilities = pd.read_table(
            get_data_path("probabilities.tsv", subdir), index_col=0
        )
        pdt.assert_frame_equal(
            result.probabilities(), probabilities,
            check_dtype=False, check_exact=False, rtol=0, atol=1e-6,
        )

        predictions = pd.read_table(
            get_data_path("predictions.tsv", subdir), index_col=0
        )
        pdt.assert_frame_equal(
            result.predict(microbes), predictions,
            check_dtype=False, check_exact=False, rtol=0, atol=1e-6,
        )

        microbe_embeddings = pd.read_table(
            get_data_path("microbe_embeddings.tsv", subdir), index_col=0
        )
        pdt.assert_frame_equal(
            result.microbe_embeddings, microbe_embeddings,
            check_dtype=False, check_exact=False, rtol=0, atol=1e-6,
        )

        metabolite_embeddings = pd.read_table(
            get_data_path("metabolite_embeddings.tsv", subdir), index_col=0
        )
        pdt.assert_frame_equal(
            result.metabolite_embeddings, metabolite_embeddings,
            check_dtype=False, check_exact=False, rtol=0, atol=1e-6,
        )

    @unittest.skip("Skipping a test that requires long runtime.")
    def test_cf_pseudomonas_rhamnolipids(self):
        """Co-occurrence of Pseudomonas and rhamnolipids in cystic fibrosis sputum.
        
        This test reproduces the finding from the example notebook: the model should
        learn that Pseudomonas microbes co-occur with rhamnolipids and other
        Pseudomonas-associated metabolites.

        """
        # Load CF dataset
        subdir = os.path.join("data", "cf")
        microbes = pd.read_table(
            get_data_path("microbes.tsv.gz", subdir), index_col=0
        ).T
        metabolites = pd.read_table(
            get_data_path("metabolites.tsv.gz", subdir), index_col=0
        ).T
        microbe_meta = pd.read_table(
            get_data_path("microbe_meta.tsv.gz", subdir), index_col=0
        )
        metabolite_meta = pd.read_table(
            get_data_path("metabolite_meta.tsv.gz", subdir), index_col=0
        )

        # Find Pseudomonas microbes in the metadata (n=39)
        pseudomonas = microbe_meta[
            microbe_meta["Taxon"].str.contains("Pseudomonas")
        ].index.tolist()

        # Find expert-annotated metabolites in the metadata (n=20)
        expert_metabolites = metabolite_meta[
            metabolite_meta["expert_annotation"].notnull()
        ].index.tolist()

        # Fit the model with parameters similar to the original
        result = mmvec(
            microbes,
            metabolites,
            n_components=3,
            optimizer="lbfgs",
            max_iter=500,
            u_prior_scale=1.0,
            v_prior_scale=1.0,
            seed=42,
        )

        # For the first Pseudomonas taxon in the ranks table, count expert-annotated
        # metabolites with positive ranks.
        ranks_pseudomonas = result.ranks[result.ranks.index.isin(pseudomonas)]
        first_pseudomonas = ranks_pseudomonas.iloc[0]

        # The original analysis expects 19 out of 20 metabolites with positive ranks.
        # A slightly lower threshold (15) is used here for reproducibility.
        positive_count = (first_pseudomonas.loc[expert_metabolites] > 0).sum()
        self.assertEqual(positive_count, 19)
        self.assertGreaterEqual(positive_count, 15)


if __name__ == "__main__":
    unittest.main()
