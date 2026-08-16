# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from unittest.mock import patch

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
from patsy import DesignMatrix, dmatrix
from scipy.stats import t

from skbio.util import get_data_path
from skbio.stats.composition import clr, rclr
from skbio.stats.composition._ancombc2 import (
    _estimate_params_dense,
    _estimate_params_sparse,
    _lstsq_fit,
    _transform_data,
    _estimate_bias_em,
    _sample_fractions,
    _format_results,
    _calc_statistics,
    _calc_pvalues,
    _adjust_pvalues,
    _init_bias_params,
    _global_test,
    _constrain_est,
    _constrain_est_identity,
    _mdfdr_dunnett,
    struc_zero,
    ancombc,
    ancombc2,
    ANCOMBCResult,
)


"""
This test module uses the HITChip Atlas dataset ("pseq_sub"), adopted and refined from
the official ANCOM-BC tutorial:

- https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/
  ANCOMBC.html

The original dataset was described in:

- Lahti, Leo, et al. "Tipping elements in the human intestinal ecosystem." Nature
  Communications 5.1 (2014): 4344.

A subset of the dataset is used for simplicity and efficiency. We followed the ANCOM-BC
tutorial to preprocess the data and aggregate taxa at the family level. A total of 300
samples were randomly selected to achieve sufficient representation of each category.
The metadata was filtered to retain attributes of interest, including age (continuous),
region (4 categories) and bmi (3 categories), according to the formula used in the
ANCOM-BC tutorial.

The reference output files were generated using the R package ANCOMBC version 2.13.2,
and re-formatted to match the Python output. The R script used for generating the
reference files is provided below.

```R
library(ANCOMBC)

set.seed(42)

table <- read.csv("pseq_sub_feature_table.csv", row.names = 1)
meta <- read.csv("pseq_sub_meta_data.csv", row.names = 1)
meta$bmi <- factor(meta$bmi, levels = c("lean", "overweight", "obese"))

res_bc <- ancombc(
    data = table,
    taxa_are_rows = FALSE,
    meta_data = meta,
    formula = "age + region + bmi",
    group = "bmi",
    p_adj_method = "holm",
    prv_cut = 0,
    lib_cut = 0,
    pseudo = 1,
    tol = 1e-5,
    max_iter = 100,
    conserve = FALSE,
    alpha = 0.05,
    global = TRUE,
    struc_zero = FALSE,
    neg_lb = FALSE,
    n_cl = 1,
    verbose = FALSE
)

write.csv(res_bc$res, "pseq_sub_ancombc_main.csv", row.names = FALSE)
write.csv(res_bc$res_global, "pseq_sub_ancombc_global.csv", row.names = FALSE)

trend_contrast <- list(
    increasing = matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
    decreasing = matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)
)
trend_node <- list(increasing = 2, decreasing = 2)

res_bc2 <- ancombc2(
    data = table,
    taxa_are_rows = FALSE,
    meta_data = meta,
    fix_formula = "age + region + bmi",
    group = "bmi",
    p_adj_method = "holm",
    prv_cut = 0,
    lib_cut = 0,
    pseudo = 0,
    pseudo_sens = TRUE,
    s0_perc = 0.05,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE,
    dunnet = TRUE,
    trend = TRUE,
    trend_control = list(contrast = trend_contrast, node = trend_node, B = 1000),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000),
    iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-05, max_iter = 100),
    struc_zero = FALSE,
    neg_lb = FALSE,
    n_cl = 1,
    verbose = FALSE
)

write.csv(res_bc2$res, "pseq_sub_ancombc2_main.csv", row.names = FALSE)
write.csv(res_bc2$res_global, "pseq_sub_ancombc2_global.csv", row.names = FALSE)
write.csv(res_bc2$res_pair, "pseq_sub_ancombc2_pair.csv", row.names = FALSE)
write.csv(res_bc2$res_dunn, "pseq_sub_ancombc2_dunn.csv", row.names = FALSE)
write.csv(res_bc2$res_trend, "pseq_sub_ancombc2_trend.csv", row.names = FALSE)
```

"""

class CoreTests(TestCase):
    def setUp(self):
        # Example 1 (sparse)
        samples = [f'S{i}' for i in range(1, 7)]
        features = [f'F{i}' for i in range(1, 8)]
        self.data1 = np.array(
            [[ 2,  1,  4,  7,  0,  0,  1],
             [ 1,  0,  0,  6,  5,  1, 10],
             [ 3,  2,  2,  9,  6,  0,  1],
             [ 0, 12,  1,  2,  0,  3,  2],
             [ 2,  8, 27,  0,  0,  7,  3],
             [10,  9,  0,  0,  4,  4,  3]])
        self.table1 = pd.DataFrame(self.data1, index=samples, columns=features)
        grouping = ["well"] * 3 + ["sick"] * 3
        self.meta1 = pd.Series(grouping, index=samples, name="status").to_frame()
        self.dmat1 = np.array(
            [[1, 1],
             [1, 1],
             [1, 1],
             [1, 0],
             [1, 0],
             [1, 0]], dtype=float)

        # Example 2 (dense, old)
        self.data2 = np.array(
            [[12, 11, 10, 10, 10, 10, 10],
             [ 9, 11, 12, 10, 10, 10, 10],
             [ 1, 11, 10, 11, 10,  5,  9],
             [22, 21,  9, 10, 10, 10, 10],
             [20, 22, 10, 10, 13, 10, 10],
             [23, 21, 14, 10, 10, 10, 10]])

        samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
        features = ["b1", "b2", "b3", "b4", "b5", "b6", "b7"]
        self.table2 = pd.DataFrame(self.data2, index=samples, columns=features)
        grouping = ["treatment"] * 3 + ["placebo"] * 3
        self.meta2 = pd.Series(grouping, index=samples, name="group").to_frame()
        self.dmat2 = np.array(
            [[1, 1],
             [1, 1],
             [1, 1],
             [1, 0],
             [1, 0],
             [1, 0]], dtype=float)

    def test_transform_data(self):
        # dense matrix, log transform
        matrix = self.data2.copy()
        obs_data, obs_mask = _transform_data(matrix)
        self.assertIsNone(obs_mask)
        exp_data = np.log(matrix)
        npt.assert_allclose(obs_data, exp_data)

        # output is a new copy
        self.assertIsNot(obs_data.base, matrix)

        # output is float type
        self.assertTrue(np.issubdtype(obs_data.dtype, np.floating))

        # original data is untouched
        npt.assert_array_equal(matrix, self.data2)

        # input is already float
        obs_data, obs_mask = _transform_data(matrix.astype(np.float64))
        npt.assert_allclose(obs_data, exp_data)
        self.assertTrue(obs_data.dtype == np.float64)

        # input is float but not float64; will be kept
        obs_data, obs_mask = _transform_data(matrix.astype(np.float32))
        npt.assert_allclose(obs_data, exp_data)
        self.assertTrue(obs_data.dtype == np.float32)

        # dense matrix, CLR
        obs_data, obs_mask = _transform_data(matrix, center=True)
        exp_data = clr(matrix, axis=0, validate=False)
        npt.assert_allclose(obs_data, exp_data)

        # sparse matrix, with pseudocount
        matrix = self.data1
        exp_mask = matrix == 0
        obs_data, obs_mask = _transform_data(matrix, pseudo=1)
        self.assertIsNone(obs_mask)
        exp_data = np.log(matrix + 1.0)
        npt.assert_allclose(obs_data, exp_data)
        self.assertIsNot(obs_data.base, matrix)
        self.assertTrue(np.issubdtype(obs_data.dtype, np.floating))

        # sparse matrix, log on observed data
        obs_data, obs_mask = _transform_data(matrix)
        npt.assert_array_equal(obs_mask, exp_mask)
        exp_data = np.ones_like(matrix, dtype=float)
        exp_data[~obs_mask] = np.log(matrix[~obs_mask])
        exp_data[obs_mask] = np.nan
        npt.assert_allclose(obs_data, exp_data)

        # sparse matrix, RCLR
        obs_data, obs_mask = _transform_data(matrix, center=True)
        npt.assert_array_equal(obs_mask, exp_mask)
        exp_data = rclr(matrix, axis=0, validate=False)
        npt.assert_allclose(obs_data, exp_data)

    def test_estimate_params_dense(self):
        # NOTE: Numerical accuracy is evaluated up to 5 decimal places. This is because
        # occassionally slightly different results will be generated during the CI
        # workflow. Although SciPy optimizers should be deterministic, this happens in
        # some cases. The initial estimation of parameters is usually precise, but the
        # subsequent iterative optimization is prone to this problem.

        # Example 1 (sparse, +1 before log)
        # By default (groups is True), the full covariance matrix is calculated.
        data_tr = np.log1p(self.data1)
        obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
            data_tr, self.dmat1)
        exp_var = np.array(
            [[0.24374, 0.26827],
             [0.0414 , 0.10443],
             [0.5995 , 0.78012],
             [0.16891, 0.17501],
             [0.17638, 0.36711],
             [0.00728, 0.04665],
             [0.00178, 0.21715]])
        exp_beta = np.array(
            [[ 1.1655 , -0.10615],
             [ 2.35492, -1.75767],
             [ 1.34178, -0.4391 ],
             [ 0.3662 ,  1.74311],
             [ 0.53648,  0.70941],
             [ 1.69172, -1.46068],
             [ 1.2904 , -0.029  ]]).T
        exp_theta = np.array(
            [-0.17616,  0.01642,  0.15975, -0.2722 ,  0.19239,  0.07981])
        exp_cov = np.array(
            [[[ 0.24374, -0.24374], [-0.24374,  0.26827]],
             [[ 0.0414 , -0.0414 ], [-0.0414 ,  0.10443]],
             [[ 0.5995 , -0.5995 ], [-0.5995 ,  0.78012]],
             [[ 0.16891, -0.16891], [-0.16891,  0.17501]],
             [[ 0.17638, -0.17638], [-0.17638,  0.36711]],
             [[ 0.00728, -0.00728], [-0.00728,  0.04665]],
             [[ 0.00178, -0.00178], [-0.00178,  0.21715]]])
        npt.assert_array_equal(obs_var.round(5), exp_var)
        npt.assert_array_equal(obs_beta.round(5), exp_beta)
        npt.assert_array_equal(obs_theta.round(5), exp_theta)
        npt.assert_array_equal(obs_cov.round(5), exp_cov)
        self.assertTrue(obs_var.flags.f_contiguous)

        # Example 2 (dense, just log)
        data_tr = np.log(self.data2)
        obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
            data_tr.copy(), self.dmat2, groups=True)
        exp_var = np.array(
            [[0.00126, 0.27087],
             [0.00036, 0.01524],
             [0.00889, 0.02048],
             [0.00028, 0.02215],
             [0.00495, 0.01983],
             [0.00028, 0.00474],
             [0.00028, 0.00898]])
        exp_beta = np.array(
            [[ 3.07409, -1.51338],
             [ 3.06003, -0.66213],
             [ 2.37962, -0.01626],
             [ 2.30259,  0.03177],
             [ 2.39004, -0.08745],
             [ 2.30259, -0.23105],
             [ 2.30259, -0.03512]]).T
        exp_theta = np.array(
            [ 0.15683,  0.14178, -0.29861, -0.03834,  0.00722,  0.03113])
        exp_cov = np.array(
            [[[ 0.00126, -0.00126], [-0.00126,  0.27087]],
             [[ 0.00036, -0.00036], [-0.00036,  0.01524]],
             [[ 0.00889, -0.00889], [-0.00889,  0.02048]],
             [[ 0.00028, -0.00028], [-0.00028,  0.02215]],
             [[ 0.00495, -0.00495], [-0.00495,  0.01983]],
             [[ 0.00028, -0.00028], [-0.00028,  0.00474]],
             [[ 0.00028, -0.00028], [-0.00028,  0.00898]]])

        npt.assert_array_equal(obs_var.round(5), exp_var)
        npt.assert_array_equal(obs_beta.round(5), exp_beta)
        npt.assert_array_equal(obs_theta.round(5), exp_theta)
        npt.assert_array_equal(obs_cov.round(5), exp_cov)

        exp_var, exp_beta, exp_theta, exp_cov = (
            obs_var, obs_beta, obs_theta, obs_cov)

        # Grouping is None. Only the diagonal of the covariance matrix is calculated,
        # which is sufficient for estimating the same parameters. No covariance is
        # returned.
        obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
            data_tr.copy(), self.dmat2, groups=None)
        self.assertIsNone(obs_cov)
        npt.assert_allclose(obs_var, exp_var)
        npt.assert_allclose(obs_beta, exp_beta)
        npt.assert_allclose(obs_theta, exp_theta)
        self.assertTrue(obs_var.flags.f_contiguous)

        # Grouping specified. Only covariance submatrices of the relevant coefficients
        # will be calculated and returned.
        groups=np.array([1])
        obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
            data_tr.copy(), self.dmat2, groups=groups)
        npt.assert_allclose(obs_var, exp_var)
        npt.assert_allclose(obs_beta, exp_beta)
        npt.assert_allclose(obs_theta, exp_theta)
        self.assertTupleEqual(obs_cov.shape, (
            data_tr.shape[1], len(groups), len(groups)))
        npt.assert_allclose(obs_cov, exp_cov[:, groups][:, :, groups])

    def test_estimate_params_grouped(self):
        # Covariance submatrices from grouped calculation should match the full matrix.
        rng = np.random.default_rng(42)
        data = rng.normal(size=(20, 30))
        dmat = np.column_stack([np.ones(20), rng.normal(size=(20, 5))])
        groups = np.array([1, 3, 4])

        obs = _estimate_params_dense(data.copy(), dmat, groups=groups)
        exp = _estimate_params_dense(data.copy(), dmat, groups=True)

        npt.assert_allclose(obs[0], exp[0])
        npt.assert_allclose(obs[1], exp[1])
        npt.assert_allclose(obs[2], exp[2])
        npt.assert_allclose(obs[3], exp[3][:, groups][:, :, groups])
        self.assertEqual(obs[3].shape, (data.shape[1], 3, 3))

    def test_lstsq_fit(self):
        dmat = np.array(
            [[1, 0, 1],
             [1, 1, 1],
             [1, 2, 1],
             [1, 3, 1]], dtype=float)
        obs_dmat_inv, obs_gram_sum = _lstsq_fit(dmat, gram=True)
        exp_dmat_inv = np.array(
            [[ 0.35,  0.2 ,  0.05, -0.1 ],
             [-0.3 , -0.1 ,  0.1 ,  0.3 ],
             [ 0.35,  0.2 ,  0.05, -0.1 ]])
        exp_gram_sum = np.array([ 0.2, -0.1,  0.2])
        npt.assert_allclose(obs_dmat_inv, exp_dmat_inv)
        npt.assert_allclose(obs_gram_sum, exp_gram_sum)

        # Compare with direct math
        exp_dmat_inv = np.linalg.pinv(dmat)
        exp_gram_sum = np.linalg.pinv(dmat.T @ dmat).sum(axis=0)
        npt.assert_allclose(obs_dmat_inv, exp_dmat_inv)
        npt.assert_allclose(obs_gram_sum, exp_gram_sum)

        # Omit inverse Gram matrix sum
        obs_dmat_inv, obs_gram_sum = _lstsq_fit(dmat, gram=False)
        npt.assert_allclose(obs_dmat_inv, exp_dmat_inv)
        self.assertIsNone(obs_gram_sum)

    def test_estimate_params_unbalanced(self):
        """Unbalanced model and fallback compute test."""
        # An unbalanced, three-covariate model
        data = np.log(np.array(
            [[ 2,  3,  5],
             [ 7, 11, 13],
             [17, 19, 23],
             [29, 31, 37],
             [41, 43, 47]], dtype=float))
        dmat = np.array(
            [[1, 0, 2],
             [1, 1, 1],
             [1, 2, 0],
             [1, 3, 1],
             [1, 5, 4]], dtype=float)
        obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
            data.copy(), dmat)

        # Directly compute per-feature parameters using independent least-squares
        # (slow and fallback) and compare with the optimized code's output.
        exp_beta = np.linalg.pinv(dmat) @ data
        diff = data - dmat @ exp_beta
        exp_theta = np.mean(diff, axis=1)
        gmat_inv = np.linalg.pinv(dmat.T @ dmat)
        exp_cov = np.empty((data.shape[1], dmat.shape[1], dmat.shape[1]))
        for i in range(data.shape[1]):
            exp_cov[i] = (
                gmat_inv @ (dmat.T * (diff[:, i] - exp_theta) ** 2) @ dmat @ gmat_inv
            )
        exp_var = np.diagonal(exp_cov, axis1=1, axis2=2)

        npt.assert_allclose(obs_var, exp_var)
        npt.assert_allclose(obs_beta, exp_beta)
        npt.assert_allclose(obs_theta, exp_theta)
        npt.assert_allclose(obs_cov, exp_cov)

    def test_estimate_params_singular(self):
        """Rank-deficient and underdetermined design tests."""
        # Both cases require a pseudoinverse instead of a regular matrix inverse.
        data = np.log(np.array(
            [[ 2,  3,  5],
             [ 7, 11, 13],
             [17, 19, 23]], dtype=float))

        # Covariates 0 and 2 are identical and constant
        dmat1 = np.array(
            [[1, 0, 1],
             [1, 1, 1],
             [1, 2, 1]], dtype=float)

        # More covariates than samples; covariate 3 is linearly dependent on 0 and 1
        dmat2 = np.array(
            [[1, 0, 1, 2, 0],
             [1, 1, 1, 3, 1],
             [1, 2, 1, 4, 4]], dtype=float)

        for dmat in (dmat1, dmat2):
            obs_var, obs_beta, obs_theta, obs_cov = _estimate_params_dense(
                data.copy(), dmat)
            exp_beta = np.linalg.pinv(dmat) @ data
            diff = data - dmat @ exp_beta
            exp_theta = np.mean(diff, axis=1)

            npt.assert_allclose(obs_beta, exp_beta)
            npt.assert_allclose(obs_theta, exp_theta, atol=1e-14)
            self.assertTrue(np.isfinite(obs_var).all())
            self.assertTrue(np.isfinite(obs_cov).all())

    def test_estimate_params_constant(self):
        """Constant feature abundance test."""
        # A feature with no variation has no estimated residual variance
        data = np.full((4, 1), np.log(5.0))
        dmat = np.array([[1, 0], [1, 1], [1, 3], [1, 4]], dtype=float)
        var_hat, beta, theta, beta_covmat = _estimate_params_dense(data, dmat)

        # Check array shapes and confirm no NaN when there is only one feature.
        self.assertEqual(var_hat.shape, (1, 2))
        self.assertEqual(beta.shape, (2, 1))
        self.assertEqual(theta.shape, (4,))
        self.assertEqual(beta_covmat.shape, (1, 2, 2))
        self.assertTrue(np.isfinite(var_hat).all())
        self.assertTrue(np.isfinite(beta).all())
        self.assertTrue(np.isfinite(theta).all())
        self.assertTrue(np.isfinite(beta_covmat).all())
        npt.assert_allclose(var_hat, 0.0, atol=1e-14)
        npt.assert_allclose(beta_covmat, 0.0, atol=1e-14)

    def test_estimate_params_sparse(self):
        # Example 1 (RCLR transform)
        missing = self.data1 == 0
        data_tr = rclr(self.data1, axis=0, validate=False)
        obs_var_hat, obs_beta, obs_theta, obs_beta_covmat = _estimate_params_sparse(
            data_tr, self.dmat1, missing)
        exp_var_hat = np.array(
            [[0.14683, 0.26361],
             [0.10253, 0.12889],
             [0.33717, 0.4007 ],
             [0.     , 0.04588],
             [0.     , 0.04129],
             [0.02128, 0.0435 ],
             [0.013  , 0.31737]])
        exp_beta = np.array(
            [[ 0.266  , -0.6268 ],
             [ 0.76291, -1.81436],
             [ 0.43704, -0.64821],
             [-0.41512,  0.73585],
             [-0.47541,  0.4906 ],
             [ 0.36923, -1.66446],
             [ 0.09796, -0.19649]]).T
        exp_theta = np.array(
            [-0.17764,  0.18778, -0.00845, -0.54944,  0.28287,  0.26657])
        exp_beta_covmat = np.array(
            [[[ 0.14683, -0.14683], [-0.14683,  0.26361]],
             [[ 0.10253, -0.10253], [-0.10253,  0.12889]],
             [[ 0.33717, -0.33717], [-0.33717,  0.4007 ]],
             [[ 0.     , -0.     ], [-0.     ,  0.04588]],
             [[ 0.     , -0.     ], [-0.     ,  0.04129]],
             [[ 0.02128, -0.02128], [-0.02128,  0.0435 ]],
             [[ 0.013  , -0.013  ], [-0.013  ,  0.31737]]])
        npt.assert_array_equal(obs_var_hat.round(5), exp_var_hat)
        npt.assert_array_equal(obs_beta.round(5), exp_beta)
        npt.assert_array_equal(obs_theta.round(5), exp_theta)
        npt.assert_array_equal(obs_beta_covmat.round(5), exp_beta_covmat)

        # Should match `_estimate_params` on non-zero data
        data_tr = np.log1p(self.data1)
        obs_var_hat, obs_beta, obs_theta, obs_beta_covmat = _estimate_params_sparse(
            data_tr, self.dmat1, np.full(self.data1.shape, False))
        exp_var_hat, exp_beta, exp_theta, exp_beta_covmat = _estimate_params_dense(
            data_tr, self.dmat1)
        npt.assert_allclose(obs_var_hat, exp_var_hat, atol=1e-5)
        npt.assert_allclose(obs_beta, exp_beta, atol=1e-5)
        npt.assert_allclose(obs_theta, exp_theta, atol=1e-5)
        npt.assert_allclose(obs_beta_covmat, exp_beta_covmat, atol=1e-5)

        # Full and diagonal-only covariance paths agree for missing-value data too.
        data_tr = rclr(self.data1, axis=0, validate=False)
        full = _estimate_params_sparse(data_tr, self.dmat1, self.data1 == 0)
        diag = _estimate_params_sparse(
            data_tr, self.dmat1, self.data1 == 0, groups=None
        )
        npt.assert_allclose(diag[0], full[0])
        npt.assert_allclose(diag[1], full[1])
        npt.assert_allclose(diag[2], full[2])
        self.assertIsNone(diag[3])
        self.assertTrue(diag[0].flags.f_contiguous)

    def test_estimate_params_sparse_grouped(self):
        rng = np.random.default_rng(43)
        n_samples, n_features, n_covariates = 20, 25, 6
        data = rng.normal(size=(n_samples, n_features))
        dmat = np.column_stack(
            [np.ones(n_samples), rng.normal(size=(n_samples, n_covariates - 1))]
        )
        zero_mask = rng.random(data.shape) < 0.2
        for j in range(n_features):
            if (~zero_mask[:, j]).sum() < n_covariates + 2:
                zero_mask[: n_covariates + 2, j] = False
        data[zero_mask] = np.nan
        groups = np.array([1, 2, 4])

        for solver in (False, True):
            full = _estimate_params_sparse(
                data.copy(),
                dmat,
                zero_mask,
                max_iter=10,
                groups=True,
                batched=solver,
            )
            subset = _estimate_params_sparse(
                data.copy(),
                dmat,
                zero_mask,
                max_iter=10,
                groups=groups,
                batched=solver,
            )
            npt.assert_allclose(subset[0], full[0])
            npt.assert_allclose(subset[1], full[1])
            npt.assert_allclose(subset[2], full[2])
            npt.assert_allclose(subset[3], full[3][:, groups][:, :, groups])

    def test_estimate_params_sparse_solvers(self):
        """Chunked compact solver agrees with the retained legacy SVD route."""
        data_tr, zero_mask = _transform_data(self.data1.astype(float), 0, True)

        legacy = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, batched="legacy"
        )
        # Exercise several block boundaries, including one feature per SVD.
        for batch_size in (1, 3, None):
            batched = _estimate_params_sparse(
                data_tr,
                self.dmat1,
                zero_mask,
                batched="batched",
                batch=batch_size,
            )
            for observed, expected in zip(batched, legacy):
                npt.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)

        # The diagonal-only covariance route must remain solver-independent too.
        legacy_diag = _estimate_params_sparse(
            data_tr,
            self.dmat1,
            zero_mask,
            batched="legacy",
            groups=None,
        )
        batched_diag = _estimate_params_sparse(
            data_tr,
            self.dmat1,
            zero_mask,
            batched="batched",
            batch=2,
            groups=None,
        )
        for observed, expected in zip(batched_diag[:3], legacy_diag[:3]):
            npt.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)
        self.assertIsNone(legacy_diag[3])
        self.assertIsNone(batched_diag[3])

        # Deliberately near-collinear designs are where a compact X.T X-like
        # application can lose precision. The batched solver detects these features
        # from the masked-design SVD and retains their stable SVD pseudoinverses.
        rng = np.random.default_rng(42)
        n_samples, n_features = 30, 12
        x = rng.normal(size=n_samples)
        dmat = np.column_stack(
            [np.ones(n_samples), x, x + 1e-6 * rng.normal(size=n_samples)]
        )
        data = rng.normal(size=(n_samples, n_features))
        zero_mask = rng.random(data.shape) < 0.2
        data[zero_mask] = np.nan

        legacy = _estimate_params_sparse(
            data,
            dmat,
            zero_mask,
            batched="legacy",
            tol=0.0,
            max_iter=10,
            groups=None,
        )
        batched = _estimate_params_sparse(
            data,
            dmat,
            zero_mask,
            batched="batched",
            batch=4,
            tol=0.0,
            max_iter=10,
            groups=None,
        )
        for observed, expected in zip(batched[:3], legacy[:3]):
            npt.assert_allclose(observed, expected)

    def test_estimate_params_sparse_direct(self):
        # The direct solve reaches the same fixed point as a tightly converged
        # version of the alternating ANCOM-BC2 update.
        data_tr = rclr(self.data1, axis=0, validate=False)
        zero_mask = self.data1 == 0
        observed = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, tol=1e-12, max_iter=1000
        )
        direct = _estimate_params_sparse(data_tr, self.dmat1, zero_mask, direct=True)
        direct_legacy = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, direct=True, batched="legacy"
        )

        for observed_array, direct_array in zip(observed, direct):
            npt.assert_allclose(direct_array, observed_array, atol=1e-10)
        for batched_array, legacy_array in zip(direct, direct_legacy):
            npt.assert_allclose(batched_array, legacy_array, rtol=1e-12, atol=1e-12)

    def test_init_bias_params(self):
        # regular case
        beta = np.array([0.2, 0.75, 1.15, 1.4, 1.85, 2.05, 2.3, 3.2])
        obs = _init_bias_params(beta)
        self.assertTupleEqual(obs, (1.6125, 0.2, 3.2, 1.0, 1.0))

        # no data point is between q1 and q3
        beta = np.array([0, 1])
        obs = _init_bias_params(beta)
        self.assertTupleEqual(obs, (0.5, 0, 1, 1, 1))

        # no data point falls below quantile=0.125
        beta = np.array([0, 0, 0, 0, 1])
        obs = _init_bias_params(beta)
        self.assertTupleEqual(obs, (0, 0, 1, 1, 1))

        # no data point falls above quantile=0.875
        beta = np.array([0, 1, 1, 1, 1])
        obs = _init_bias_params(beta)
        self.assertTupleEqual(obs, (1, 0, 1, 1, 1))

        # variance of data above quantile=0.75 is 0
        beta = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10])
        obs = _init_bias_params(beta)
        self.assertTupleEqual(obs, (0, 0, 10, 1, 1))

    def test_estimate_bias_em(self):
        # Example 1 (sparse): log1p, no NaN
        data = np.log1p(self.data1)
        var_hat, beta, *_ = _estimate_params_dense(data, self.dmat1)
        obs = np.vstack(list(map(_estimate_bias_em, beta, var_hat.T)))
        exp = np.array([[1.28803, 1.28828, 0.00168],
                        [0.06374, 0.06394, 0.02038]])
        npt.assert_array_equal(obs.round(5), exp)

        # RCLR transform, has NaN
        data = rclr(self.data1, axis=0)
        var_hat, beta, *_ = _estimate_params_sparse(data, self.dmat1, self.data1 == 0)
        obs = np.vstack(list(map(_estimate_bias_em, beta, var_hat.T)))
        exp = np.array([[ 0.19357,  0.19545,  0.00043],
                        [-0.50523, -0.50468,  0.0116 ]])
        npt.assert_array_equal(obs.round(5), exp)

        # Example 2 (dense, log)
        data = np.log(self.data2)
        var_hat, beta, *_ = _estimate_params_dense(data, self.dmat2)
        obs = np.vstack(list(map(_estimate_bias_em, beta, var_hat.T)))
        exp = np.array([[ 2.30495,  2.30494,  0.00007],
                        [-0.1981 , -0.19314,  0.00187]])
        npt.assert_array_equal(obs.round(5), exp)

        # CLR transform
        data = clr(self.data2, axis=0)
        var_hat, beta, *_ = _estimate_params_dense(data, self.dmat2)
        obs = np.vstack(list(map(_estimate_bias_em, beta, var_hat.T)))
        exp = np.array([[ 0.00331,  0.01708,  0.00024],
                        [-0.1981 , -0.19314,  0.00187]])
        npt.assert_array_equal(obs.round(5), exp)

    def test_sample_bias(self):
        data = np.log1p(self.table2.to_numpy())
        dmat = dmatrix("group", self.meta2)
        var_hat, beta, _, _ = _estimate_params_dense(data.copy(), dmat)
        bias = np.empty((2, 3))
        for i in range(2):
            bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], max_iter=1)
        delta_em = bias[:, 0]
        # beta_hat = beta.T - delta_em

        obs = _sample_fractions(data, dmat, beta, delta_em)
        exp = np.array(
            [2.43809627, 2.42448053, 2.08291958, 2.36465192, 2.40607366, 2.42865545]
        )
        npt.assert_allclose(obs, exp, atol=1e-5)

    def test_calc_statistics(self):
        data = np.log1p(self.table2.to_numpy())
        dmat = dmatrix("group", self.meta2)
        var_hat, beta, _, _ = _estimate_params_dense(data, dmat)
        bias = np.empty((2, 3))
        for i in range(2):
            res = _estimate_bias_em(beta[i], var_hat[:, i], max_iter=1)
            bias[i] = res
        delta_em = bias[:, 0]
        beta_hat = beta.T - delta_em

        obs = _calc_statistics(beta_hat, var_hat, 0.05, "holm")

        exp_se_hat = np.array([[0.03349832, 0.38489241],
                               [0.01784618, 0.09653226],
                               [0.086538  , 0.12040718],
                               [0.01530208, 0.11952252],
                               [0.06485124, 0.11491594],
                               [0.01530208, 0.07187641],
                               [0.01530208, 0.07062619]])
        exp_W = np.array([[ 2.14805775e+01, -3.06882664e+00],
                          [ 3.95639779e+01, -5.55591315e+00],
                          [ 8.05617300e-01,  5.70547350e-01],
                          [-1.24061800e-01,  9.50716360e-01],
                          [ 1.21029210e+00,  3.69040700e-02],
                          [-1.24061800e-01, -1.63359659e+00],
                          [-1.24061800e-01,  7.48421510e-01]])
        exp_pval = np.array([[2.36547636e-102, 2.14901260e-003],
                             [0.00000000e+000, 2.76164163e-008],
                             [4.20463527e-001, 5.68306514e-001],
                             [9.01266346e-001, 3.41748381e-001],
                             [2.26166816e-001, 9.70561497e-001],
                             [9.01266346e-001, 1.02343585e-001],
                             [9.01266346e-001, 4.54205953e-001]])
        exp_qval = np.array([[1.41928582e-101, 1.28940756e-002],
                             [0.00000000e+000, 1.93314914e-007],
                             [1.00000000e+000, 1.00000000e+000],
                             [1.00000000e+000, 1.00000000e+000],
                             [1.00000000e+000, 1.00000000e+000],
                             [1.00000000e+000, 5.11717926e-001],
                             [1.00000000e+000, 1.00000000e+000]])
        exp_reject = np.array([[ True,  True],
                               [ True,  True],
                               [False, False],
                               [False, False],
                               [False, False],
                               [False, False],
                               [False, False]])

        npt.assert_allclose(obs[0], exp_se_hat, atol=1e-5)
        npt.assert_allclose(obs[1], exp_W, atol=1e-3)
        npt.assert_allclose(obs[2], exp_pval, atol=1e-5)
        npt.assert_allclose(obs[3], exp_qval, atol=1e-5)
        npt.assert_array_equal(obs[4], exp_reject)

        # The memory-efficient DataFrame path should contain the same statistics in the
        # same feature-major / covariate-minor order.
        observed = _format_results(
            beta_hat,
            var_hat,
            self.table2.columns,
            ["Intercept", "group[T.treatment]"],
            0.05,
            "holm",
        )
        self.assertListEqual(
            list(observed.columns),
            ["Log2(FC)", "SE", "W", "pvalue", "qvalue", "Signif"],
        )
        self.assertListEqual(
            observed.index.names, ["FeatureID", "Covariate"]
        )
        npt.assert_allclose(
            observed["Log2(FC)"].to_numpy().reshape(beta_hat.shape), beta_hat
        )
        npt.assert_allclose(
            observed["SE"].to_numpy().reshape(beta_hat.shape), exp_se_hat, atol=1e-5
        )
        npt.assert_allclose(
            observed["W"].to_numpy().reshape(beta_hat.shape), exp_W, atol=1e-3
        )
        npt.assert_allclose(
            observed["pvalue"].to_numpy().reshape(beta_hat.shape), exp_pval, atol=1e-5
        )
        npt.assert_allclose(
            observed["qvalue"].to_numpy().reshape(beta_hat.shape), exp_qval, atol=1e-5
        )
        npt.assert_array_equal(
            observed["Signif"].to_numpy().reshape(beta_hat.shape), exp_reject
        )
        self.assertEqual(str(observed["Signif"].dtype), "boolean")

    def test_calc_statistics_nan_dof(self):
        beta_hat = np.array([[1.0, -1.0], [2.0, -2.0]])
        var_hat = np.ones_like(beta_hat)
        dof = np.array([10.0, np.nan])

        pval = _calc_pvalues(beta_hat, dof)
        npt.assert_array_equal(np.isnan(pval), [[False, False], [True, True]])

        _, _, obs_pval, qval, reject = _calc_statistics(
            beta_hat, var_hat, 0.05, "holm", dof
        )
        npt.assert_array_equal(obs_pval, pval)
        npt.assert_array_equal(np.isnan(qval), [[False, False], [True, True]])
        npt.assert_array_equal(reject, [[False, False], [False, False]])

    def test_adjust_pvalues_nan(self):
        pval = np.array([[0.01, 0.4], [np.nan, 0.2], [0.2, np.nan]])
        obs = _adjust_pvalues(pval, "holm")
        exp = np.array([[0.02, 0.4], [np.nan, 0.4], [0.2, np.nan]])
        npt.assert_allclose(obs, exp, equal_nan=True)

        # The result-construction path reuses the p-value array for q-values. Verify
        # that aliasing ``out`` with the input preserves NaNs and gives the same result.
        work = pval.copy()
        returned = _adjust_pvalues(work, "holm", out=work)
        self.assertIs(returned, work)
        npt.assert_allclose(work, exp, equal_nan=True)

        # A distinct output array should be supported too.
        work = np.empty_like(pval)
        returned = _adjust_pvalues(pval, "holm", out=work)
        self.assertIs(returned, work)
        npt.assert_allclose(work, exp, equal_nan=True)

        # Common aliases use the optimized Benjamini-Hochberg implementation, while
        # other methods continue to fall back to the generic adjustment function.
        pval = np.array([0.01, 0.04, 0.03, 0.2])
        exp = np.array([0.04, 0.05333333, 0.05333333, 0.2])
        for method in ("bh", "fdr_bh", "benjamini-hochberg"):
            obs = _adjust_pvalues(pval, method)
            npt.assert_allclose(obs, exp)

        obs = _adjust_pvalues(pval, "bonferroni")
        npt.assert_allclose(obs, [0.04, 0.16, 0.12, 0.8])

        obs = _adjust_pvalues(pval, None)
        npt.assert_array_equal(obs, pval)

    def test_post_hoc_methods_recalculate(self):
        table = pd.DataFrame(
            np.arange(1, 73, dtype=float).reshape(9, 8),
            index=[f"S{i}" for i in range(9)],
        )
        metadata = pd.DataFrame(
            {
                "grouping": pd.Categorical(["a"] * 3 + ["b"] * 3 + ["c"] * 3),
                "age": np.arange(9, dtype=float),
            },
            index=table.index,
        )
        res = ancombc(
            table, metadata, "grouping + age", grouping="grouping", max_iter=2
        )

        first = res.global_test()
        second = res.global_test()
        self.assertIsNot(first, second)
        pdt.assert_frame_equal(first, second)

    def test_post_hoc_methods_inherit_fit_settings(self):
        table = pd.DataFrame(
            np.arange(1, 73, dtype=float).reshape(9, 8),
            index=[f"S{i}" for i in range(9)],
        )
        metadata = pd.DataFrame(
            {
                "grouping": pd.Categorical(["a"] * 3 + ["b"] * 3 + ["c"] * 3),
                "age": np.arange(9, dtype=float),
            },
            index=table.index,
        )
        res = ancombc(
            table,
            metadata,
            "grouping + age",
            grouping="grouping",
            alpha=0.1,
            p_adjust="bh",
            max_iter=2,
        )

        inherited = res.global_test()
        explicit = res.global_test(alpha=0.1, p_adjust="bh")
        pdt.assert_frame_equal(inherited, explicit)


class AncombcTests(TestCase):

    def test_ancombc(self):
        table = pd.DataFrame(
            [
                [12, 11, 10, 10, 10, 10, 10],
                [9, 11, 12, 10, 10, 10, 10],
                [1, 11, 10, 11, 10, 5, 9],
                [22, 21, 9, 10, 10, 10, 10],
                [20, 22, 10, 10, 13, 10, 10],
                [23, 21, 14, 10, 10, 10, 10],
            ],
            index=["s1", "s2", "s3", "s4", "s5", "s6"],
            columns=["b1", "b2", "b3", "b4", "b5", "b6", "b7"],
        )
        metadata = pd.Series(
            ["treatment", "treatment", "treatment", "placebo", "placebo", "placebo"],
            index=["s1", "s2", "s3", "s4", "s5", "s6"],
            name="grouping",
        ).to_frame()

        # run ANCOM-BC
        res = ancombc(table + 1, metadata, "grouping")

        # check "method" attribute in result
        self.assertEqual(res.method, "ANCOM-BC")

        # check differential abundance of intercept and grouping
        obs = res.res["Signif"].to_numpy()
        exp = np.array([
            [1.0, 1.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
        ]).flatten()
        npt.assert_array_equal(obs, exp)

        # input as numpy array
        res = ancombc(table.to_numpy() + 1, metadata, "grouping")
        obs = res.res["Signif"].to_numpy()
        npt.assert_array_equal(obs, exp)

        # invalid alpha parameter
        for alpha in (-1, 1.1):
            with self.assertRaises(ValueError):
                ancombc(table + 1, metadata, "grouping", alpha=alpha)

    def test_grouping_validation(self):
        table = pd.DataFrame(
            np.arange(1, 73, dtype=float).reshape(9, 8),
            index=[f"S{i}" for i in range(9)],
        )
        metadata = pd.DataFrame(
            {
                "grouping": pd.Categorical(["a"] * 3 + ["b"] * 3 + ["c"] * 3),
                "age": np.arange(9, dtype=float),
                "binary": pd.Categorical(["a"] * 5 + ["b"] * 4),
            },
            index=table.index,
        )

        with self.assertRaisesRegex(ValueError, "not a metadata column"):
            ancombc(table, metadata, "grouping + age", grouping="missing")
        with self.assertRaisesRegex(ValueError, "must be a term in `formula`"):
            ancombc(table, metadata, "age", grouping="grouping")
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            ancombc(table, metadata, "binary + age", grouping="binary")

    # def test_ancombc_pseq(self):
    #     """Test ANCOM-BC on the HITChip Atlas dataset."""
    #     table = pd.read_csv(
    #         get_data_path("raw/pseq_feature_table_subset.csv.gz"), index_col=0
    #     )
    #     meta_data = pd.read_csv(
    #         get_data_path("raw/pseq_meta_data_subset.csv.gz"), index_col=0
    #     )
    #     meta_data = meta_data.dropna(axis=1, how="any")
    #     meta_data["bmi"] = pd.Categorical(
    #         meta_data["bmi"], categories=["obese", "overweight", "lean"]
    #     )

    #     # run ancom-bc for the HITChip Atlas dataset
    #     res = ancombc(table + 1, meta_data, "age + region + bmi").res

    #     # format multi-index dataframe
    #     obs = res["Signif"].unstack()
    #     obs.columns.name = None
    #     obs.index.name = "taxon"
    #     obs = obs.rename(columns={"Intercept": "(Intercept)"})
    #     for c in obs.columns:
    #         obs = obs.rename(columns={c: c.replace("[T.", "").replace("]", "")})

    #     # load ancom-bc results generated by the R package ANCOMBC
    #     exp = pd.read_csv(
    #         get_data_path("raw/pseq_subset_out_res_diff_abn.csv"), index_col="taxon"
    #     ).drop("Unnamed: 0", axis=1)

    #     similarity = exp.eq(obs).sum().sum() / exp.size
    #     npt.assert_equal(similarity, 1.0)

    def test_ancombc_pseq_sub(self):
        wkdir = '/home/drz/Desktop'

        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc(
            table + 1, meta, formula="age + region + bmi", grouping="bmi"
        )
        obs = res.res
        exp = pd.read_table(get_data_path("pseq_sub_ancombc_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp, atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc_main.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc_main.csv", index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp)

        # global test
        obs = res.global_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc_global.tsv"), index_col=0)
        pdt.assert_frame_equal(obs, exp, atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc_global.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc_global.csv", index_col=0)
        pdt.assert_frame_equal(obs, exp)

class Ancombc2Tests(TestCase):
    def setUp(self):
        self.table = pd.DataFrame(
            [
                [12, 11, 10, 10, 10, 10, 10],
                [9, 11, 12, 10, 10, 10, 10],
                [1, 11, 10, 11, 10, 5, 9],
                [22, 21, 9, 10, 10, 10, 10],
                [20, 22, 10, 10, 13, 10, 10],
                [23, 21, 14, 10, 10, 10, 10],
            ],
            index=["s1", "s2", "s3", "s4", "s5", "s6"],
            columns=["b1", "b2", "b3", "b4", "b5", "b6", "b7"],
        )
        self.grouping = pd.Series(
            ["treatment", "treatment", "treatment", "placebo", "placebo", "placebo"],
            index=["s1", "s2", "s3", "s4", "s5", "s6"],
            name="grouping",
        )


    def test_ancombc2(self):
        # ancom-bc2 results of test dataset
        res = ancombc2(self.table, self.grouping.to_frame(), "grouping")
        self.assertEqual(res.method, "ANCOM-BC2")
        self.assertIsInstance(res._dmat, DesignMatrix)
        self.assertEqual(res._dmat.design_info.column_names, ["Intercept", "grouping[T.treatment]"])
        self.assertFalse(res.has_covariance)
        self.assertIsNone(res._vcov_hat)
        self.assertIsNone(res.grouping)
        with self.assertRaisesRegex(ValueError, "grouping=<metadata column>"):
            res.global_test()

        # A two-level factor is valid in the primary model but cannot be selected as
        # the post-hoc grouping, which requires at least three observed groups.
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            ancombc2(
                self.table,
                self.grouping.to_frame(),
                "grouping",
                grouping="grouping",
            )
        obs = res.res["Signif"].to_numpy()

        # expected differential abundance of intercept and grouping
        exp = np.array([
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
        ]).flatten()
        npt.assert_array_equal(obs, exp)

    def test_grouping_controls_posthoc_availability(self):
        res = ancombc2(self.table, self.grouping.to_frame(), "grouping")
        self.assertFalse(res.has_covariance)
        for method in (
            res.global_test, res.pairwise_test, res.dunnett_test, res.trend_test
        ):
            with self.assertRaisesRegex(ValueError, "grouping=<metadata column>"):
                method()

        # Three groups retain only the grouping covariance submatrix.
        table = pd.concat([self.table, self.table.iloc[:3]], ignore_index=True)
        table.index = [f"s{i}" for i in range(len(table))]
        metadata = pd.DataFrame(
            {
                "grouping": pd.Categorical(["a"] * 3 + ["b"] * 3 + ["c"] * 3),
                "age": np.arange(9, dtype=float),
            },
            index=table.index,
        )
        res = ancombc2(
            table, metadata, "grouping + age", grouping="grouping", max_iter=2
        )
        self.assertEqual(res.grouping, "grouping")
        self.assertTrue(res.has_covariance)
        self.assertEqual(res._vcov_hat.shape, (table.shape[1], 2, 2))
        self.assertEqual(res.global_test().shape[0], table.shape[1])
        self.assertEqual(
            res.dunnett_test(bootstraps=5, seed=123)
            .index.get_level_values("FeatureID")
            .nunique(),
            table.shape[1],
        )


    def test_ancombc2_pseq_sub(self):
        wkdir = '/home/drz/Desktop'

        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc2(
            table, meta, formula="age + region + bmi", grouping="bmi"
        )
        obs = res.res
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc2_main.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc2_main.csv", index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp)

        # global test
        obs = res.global_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_global.tsv"), index_col=0)
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc2_global.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc2_global.csv", index_col=0)
        pdt.assert_frame_equal(obs, exp)

        # pairwise test
        obs = res.pairwise_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_pair.tsv"), index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc2_pair.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc2_pair.csv", index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp)

        # dunnett test
        obs = res.dunnett_test(seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_dunn.tsv"), index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc2_dunn.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc2_dunn.csv", index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp)

        # trend test
        obs = res.trend_test(seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_trend.tsv"), index_col=0)
        pdt.assert_frame_equal(obs[["W", "Signif"]], exp[["W", "Signif"]], atol=1e-3)
        # NOTE: Trend test is highly stochastic, therefore we cannot directly compare
        # p- and q-values. See its documentation.

        # obs.to_csv(f"{wkdir}/pseq_sub_ancombc2_trend.csv")
        exp = pd.read_csv(f"{wkdir}/pseq_sub_ancombc2_trend.csv", index_col=0)
        pdt.assert_frame_equal(obs, exp)

    def test_ancombc2_aggregator(self):
        table = self.table
        metadata = self.grouping.to_frame()
        mapping = {
            "b1": "first",
            "b2": "first",
            "b3": "second",
            "b4": "second",
            "b5": "second",
            "b6": "second",
            "b7": "second",
        }
        expected_features = ["first", "second"]

        for aggregator in (
            mapping,
            pd.Series(mapping),
            lambda feature: mapping[feature],
            [mapping[feature] for feature in table.columns],
        ):
            res = ancombc2(
                table,
                metadata,
                "grouping",
                aggregator=aggregator,
            )
            self.assertEqual(
                res.res.index.get_level_values("FeatureID").unique().tolist(),
                expected_features,
            )

        res = ancombc2(
            table,
            metadata,
            "grouping",
            aggregator=["first", "first"] + ["second"] * 5,
        )
        self.assertEqual(
            res.res.index.get_level_values("FeatureID").unique().tolist(),
            expected_features,
        )


class PostHocTests(TestCase):

    def test_constrain_est_identity(self):
        beta_hat = np.array(
            [[-1.0, 2.0, -3.0], [2.0, 1.0, 3.0], [1.0, 3.0, 2.0]]
        )
        contrast = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, -1.0, 1.0]]
        )
        expected = np.array(
            [_constrain_est(beta, np.eye(3), contrast) for beta in beta_hat]
        )

        observed = _constrain_est_identity(beta_hat, contrast)

        npt.assert_allclose(observed, expected, atol=1e-8)

    @patch("skbio.stats.composition._ancombc2._dunn_global")
    def test_mdfdr_dunnett(self, mock_dunn_global):
        mock_dunn_global.return_value = pd.DataFrame(
            {"reject": [True, False, False]}
        )
        W = np.array([[4.0, 0.0], [3.0, 2.0], [1.0, 2.0]])

        obs_pval, obs_qval = _mdfdr_dunnett(
            W=W,
            dof=10.0,
            fwer_ctrl="holm",
            dmat=None,
            group="group",
            bootstraps=100,
            rng=np.random.default_rng(123),
            alpha=0.05,
        )

        exp_pval = 2 * t.sf(np.abs(W[0]), 10.0)
        exp_qval = np.ones_like(W)
        exp_qval[0] = [exp_pval[0] * 6, 1.0]
        npt.assert_allclose(obs_pval[0], exp_pval)
        npt.assert_array_equal(obs_pval[1:], 1.0)
        npt.assert_allclose(obs_qval, exp_qval)

    # def test_global_test(self):
    #     table = pd.read_csv(
    #         get_data_path("raw/pseq_feature_table_subset.csv.gz"), index_col=0
    #     )
    #     meta_data = pd.read_csv(
    #         get_data_path("raw/pseq_meta_data_subset.csv.gz"), index_col=0
    #     )
    #     meta_data = meta_data.dropna(axis=1, how="any")
    #     meta_data["bmi"] = pd.Categorical(
    #         meta_data["bmi"], categories=["obese", "overweight", "lean"]
    #     )
    #     feature_table = np.log1p(table.to_numpy())
    #     dmat = dmatrix("age + region + bmi", meta_data)
    #     covars = dmat.design_info.column_names
    #     n_covars = len(covars)

    #     var_hat, beta, _, vcov_hat = _estimate_params_dense(feature_table, dmat)

    #     bias = np.empty((n_covars, 3))
    #     for i in range(n_covars):
    #         bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], tol=1e-5, max_iter=100)
    #     delta_em = bias[:, 0]

    #     beta_hat = beta.T - delta_em

    #     obs = _global_test(dmat, "bmi", beta_hat, vcov_hat, 0.05, "holm")[-1]
    #     exp = np.array([False,  True, False, False,  True, False, False,  True,  True,
    #                     False, False, False, False, False, False,  True, False, False,
    #                     False, False,  True])
    #     npt.assert_array_equal(obs, exp)


# class StrucZeroTests(TestCase):

#     def test_struc_zero(self):
#         # This test returns all False results (i.e., none of the features have
#         # structural zeros). Please see the doctest for an example that generates both
#         # False and True results. Also, the original (un-subsampled) HITChip Atlas
#         # dataset should produce some True results.
#         table = pd.read_csv(
#             get_data_path("raw/pseq_feature_table_subset.csv.gz"), index_col=0
#         )
#         features = table.columns

#         meta_data = pd.read_csv(
#             get_data_path("raw/pseq_meta_data_subset.csv.gz"), index_col=0
#         )
#         meta_data = meta_data.dropna(axis=1, how="any")
#         categories = ["obese", "overweight", "lean"]
#         meta_data["bmi"] = pd.Categorical(
#             meta_data["bmi"], categories=["obese", "overweight", "lean"]
#         )

#         obs = struc_zero(table, meta_data, "bmi", neg_lb=False)
#         exp = np.zeros((len(features), len(categories)), dtype=bool)

#         # note: groups are sorted alphabetically
#         exp = pd.DataFrame(exp, index=features, columns=["lean", "obese", "overweight"])
#         pdt.assert_frame_equal(obs, exp)

#         obs = struc_zero(table, meta_data, "bmi", neg_lb=True)
#         pdt.assert_frame_equal(obs, exp)


if __name__ == "__main__":
    main()
