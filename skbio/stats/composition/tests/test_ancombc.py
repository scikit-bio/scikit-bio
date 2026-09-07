# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from unittest.mock import Mock, patch

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
from patsy import DesignMatrix, dmatrix
from scipy.stats import t

from skbio.util import get_data_path
from skbio.stats.composition import clr, rclr
from skbio.stats.composition._ancombc import (
    _calc_residual,
    _calc_residual_sparse,
    _estimate_params,
    _estimate_params_dense,
    _estimate_params_sparse,
    _lstsq_dense,
    _lstsq_sparse,
    _lstsq_sparse_batch,
    _apply_pinv,
    _calc_variance,
    _calc_covariance,
    _calc_var_cov,
    _adjust_variance,
    _transform_data,
    _validate_grouping,
    _estimate_bias_em,
    _sample_fractions,
    _calc_statistics,
    _calc_pvalues,
    _init_bias_params,
    _r_fit_info,
    _r_rebase_categorical,
    _group_covmat,
    _var_diff,
    _safe_inverse_spd,
    _global_test,
    _global_stats,
    _constrain_est,
    _constrain_est_identity,
    _prep_trend_projection,
    _dunn_global,
    _trend_test,
    _mdfdr_dunnett,
    _mdfdr_pairwise,
    _ancombc_core,
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

        # input is float32; will be kept. Build the expected result from the float32
        # input because the log is intentionally evaluated at float32 precision.
        matrix32 = matrix.astype(np.float32)
        obs_data, obs_mask = _transform_data(matrix32)
        exp_data32 = np.log(matrix32)
        npt.assert_array_equal(obs_data, exp_data32)
        self.assertTrue(obs_data.dtype == np.float32)

        # float16 is promoted to float32 because NumPy linear algebra does not support
        # float16. The transform is then evaluated at float32 precision.
        matrix16 = matrix.astype(np.float16)
        obs_data, obs_mask = _transform_data(matrix16)
        exp_data32 = np.log(matrix16.astype(np.float32))
        npt.assert_array_equal(obs_data, exp_data32)
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

        # Sparse float16 data are likewise promoted to float32. Compare with RCLR on
        # the promoted input so the reference uses the same intended precision.
        matrix16 = matrix.astype(np.float16)
        obs_data, obs_mask = _transform_data(matrix16, center=True)
        npt.assert_array_equal(obs_mask, exp_mask)
        exp_data32 = rclr(matrix16.astype(np.float32), axis=0, validate=False)
        npt.assert_allclose(obs_data, exp_data32, rtol=1e-6, atol=1e-7)
        self.assertTrue(obs_data.dtype == np.float32)

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

    def test_estimate_params_dense_corrected(self):
        # Final fits retain the supplied sampling correction. For two samples per
        # group, the sandwich variances of the group means are 0.5 and 2.
        data = np.array([[1., 2.], [3., 6.], [5., 4.], [7., 8.]])
        dmat = np.array([[1., 0.], [1., 0.], [1., 1.], [1., 1.]])
        expected_cov = np.array([[[0.5, -0.5], [-0.5, 1.0]],
                                 [[2.0, -2.0], [-2.0, 4.0]]])
        for keep_data in (True, False):
            with self.subTest(keep_data=keep_data):
                work = data.copy()
                var, beta, theta, cov, _, _ = _estimate_params(
                    work, dmat, True, None, biased=False, keep_data=keep_data
                )
                npt.assert_allclose(beta, [[2., 4.], [4., 2.]])
                npt.assert_allclose(var, [[0.5, 1.0], [2.0, 4.0]])
                npt.assert_allclose(cov, expected_cov)
                npt.assert_array_equal(theta, np.zeros(4))
                if keep_data:
                    npt.assert_array_equal(work, data)

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

    def test_lstsq_dense(self):
        dmat = np.array(
            [[1, 0, 1],
             [1, 1, 1],
             [1, 2, 1],
             [1, 3, 1]], dtype=float)
        obs_dmat_inv, obs_gram_sum = _lstsq_dense(dmat, gram=True)
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
        obs_dmat_inv, obs_gram_sum = _lstsq_dense(dmat, gram=False)
        npt.assert_allclose(obs_dmat_inv, exp_dmat_inv)
        self.assertIsNone(obs_gram_sum)

    def test_lstsq_sparse(self):
        data = np.log1p(self.data1)
        missing = np.zeros_like(data, dtype=bool)
        dmat_inv, _ = _lstsq_dense(self.dmat1)
        exp_beta = (dmat_inv @ data).T
        exp_theta = np.mean(data - self.dmat1 @ exp_beta.T, axis=1)

        for direct in (False, True):
            obs_theta, obs_beta, _, _ = _lstsq_sparse(
                data, self.dmat1, missing, direct)
            npt.assert_allclose(obs_theta, exp_theta)
            npt.assert_allclose(obs_beta, exp_beta)

    def test_lstsq_sparse_batch(self):
        data = np.log1p(self.data1)
        missing = np.zeros_like(data, dtype=bool)
        dmat_inv, _ = _lstsq_dense(self.dmat1)
        exp_beta = (dmat_inv @ data).T
        exp_theta = np.mean(data - self.dmat1 @ exp_beta.T, axis=1)

        for direct in (False, True):
            obs_theta, obs_beta, _, _ = _lstsq_sparse_batch(
                data, self.dmat1, missing, direct, batch=1
            )
            npt.assert_allclose(obs_theta, exp_theta)
            npt.assert_allclose(obs_beta, exp_beta)

    def test_lstsq_sparse_local_references(self):
        # R/Patsy treatment coding chooses a new local reference when zero omission
        # removes the global reference level but leaves at least two levels. Verify
        # that the R-compatible sparse route reproduces a direct feature-local Patsy
        # fit for multiple additive categorical terms plus a numeric covariate.
        rng = np.random.default_rng(42)
        a = np.repeat(["a0", "a1", "a2"], 12)
        b = np.tile(np.repeat(["b0", "b1", "b2"], 4), 3)
        metadata = pd.DataFrame({"a": a, "b": b, "x": rng.normal(size=36)})
        dmat = dmatrix("a + b + x", metadata)
        data = rng.normal(size=(36, 3))
        missing = np.zeros(data.shape, dtype=bool)
        missing[a == "a0", 0] = True
        missing[(a == "a0") | (b == "b0"), 1] = True
        data[missing] = np.nan

        _, beta, estimable, _ = _lstsq_sparse_batch(
            data, dmat, missing, False, batch=2, biased=False, match_r=True,
        )
        self.assertIsNone(estimable)

        names = dmat.design_info.column_names
        for f in range(data.shape[1]):
            observed = ~missing[:, f]
            local = dmatrix("a + b + x", metadata.loc[observed])
            coef = np.linalg.lstsq(
                np.asarray(local), data[observed, f], rcond=None
            )[0]
            expected = np.zeros(len(names))
            for name, value in zip(local.design_info.column_names, coef):
                expected[names.index(name)] = value
            npt.assert_allclose(beta[f], expected, atol=1e-12)

        # In the first feature, a0 is absent and a1 becomes the local reference. In
        # the second, both a0 and b0 are absent, so a1 and b1 become local references.
        self.assertEqual(beta[0, names.index("a[T.a1]")], 0.0)
        self.assertEqual(beta[1, names.index("a[T.a1]")], 0.0)
        self.assertEqual(beta[1, names.index("b[T.b1]")], 0.0)

    def test_apply_pinv(self):
        data, missing = _transform_data(self.data1.astype(float), 0, True)
        dmat = self.dmat1
        n_feats = data.shape[1]

        # Build compact operators in the same way as the batched sparse route.
        W = 1.0 - missing.T
        X_w = dmat[None, :, :] * W[:, :, None]
        U, S, Vh = np.linalg.svd(X_w, full_matrices=False)
        cutoff = 1e-15 * np.max(S, axis=1, keepdims=True)
        S_inv = np.divide(1.0, S, out=np.zeros_like(S), where=S > cutoff)

        resp = np.zeros_like(data)
        np.copyto(resp, data, where=~missing)

        obs = _apply_pinv(dmat, resp, Vh, S_inv)

        # Reusable workspaces should produce the same result without allocating the
        # intermediate projection on each call.
        rhs = np.empty((dmat.shape[1], n_feats), dtype=data.dtype)
        tmp = np.empty_like(S_inv)
        out = np.empty_like(obs)
        returned = _apply_pinv(
            dmat, resp, Vh, S_inv, out=out, rhs=rhs, illed={}, tmp=tmp
        )
        self.assertIs(returned, out)
        npt.assert_allclose(out, obs)

        exp = np.empty_like(obs)
        for i in range(n_feats):
            x_i = dmat * (~missing[:, i])[:, None]
            exp[i] = np.linalg.pinv(x_i) @ resp[:, i]
        npt.assert_allclose(obs, exp)

        # Stable fallback should override selected rows exactly.
        illed = {0: (np.array([0]), np.zeros((1, dmat.shape[1], dmat.shape[0])))}
        obs_fallback = _apply_pinv(dmat, resp, Vh, S_inv, illed=illed)
        npt.assert_allclose(obs_fallback[0], 0.0)
        npt.assert_allclose(obs_fallback[1:], obs[1:])

    def test_calc_variance(self):
        res2 = np.array(
            [[ 1,  2,  3],
             [ 4,  5,  6],
             [ 7,  8,  9],
             [10, 11, 12]], dtype=float)
        hmat = np.array(
            [[1, 2],
             [3, 4],
             [5, 6],
             [7, 8]], dtype=float)
        obs = _calc_variance(res2, hmat)
        exp = np.sum(res2[:, :, None] * hmat[:, None, :] ** 2, axis=0)
        npt.assert_allclose(obs, exp)
        self.assertTrue(obs.flags.f_contiguous)

    def test_calc_covariance(self):
        res2 = np.array(
            [[ 1,  2,  3],
             [ 4,  5,  6],
             [ 7,  8,  9],
             [10, 11, 12]], dtype=float)
        hmat = np.array(
            [[1, 2],
             [3, 4],
             [5, 6],
             [7, 8]], dtype=float)

        obs = _calc_covariance(res2, hmat)
        exp = np.sum(
            res2[:, :, None, None]
            * hmat[:, None, :, None]
            * hmat[:, None, None, :],
            axis=0,
        )
        npt.assert_allclose(obs, exp)

    def test_calc_var_cov(self):
        res2 = np.array(
            [[ 1,  2,  3],
             [ 4,  5,  6],
             [ 7,  8,  9],
             [10, 11, 12]], dtype=float)
        hmat = np.array(
            [[ 1,  2,  3],
             [ 4,  5,  6],
             [ 7,  8,  9],
             [10, 11, 12]], dtype=float)
        groups = np.array([0, 2])

        obs_var, obs_cov = _calc_var_cov(res2, hmat, groups)
        exp_var = np.sum(res2[:, :, None] * hmat[:, None, :] ** 2, axis=0)
        exp_cov = np.sum(
            res2[:, :, None, None]
            * hmat[:, None, groups, None]
            * hmat[:, None, None, groups],
            axis=0,
        )
        npt.assert_allclose(obs_var, exp_var)
        npt.assert_allclose(obs_cov, exp_cov)
        self.assertTrue(obs_var.flags.f_contiguous)

    def test_adjust_variance(self):
        var_hat = np.array([[1.0, 4.0, 9.0], [16.0, 25.0, 36.0]])
        vcov_hat = np.full((2, 2, 2), -1.0)
        var_delta = np.array([0.25, 1.0, 4.0])
        groups = np.array([2, 0])

        _adjust_variance(var_hat, vcov_hat, var_delta.copy(), 0.5, groups)

        exp_var = np.array([[13.5, 31.5, 69.5], [31.5, 58.5, 108.5]])
        exp_cov = np.array(
            [[[69.5, -1.0], [-1.0, 13.5]], [[108.5, -1.0], [-1.0, 31.5]]]
        )
        npt.assert_allclose(var_hat, exp_var)
        npt.assert_allclose(vcov_hat, exp_cov)

        # A zero quantile disables the variance-stabilizing offset.
        var_hat = np.array([[1.0, 4.0, 9.0], [16.0, 25.0, 36.0]])
        _adjust_variance(var_hat, None, var_delta.copy(), 0)
        exp_var = np.array([[2.25, 9.0, 25.0], [20.25, 36.0, 64.0]])
        npt.assert_allclose(var_hat, exp_var)

        # Missing variances are excluded from the quantile offset.
        var_hat = np.array([[1.0, np.nan, 9.0], [16.0, 25.0, 36.0]])
        _adjust_variance(var_hat, None, var_delta.copy(), 0.5)
        exp_var = np.array([[13.5, np.nan, 69.5], [31.5, 72.0, 108.5]])
        npt.assert_allclose(var_hat, exp_var, equal_nan=True)

        # A full covariance matrix retains every adjusted variance on its diagonal.
        var_hat = np.array([[1.0, 4.0], [9.0, 16.0]])
        vcov_hat = np.zeros((2, 2, 2))
        _adjust_variance(var_hat, vcov_hat, np.array([1.0, 4.0]), 0)
        npt.assert_allclose(np.diagonal(vcov_hat, axis1=1, axis2=2), var_hat)

    def test_statistical_helpers(self):
        groups = np.array([0, 2])
        beta = np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])
        vcov = np.broadcast_to(np.eye(3), (2, 3, 3)).copy()

        npt.assert_array_equal(_group_covmat(vcov, groups), vcov[:, groups][:, :, groups])
        subcov = vcov[:, groups][:, :, groups]
        self.assertIs(_group_covmat(subcov, groups), subcov)
        self.assertEqual(_var_diff([[2.0, 0.5], [0.5, 3.0]]), 4.0)
        self.assertEqual(_var_diff(np.eye(3)), 3.0)

        observed = _global_test(groups, beta, vcov, p_adjust="holm")
        self.assertEqual(observed[0].shape, (2,))
        self.assertEqual(observed[1].shape, (2,))

        W_global, pval = _global_stats(
            groups, beta, vcov, estimable=np.zeros(beta.shape[0], dtype=bool)
        )
        npt.assert_array_equal(np.isnan(W_global), True)
        npt.assert_array_equal(pval, 1.0)

        W = np.array([[1.0, -2.0], [0.5, -0.5]])
        pval, qval = _mdfdr_dunnett(
            W, None, "holm", 4, 0.05, np.random.default_rng(0)
        )
        self.assertEqual(pval.shape, W.shape)
        self.assertEqual(qval.shape, W.shape)

    def test_numerical_fallbacks(self):
        beta = np.array([1.0, -1.0])
        contrast = np.array([[1.0, -1.0]])
        with patch("skbio.stats.composition._ancombc.minimize", side_effect=RuntimeError):
            npt.assert_array_equal(_constrain_est(beta, np.eye(2), contrast), [0.0, 0.0])

        singular = np.array([[1.0, 1.0], [1.0, 1.0]])
        npt.assert_allclose(
            _safe_inverse_spd(singular), np.linalg.inv(singular + 1e-8 * np.eye(2))
        )

    def test_core_parameter_validation(self):
        for kwargs, message in (
            ({"pseudocount": -1}, "Pseudocount"),
            ({"var_quantile": -0.1}, "var_quantile"),
        ):
            with self.subTest(kwargs=kwargs):
                with self.assertRaisesRegex(ValueError, message):
                    ancombc2(self.table2, self.meta2, "group", **kwargs)

        with self.assertRaisesRegex(TypeError, "metadata column name"):
            _validate_grouping(self.meta2, dmatrix("group", self.meta2), 1)

    def test_core_estimability_guards(self):
        n_features, n_covariates = self.data2.shape[1], self.dmat2.shape[1]
        params = (
            np.ones((n_features, n_covariates)),
            np.zeros((n_covariates, n_features)),
            np.zeros(self.data2.shape[0]),
            None,
            np.ones((n_features, n_covariates), dtype=bool),
            np.full(n_features, n_covariates),
        )
        with (
            patch("skbio.stats.composition._ancombc._estimate_params", side_effect=[params, params]),
            patch("skbio.stats.composition._ancombc._estimate_bias_em", return_value=(0., 0., 0.)),
            patch("skbio.stats.composition._ancombc._sample_fractions", return_value=np.zeros(6)),
        ):
            result = ancombc2(self.table2, self.meta2, "group")
        self.assertIsNone(result._estimable)

        failed = (*params[:4], np.zeros((n_features, n_covariates), dtype=bool), params[5])
        with patch("skbio.stats.composition._ancombc._estimate_params", return_value=failed):
            with self.assertRaisesRegex(ValueError, "No estimable features"):
                ancombc2(self.table2, self.meta2, "group")

    def test_calc_residual(self):
        data = np.arange(20, dtype=float).reshape(4, 5)
        dmat = np.array(
            [[1, 0],
             [1, 1],
             [1, 2],
             [1, 3]], dtype=float)
        beta = np.array(
            [[1, 2, 3, 4, 5],
             [2, 3, 4, 5, 6]], dtype=float)

        # Compare result with full-matrix math
        _calc_residual(obs := data.copy(), dmat, beta)
        exp = data - dmat @ beta
        npt.assert_allclose(obs, exp)

        # Change memory size
        exp = obs
        _calc_residual(obs := data.copy(), dmat, beta, target_bytes=32)
        npt.assert_allclose(obs, exp)

        # An empty feature matrix does not require an allocation or calculation.
        empty = np.empty((4, 0))
        _calc_residual(empty, dmat, np.empty((2, 0)))
        self.assertEqual(empty.shape, (4, 0))

    def test_sparse_failed_fit_paths(self):
        metadata = pd.DataFrame({"group": ["a"] * 2 + ["b"] * 2 + ["c"] * 2})
        dmat = dmatrix("group", metadata)
        data = np.arange(18, dtype=float).reshape(6, 3)
        missing = np.zeros(data.shape, dtype=bool)
        missing[:4, 0] = True
        data[missing] = np.nan

        # The first feature has only one observed categorical level, matching R's
        # failed-fit behavior in both sparse solvers and direct fixed-point routes.
        for batch in (None, 2):
            with self.subTest(batch=batch):
                var, beta, theta, cov, estimable, _ = _estimate_params_sparse(
                    data.copy(), dmat, missing, True, direct=True, batch=batch
                )
                npt.assert_array_equal(estimable[0], False)
                npt.assert_allclose(var[0], 0.1 * 6 * _lstsq_dense(dmat, True)[1] ** 2)
                self.assertTrue(np.isfinite(beta[1:]).all())
                self.assertTrue(np.isfinite(theta).all())
                npt.assert_allclose(np.diagonal(cov[0]), var[0])

        # The unbiased route omits sample-effect fitting altogether.
        theta, _, _, _ = _lstsq_sparse(
            data, dmat, missing, direct=False, biased=False
        )
        npt.assert_array_equal(theta, np.zeros(data.shape[0]))

        # A deliberately strict threshold routes the batched direct solver through
        # its stable full-SVD fallback.
        theta, beta, _, _ = _lstsq_sparse_batch(
            data, dmat, missing, direct=True, batch=2, max_cond=1
        )
        self.assertTrue(np.isfinite(theta).all())
        self.assertTrue(np.isfinite(beta[1:]).all())

    def test_r_fit_info_and_rebase_filtering(self):
        metadata = pd.DataFrame({"group": ["a"] * 2 + ["b"] * 2 + ["c"] * 2})
        dmat = dmatrix("group", metadata)
        missing = np.zeros((6, 2), dtype=bool)
        missing[:2, 0] = True
        missing[:4, 1] = True

        valid, rebases = _r_fit_info(dmat, missing)
        npt.assert_array_equal(valid, [True, False])
        self.assertEqual(len(rebases), 1)

        beta = np.array([[1., 2., 3.], [4., 5., 6.]])
        _r_rebase_categorical(beta, rebases, valid)
        npt.assert_allclose(beta[0], [3., 0., 1.])
        npt.assert_allclose(beta[1], [4., 5., 6.])

        # Full-rank categorical coding takes the conservative non-rebasing path.
        full_rank = dmatrix("group - 1", metadata)
        valid, rebases = _r_fit_info(full_rank, missing)
        npt.assert_array_equal(valid, [True, False])
        self.assertEqual(rebases, ())

        # A categorical main effect used in an interaction is checked but not rebased.
        metadata["score"] = np.arange(len(metadata), dtype=float)
        interaction = dmatrix("group * score", metadata)
        valid, rebases = _r_fit_info(interaction, np.zeros((6, 1), dtype=bool))
        self.assertIsNone(valid)
        self.assertEqual(rebases, ())

        # An interaction-only categorical term has no standalone main effect.
        interaction = dmatrix("group:score", metadata)
        valid, rebases = _r_fit_info(interaction, np.zeros((6, 1), dtype=bool))
        self.assertIsNone(valid)
        self.assertEqual(rebases, ())

    def test_trend_test_default_rng(self):
        groups = np.array([0, 1])
        beta = np.array([[1.0, -1.0], [0.5, -0.5]])
        var = np.ones_like(beta)
        vcov = np.broadcast_to(np.eye(2), (2, 2, 2)).copy()

        observed = _trend_test(groups, beta, var, vcov, bootstraps=2)
        self.assertEqual(observed[0].shape, beta.shape)
        self.assertEqual(observed[2].shape, (2,))

    def test_calc_residual_sparse(self):
        data = np.arange(20, dtype=float).reshape(4, 5)
        dmat = np.array(
            [[1, 0],
             [1, 1],
             [1, 2],
             [1, 3]], dtype=float)
        beta = np.array(
            [[1, 2, 3, 4, 5],
             [2, 3, 4, 5, 6]], dtype=float)
        theta = np.array([0.5, -0.5, 1.0, -1.0])

        # Dense matrix: compare result with the dense function
        _calc_residual(exp := data.copy(), dmat, beta)
        exp -= theta[:, None]
        np.square(exp, out=exp)
        obs = _calc_residual_sparse(
            data, dmat, beta.T, theta, np.zeros(data.shape, dtype=bool)
        )
        npt.assert_allclose(obs, exp)

        # Sparse matrix: missing residuals are zeroed after squaring.
        missing = np.zeros(data.shape, dtype=bool)
        missing[[0, 2, 3], [1, 3, 4]] = True
        data[missing] = np.nan
        _calc_residual(exp := data.copy(), dmat, beta)
        exp -= theta[:, None]
        np.square(exp, out=exp)
        exp[missing] = 0.0
        obs = _calc_residual_sparse(data, dmat, beta.T, theta, missing)
        npt.assert_allclose(obs, exp)

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
        obs_var_hat, obs_beta, obs_theta, obs_covmat, obs_esti, obs_rank = (
            _estimate_params_sparse(data_tr, self.dmat1, missing)
        )
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
        exp_covmat = np.array(
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
        npt.assert_array_equal(obs_covmat.round(5), exp_covmat)
        self.assertIsNone(obs_esti)
        npt.assert_array_equal(obs_rank, [2] * 7)

        # Should match `_estimate_params_dense` on non-zero data
        data_tr = np.log1p(self.data1)
        obs_var_hat, obs_beta, obs_theta, obs_covmat, _, _ = (
            _estimate_params_sparse(
                data_tr, self.dmat1, np.full(self.data1.shape, False)
            )
        )
        exp_var_hat, exp_beta, exp_theta, exp_covmat = _estimate_params_dense(
            data_tr, self.dmat1)
        npt.assert_allclose(obs_var_hat, exp_var_hat, atol=1e-5)
        npt.assert_allclose(obs_beta, exp_beta, atol=1e-5)
        npt.assert_allclose(obs_theta, exp_theta, atol=1e-5)
        npt.assert_allclose(obs_covmat, exp_covmat, atol=1e-5)

        # Full and diagonal-only covariance paths agree for missing-value data too.
        data_tr = rclr(self.data1, axis=0, validate=False)
        full = _estimate_params_sparse(data_tr, self.dmat1, self.data1 == 0)
        diag = _estimate_params_sparse(
            data_tr, self.dmat1, self.data1 == 0, None
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
                data.copy(), dmat, zero_mask, True, batch=solver, max_iter=10
            )
            subset = _estimate_params_sparse(
                data.copy(), dmat, zero_mask, groups, batch=solver, max_iter=10
            )
            npt.assert_allclose(subset[0], full[0])
            npt.assert_allclose(subset[1], full[1])
            npt.assert_allclose(subset[2], full[2])
            npt.assert_allclose(subset[3], full[3][:, groups][:, :, groups])

    def test_estimate_params_sparse_solvers(self):
        """Chunked compact solver agrees with the retained legacy SVD route."""
        data_tr, zero_mask = _transform_data(self.data1.astype(float), 0, True)

        legacy = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, batch=None
        )
        # Exercise several block boundaries, including one feature per SVD.
        for batch_size in (1, 3, None):
            batched = _estimate_params_sparse(
                data_tr, self.dmat1, zero_mask, batch=batch_size
            )
            for observed, expected in zip(batched[:4], legacy[:4]):
                if observed is None:
                    self.assertIsNone(expected)
                else:
                    npt.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)
            npt.assert_array_equal(batched[5], legacy[5])

        # The diagonal-only covariance route must remain solver-independent too.
        legacy_diag = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, None, batch=None
        )
        batched_diag = _estimate_params_sparse(
            data_tr, self.dmat1, zero_mask, None, batch=2
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
            data, dmat, zero_mask, None, batch=None, tol=0.0, max_iter=10
        )
        batched = _estimate_params_sparse(
            data, dmat, zero_mask, None, batch=4, tol=0.0, max_iter=10
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
            data_tr, self.dmat1, zero_mask, direct=True, batch=None
        )

        for observed_array, direct_array in zip(observed[:4], direct[:4]):
            if observed_array is None:
                self.assertIsNone(direct_array)
            else:
                npt.assert_allclose(direct_array, observed_array, atol=1e-10)
        for batched_array, legacy_array in zip(direct[:4], direct_legacy[:4]):
            if batched_array is None:
                self.assertIsNone(legacy_array)
            else:
                npt.assert_allclose(
                    batched_array, legacy_array, rtol=1e-12, atol=1e-12
                )
        npt.assert_array_equal(direct[5], direct_legacy[5])

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
            bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], max_iter=1)
        delta_em = bias[:, 0]
        beta_hat = beta.T - delta_em

        obs = _calc_statistics(beta_hat, var_hat, 0.05, "holm")

        exp_lfc = beta_hat
        exp_se = np.array([[0.03349832, 0.38489241],
                           [0.01784618, 0.09653226],
                           [0.086538  , 0.12040718],
                           [0.01530208, 0.11952252],
                           [0.06485124, 0.11491594],
                           [0.01530208, 0.07187641],
                           [0.01530208, 0.07062619]])
        exp_W = np.array([[ 21.4805775, -3.06882664],
                          [ 39.5639779, -5.55591315],
                          [  0.8056173,  0.57054735],
                          [ -0.1240618,  0.95071636],
                          [  1.2102921,  0.03690407],
                          [ -0.1240618, -1.63359659],
                          [ -0.1240618,  0.74842151]])
        exp_pval = np.array([[2.37e-102, 0.002149 ],
                             [0.       , 2.7616e-8],
                             [0.4204635, 0.5683065],
                             [0.9012663, 0.3417484],
                             [0.2261668, 0.9705615],
                             [0.9012663, 0.1023436],
                             [0.9012663, 0.454206 ]])
        exp_qval = np.array([[1.4e-101, 1.289e-2],
                             [0.      , 1.933e-7],
                             [1.      , 1.      ],
                             [1.      , 1.      ],
                             [1.      , 1.      ],
                             [1.      , 0.511718],
                             [1.      , 1.      ]])
        exp_reject = np.array([[ True,  True],
                               [ True,  True],
                               [False, False],
                               [False, False],
                               [False, False],
                               [False, False],
                               [False, False]])

        npt.assert_allclose(obs[0], exp_lfc)
        npt.assert_allclose(obs[1], exp_se, atol=1e-5)
        npt.assert_allclose(obs[2], exp_W, atol=1e-5)
        npt.assert_allclose(obs[3], exp_pval, atol=1e-5)
        npt.assert_allclose(obs[4], exp_qval, atol=1e-5)
        npt.assert_array_equal(obs[5], exp_reject)

        # dof has NaN
        beta_hat = np.array([[1.0, -1.0], [2.0, -2.0]])
        var_hat = np.ones_like(beta_hat)
        dof = np.array([10.0, np.nan])
        _, _, _, pval, qval, reject = _calc_statistics(
            beta_hat, var_hat, 0.05, "holm", dof
        )
        exp = [[False, False], [True, True]]
        npt.assert_array_equal(np.isnan(pval), exp)
        npt.assert_array_equal(np.isnan(qval), exp)
        npt.assert_array_equal(reject, np.full((2, 2), False))

    def test_calc_statistics_inplace(self):
        for dtype in (np.float32, np.float64):
            for order in ("C", "F"):
                beta = np.array([[1., -2.], [3., 4.]], dtype=dtype, order=order)
                var = np.array([[4., 9.], [16., 25.]], dtype=dtype, order=order)
                for mask in (None, np.array([[True, False], [True, True]])):
                    expected = _calc_statistics(beta, var, .05, "holm",
                                                estimable=mask)
                    npt.assert_array_equal(beta, [[1., -2.], [3., 4.]])
                    npt.assert_array_equal(var, [[4., 9.], [16., 25.]])
                    work_beta, work_var = beta.copy(order=order), var.copy(order=order)
                    observed = _calc_statistics(work_beta, work_var, .05, "holm",
                                                estimable=mask, inplace=True)
                    self.assertIs(observed[0], work_beta)
                    self.assertIs(observed[1], work_var)
                    for obs, exp in zip(observed, expected):
                        npt.assert_array_equal(obs, exp)

    def test_posthoc_retention(self):
        rng = np.random.default_rng(12)
        metadata = pd.DataFrame({"group": ["a"] * 4 + ["b"] * 4 + ["c"] * 4})
        table = rng.integers(1, 20, size=(12, 20))
        for fit in (ancombc, ancombc2):
            plain = fit(table, metadata, "group")
            grouped = fit(table, metadata, "group", grouping="group")
            pdt.assert_frame_equal(plain.result, grouped.result)
            for name in ("_beta_hat", "_var_hat", "_vcov_hat", "_dmat", "_dof",
                         "_estimable", "_groups", "_global_cache"):
                self.assertIsNone(getattr(plain, name))
            for name in ("global_test", "pairwise_test", "dunnett_test", "trend_test"):
                with self.assertRaisesRegex(ValueError, "requires a post-hoc grouping"):
                    getattr(plain, name)()
            # Editing the primary table cannot change retained model arrays.
            beta, var = grouped._beta_hat.copy(), grouped._var_hat.copy()
            expected = grouped.pairwise_test()
            grouped.result.loc[:, "Log(FC)"] = 0.
            grouped.result.loc[:, "SE"] = 0.
            npt.assert_array_equal(grouped._beta_hat, beta)
            npt.assert_array_equal(grouped._var_hat, var)
            pdt.assert_frame_equal(grouped.pairwise_test(), expected)

    def test_calc_pvalues(self):
        W = np.array([[0.0, -1.0], [2.0, -3.0]])

        obs = _calc_pvalues(W)
        exp = np.array([[1.0, 0.3173105], [0.0455003, 0.0026998]])
        npt.assert_array_equal(obs.round(7), exp)

        dof = np.array([5.0, 10.0])
        obs = _calc_pvalues(W, dof)
        exp = np.array([[1.0, 0.3632175], [0.073388, 0.0133437]])
        npt.assert_array_equal(obs.round(7), exp)

        dof = 5.0
        obs = _calc_pvalues(W, dof)
        exp = np.array([[1.0, 0.3632175], [0.1019395, 0.0300992]])
        npt.assert_array_equal(obs.round(7), exp)

        dof = np.array([5.0, np.nan])
        obs = _calc_pvalues(W, dof)
        exp = np.array([[1.0, 0.3632175], [np.nan, np.nan]])
        npt.assert_array_equal(obs.round(7), exp)

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

    def test_global_stats_cache(self):
        table = np.arange(1, 73, dtype=float).reshape(9, 8)
        metadata = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3 + ["c"] * 3})
        for fit in (ancombc, ancombc2):
            fitted = fit(table, metadata, "group", grouping="group", max_iter=2)
            kwargs = {key: getattr(fitted, key) for key in fitted._private_defaults}

            def fresh():
                return ANCOMBCResult(fitted.result, fitted._method, **kwargs)

            for methods in (("global_test", "pairwise_test"),
                            ("pairwise_test", "global_test")):
                calls = [(method, alpha, adjust)
                         for alpha, adjust in ((0.05, "holm"), (0.5, "bh"),
                                               (0.1, None))
                         for method in methods]
                expected = [getattr(fresh(), method)(alpha=alpha, p_adjust=adjust)
                            for method, alpha, adjust in calls]
                result = fresh()
                self.assertIsNone(result._global_cache)
                with patch("skbio.stats.composition._ancombc._global_stats",
                           wraps=_global_stats) as calculate:
                    for (method, alpha, adjust), exp in zip(calls, expected):
                        obs = getattr(result, method)(alpha=alpha, p_adjust=adjust)
                        pdt.assert_frame_equal(obs, exp)
                        # User edits must not change cached statistics or later calls.
                        obs.loc[:, "W"] = -999.
                        obs.loc[:, "pvalue"] = -999.
                        again = getattr(result, method)(alpha=alpha, p_adjust=adjust)
                        pdt.assert_frame_equal(again, exp)
                    calculate.assert_called_once()

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

    def test_validate_grouping(self):
        # metadata with a 3-category, a 2-category and a numeric column
        metadata = pd.DataFrame({
            "group": pd.Categorical(["a"] * 3 + ["b"] * 3 + ["c"] * 3),
            "binary": pd.Categorical(["a"] * 5 + ["b"] * 4),
            "score": np.array([0, 0, 1, 1, 2, 2, 0, 1, 2], dtype=float)})

        # normal case: 3-category column
        dmat = dmatrix("group", metadata)
        obs = _validate_grouping(metadata, dmat, "group")
        npt.assert_array_equal(obs, [1, 2])

        # not in the formula
        with self.assertRaisesRegex(ValueError, "must be a term in"):
            _validate_grouping(metadata, dmat, "score")

        # numeric column is prohibited
        dmat = dmatrix("score", metadata)
        with self.assertRaisesRegex(ValueError, "at least two covariates"):
            _validate_grouping(metadata, dmat, "score")

        # 2-category column is prohibited (post-hoc analysis needs at least 3)
        dmat = dmatrix("binary", metadata)
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            _validate_grouping(metadata, dmat, "binary")

        # tricky case: 2-category column without intercept (will have 2 indices in
        # design matrix but is still prohibited)
        dmat = dmatrix("binary - 1", metadata)
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            _validate_grouping(metadata, dmat, "binary")

        # plain string column is okay
        metadata["group"] = metadata["group"].astype(object)
        dmat = dmatrix("group", metadata)
        obs = _validate_grouping(metadata, dmat, "group")
        npt.assert_array_equal(obs, [1, 2])

        # complex formula
        dmat = dmatrix("binary * group + score", metadata)
        obs = _validate_grouping(metadata, dmat, "group")
        npt.assert_array_equal(obs, [2, 3])
        with self.assertRaisesRegex(ValueError, "at least two covariates"):
            _validate_grouping(metadata, dmat, "score")
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            _validate_grouping(metadata, dmat, "binary")

        # numeric column cast into factor
        metadata["score"] = pd.Categorical(metadata["score"].astype(int))
        dmat = dmatrix("score", metadata)
        obs = _validate_grouping(metadata, dmat, "score")
        npt.assert_array_equal(obs, [1, 2])


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
        self.assertEqual(res._method, "ANCOM-BC")

        # The result object presents itself as the primary DataFrame while retaining an
        # explicit ``result`` attribute.
        pdt.assert_series_equal(res["qvalue"], res.result["qvalue"])
        selected = res[res["qvalue"] <= 0.05]
        expected = res.result[res.result["qvalue"] <= 0.05]
        pdt.assert_frame_equal(selected, expected)
        self.assertEqual(repr(res), repr(res.result))
        self.assertEqual(res._repr_html_(), res.result._repr_html_())

        # check differential abundance of intercept and grouping
        obs = res.result["Signif"].to_numpy()
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
        obs = res.result["Signif"].to_numpy()
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

    def test_ancombc_pseq_sub(self):
        """Test on the HITChip Atlas dataset."""
        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc(
            table + 1, meta, formula="age + region + bmi", grouping="bmi"
        )
        obs = res.result
        exp = pd.read_table(get_data_path("pseq_sub_ancombc_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        exp.rename(columns={"Log2(FC)": "Log(FC)"}, inplace=True)
        pdt.assert_frame_equal(obs, exp, atol=1e-3)

        # global test
        obs = res.global_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc_global.tsv"), index_col=0)
        pdt.assert_frame_equal(obs, exp, atol=1e-3)


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
        table, grouping = self.table, self.grouping.to_frame()

        res = ancombc2(table, grouping, "grouping")
        self.assertEqual(res._method, "ANCOM-BC2")
        self.assertIsNone(res._dmat)
        self.assertIsNone(res._vcov_hat)
        with self.assertRaisesRegex(ValueError, "requires a post-hoc grouping"):
            res.global_test()

        # A two-level factor is valid in the primary model but cannot be selected as
        # the post-hoc grouping, which requires at least three observed groups.
        with self.assertRaisesRegex(ValueError, "at least three observed groups"):
            ancombc2(table, grouping, "grouping", grouping="grouping")

        for var_quantile in (-0.1, 1.1):
            with self.assertRaisesRegex(ValueError, "`var_quantile`"):
                ancombc2(table, grouping, "grouping", var_quantile=var_quantile)

        obs = res.result["Signif"].to_numpy()

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

    def test_rank_deficient(self):
        # With pseudocount=0, each feature is fitted only on samples where it is
        # observed. A feature confined to one level of a categorical factor therefore
        # cannot necessarily identify every model coefficient.
        samples = [f"s{i}" for i in range(12)]
        metadata = pd.DataFrame(
            {"group": pd.Categorical(["a"] * 6 + ["b"] * 6)}, index=samples
        )
        rng = np.random.default_rng(4)
        common = rng.integers(2, 20, size=(12, 8))
        ref_only = np.array([5, 7, 6, 8, 9, 5] + [0] * 6)[:, None]
        nonref_only = np.array([0] * 6 + [4, 8, 6, 7, 9, 5])[:, None]
        features = [f"f{i}" for i in range(8)] + ["ref_only", "nonref_only"]
        table = pd.DataFrame(
            np.hstack((common, ref_only, nonref_only)), index=samples, columns=features
        )

        # The masked-design SVD identifies coefficient estimability essentially for
        # free while retaining finite Moore-Penrose coefficients for fitted values.
        data, missing = _transform_data(table.to_numpy(), center=True)
        dmat = dmatrix("group", metadata)
        for batch in (None, 1):
            var, beta, theta, cov, estimable, rank = _estimate_params_sparse(
                data, dmat, missing, groups=None, batch=batch, match_r=False
            )
            self.assertIsNone(cov)
            self.assertTrue(np.isfinite(var).all())
            self.assertTrue(np.isfinite(beta).all())
            self.assertTrue(np.isfinite(theta).all())
            npt.assert_array_equal(estimable[:-2], True)
            npt.assert_array_equal(estimable[-2:], [[True, False], [False, False]])
            npt.assert_array_equal(rank, [2] * 8 + [1, 1])

        res = _ancombc_core(table, metadata, "group", v2=True, match_r=False)

        self.assertIsNone(res._beta_hat)
        self.assertIsNone(res._estimable)

        # Only inferential output is suppressed. The reference-group intercept remains
        # estimable for a feature observed only in that group, whereas its group effect
        # does not. A feature observed only in the non-reference group cannot separate
        # the intercept from the group coefficient, so neither is estimable.
        obs = res.result.loc[["ref_only", "nonref_only"]]
        self.assertTrue(np.isfinite(obs.loc[("ref_only", "Intercept"), "Log(FC)"]))
        invalid = [
            ("ref_only", "group[T.b]"),
            ("nonref_only", "Intercept"),
            ("nonref_only", "group[T.b]"),
        ]
        for idx in invalid:
            row = obs.loc[idx]
            self.assertTrue(row[["Log(FC)", "SE", "W"]].isna().all())
            self.assertEqual(row["pvalue"], 1.0)
            self.assertEqual(row["qvalue"], 1.0)
            self.assertFalse(row["Signif"])

        # R's lm cannot construct contrasts when zero omission leaves only one factor
        # level, so the default compatibility mode suppresses the entire feature fit.
        res_r = ancombc2(table, metadata, "group")
        obs_r = res_r.result.loc[["ref_only", "nonref_only"]]
        self.assertTrue(obs_r[["Log(FC)", "SE", "W"]].isna().all().all())
        npt.assert_array_equal(obs_r["pvalue"], 1.0)
        npt.assert_array_equal(obs_r["qvalue"], 1.0)
        self.assertFalse(obs_r["Signif"].any())

        # Presence/absence evidence is deliberately separate from ordinary regression
        # inference: structural-zero analysis identifies the missing group directly.
        zero = struc_zero(table, metadata, "group")
        self.assertTrue(zero.loc["ref_only", "b"])
        self.assertTrue(zero.loc["nonref_only", "a"])

    def test_rank_deficient_posthoc(self):
        # A feature absent from one of three groups has one aliased grouping
        # coefficient. Whole-group post-hoc hypotheses are therefore unavailable even
        # though the remaining primary coefficients may still be estimable.
        samples = [f"s{i}" for i in range(15)]
        metadata = pd.DataFrame(
            {"group": pd.Categorical(["a"] * 5 + ["b"] * 5 + ["c"] * 5)},
            index=samples,
        )
        rng = np.random.default_rng(5)
        common = rng.integers(2, 20, size=(15, 8))
        partial = np.array(
            [4, 5, 6, 7, 5, 8, 9, 7, 6, 8, 0, 0, 0, 0, 0]
        )[:, None]
        table = pd.DataFrame(
            np.hstack((common, partial)),
            index=samples,
            columns=[f"f{i}" for i in range(8)] + ["partial"],
        )

        # R does not fail this feature wholesale: two factor levels remain observed.
        # The unused third-level dummy is omitted by lm and ANCOMBC fills its global
        # coefficient slot with zero. The default compatibility mode mirrors that.
        res_r = ancombc2(table, metadata, "group", grouping="group")
        self.assertIsNone(res_r._estimable)
        self.assertAlmostEqual(
            res_r.result.loc[("partial", "group[T.c]"), "Log(FC)"], 0.0
        )
        self.assertTrue(np.isfinite(res_r.global_test().loc["partial", "W"]))

        res = _ancombc_core(
            table, metadata, "group", v2=True, grouping="group", match_r=False
        )

        npt.assert_array_equal(res._estimable[-1], [True, True, False])

        row = res.global_test().loc["partial"]
        self.assertTrue(np.isnan(row["W"]))
        self.assertEqual(row["pvalue"], 1.0)
        self.assertEqual(row["qvalue"], 1.0)
        self.assertFalse(row["Signif"])

        for result in (
            res.pairwise_test().loc["partial"],
            res.dunnett_test(bootstraps=5, seed=1).loc["partial"],
        ):
            self.assertTrue(result[["Log(FC)", "SE", "W"]].isna().all().all())
            npt.assert_array_equal(result["pvalue"], 1.0)
            npt.assert_array_equal(result["qvalue"], 1.0)
            self.assertFalse(result["Signif"].any())

        row = res.trend_test(bootstraps=5, seed=1).loc["partial"]
        self.assertTrue(np.isnan(row["W"]))
        self.assertEqual(row["pvalue"], 1.0)
        self.assertEqual(row["qvalue"], 1.0)
        self.assertFalse(row["Signif"])

    def test_grouping_controls_posthoc_availability(self):
        res = ancombc2(self.table, self.grouping.to_frame(), "grouping")
        self.assertIsNone(res._vcov_hat)
        for method in (
            res.global_test, res.pairwise_test, res.dunnett_test, res.trend_test
        ):
            with self.assertRaisesRegex(ValueError, "requires a post-hoc grouping"):
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
        self.assertIsNotNone(res._vcov_hat)
        self.assertEqual(res._vcov_hat.shape, (table.shape[1], 2, 2))
        self.assertEqual(res.global_test().shape[0], table.shape[1])
        self.assertEqual(
            res.dunnett_test(bootstraps=5, seed=123)
            .index.get_level_values("FeatureID")
            .nunique(),
            table.shape[1],
        )

    def test_ancombc2_pseq_sub(self):
        """Test on the HITChip Atlas dataset."""
        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc2(table, meta, formula="age + region + bmi", grouping="bmi")
        obs = res.result
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        exp.rename(columns={"Log2(FC)": "Log(FC)"}, inplace=True)
        exp_main = exp.iloc[:, :-2]
        pdt.assert_frame_equal(obs, exp_main, atol=1e-3)

        # R may keep a fit when one level of a multi-level factor is absent. In these
        # two features, region still has other observed levels, so the missing US dummy
        # is filled with zero rather than causing a whole-feature failure.
        for feature in ("Cyanobacteria", "Spirochaetes"):
            self.assertAlmostEqual(
                obs.loc[(feature, "region[T.US]"), "Log(FC)"], 0.0
            )

        # global test
        obs = res.global_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_global.tsv"), index_col=0)
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # pairwise test
        obs = res.pairwise_test()
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_pair.tsv"), index_col=(0, 1))
        exp.rename(columns={"Log2(FC)": "Log(FC)"}, inplace=True)
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # dunnett test
        obs = res.dunnett_test(seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_dunn.tsv"), index_col=(0, 1))
        exp.rename(columns={"Log2(FC)": "Log(FC)"}, inplace=True)
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # trend test
        obs = res.trend_test(seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_trend.tsv"), index_col=0)
        pdt.assert_frame_equal(obs[["W", "Signif"]], exp[["W", "Signif"]], atol=1e-3)
        # NOTE: Trend test is highly stochastic, therefore we cannot directly compare
        # p- and q-values. See its documentation.

    def test_ancombc2_sensitivity(self):
        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        fits = [
            ancombc2(table, meta, "age + region + bmi", "bmi", pseudocount=pseudo)
            for pseudo in (0, 0.1, 0.5, 1)
        ]

        def sensitivity(results):
            result = results[0].copy()
            signif = pd.concat([x["Signif"] for x in results], axis=1)
            result["Pass"] = signif.eq(signif.iloc[:, 0], axis=0).all(axis=1)
            result["Robust"] = result["Signif"] & result["Pass"]
            return result

        # Primary result
        obs = sensitivity([fit.result for fit in fits])
        exp = pd.read_table(
            get_data_path("pseq_sub_ancombc2_main.tsv"), index_col=(0, 1)
        )
        pdt.assert_frame_equal(
            obs[["Pass", "Robust"]], exp[["Pass", "Robust"]], check_dtype=False
        )

        # Global test
        global_results = [fit.global_test() for fit in fits]
        obs_global = sensitivity(global_results)
        exp = pd.read_table(
            get_data_path("pseq_sub_ancombc2_global.tsv"), index_col=0
        )
        pdt.assert_frame_equal(
            obs_global[["Pass", "Robust"]],
            exp[["Pass", "Robust"]],
            check_dtype=False,
        )

        # Pairwise directional test
        obs = sensitivity([fit.pairwise_test() for fit in fits])
        exp = pd.read_table(
            get_data_path("pseq_sub_ancombc2_pair.tsv"), index_col=(0, 1)
        )
        pdt.assert_frame_equal(
            obs[["Pass", "Robust"]], exp[["Pass", "Robust"]], check_dtype=False
        )

        # Dunnett's test. Reuse the same seed so pseudo-count sensitivity is not
        # confounded with Monte Carlo variation.
        obs = sensitivity([fit.dunnett_test(seed=123) for fit in fits])
        exp = pd.read_table(
            get_data_path("pseq_sub_ancombc2_dunn.tsv"), index_col=(0, 1)
        )
        pdt.assert_frame_equal(
            obs[["Pass", "Robust"]], exp[["Pass", "Robust"]], check_dtype=False
        )

        # The R implementation assigns the global sensitivity decision to the trend
        # test rather than rerunning the stochastic trend test at each pseudo-count.
        obs = fits[0].trend_test(seed=123)
        obs["Pass"] = obs_global["Pass"]
        obs["Robust"] = obs["Signif"] & obs["Pass"]
        exp = pd.read_table(
            get_data_path("pseq_sub_ancombc2_trend.tsv"), index_col=0
        )
        pdt.assert_frame_equal(
            obs[["Pass", "Robust"]], exp[["Pass", "Robust"]], check_dtype=False
        )


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
                res.result.index.get_level_values("FeatureID").unique().tolist(),
                expected_features,
            )

        res = ancombc2(
            table,
            metadata,
            "grouping",
            aggregator=["first", "first"] + ["second"] * 5,
        )
        self.assertEqual(
            res.result.index.get_level_values("FeatureID").unique().tolist(),
            expected_features,
        )


class PostHocTests(TestCase):

    def test_dunn_bootstrap_counts(self):
        # Ties do not count as exceedances; invalid features remain nonsignificant.
        draws = np.array([[[1., 0.], [3., 0.], [1., 1.]],
                          [[2., 0.], [2., 0.], [0., 0.]],
                          [[0., 0.], [4., 0.], [3., 3.]]])
        W = np.array([[1., 0.], [2., 0.], [np.nan, np.nan]])
        for dof in (None, 10., np.array([10., 20., np.nan])):
            rng = Mock()
            rng.standard_normal.side_effect = draws.copy()
            rng.standard_t.side_effect = draws.copy()
            obs = _dunn_global(W, 3, dof, None, 0.05, rng,
                               estimable=np.array([True, True, False]))
            npt.assert_array_equal(obs["p_val"], [1 / 3, 2 / 3, 1.])
            npt.assert_array_equal(obs["q_val"], obs["p_val"])
            npt.assert_array_equal(obs["reject"], False)

    def test_trend_bootstrap_counts(self):
        # Projection onto the positive orthant gives exact ties at zero.
        draws = np.array([[[0., 0.], [1., 2.], [-1., -1.]],
                          [[1., 2.], [0., 0.], [0., 0.]],
                          [[-1., -1.], [1., 2.], [-1., -1.]]])
        for estimable in (None, np.array([True, True, False])):
            rng = Mock()
            rng.standard_normal.side_effect = (
                draws.copy() if estimable is None else draws[:, estimable].copy()
            )
            _, _, _, obs_pval, obs_qval, _ = _trend_test(
                np.arange(2), np.zeros((3, 2)), np.ones((3, 2)),
                np.tile(np.eye(2), (3, 1, 1)), p_adjust=None,
                trend_contrast={"positive": np.eye(2)},
                trend_node={"positive": 1}, bootstraps=3, rng=rng,
                estimable=estimable,
            )
            expected = [1 / 3, 2 / 3, 0. if estimable is None else 1.]
            npt.assert_array_equal(obs_pval, expected)
            npt.assert_array_equal(obs_qval, expected)

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

    def test_prepared_trend_projection(self):
        beta = np.array([[-1., 2.], [3., -4.], [0., 0.]])
        # Redundant constraints exercise singular active-set Gram matrices.
        for contrast in (np.eye(2), np.array([[1., 0.], [1., 0.], [0., 1.]])):
            projections = _prep_trend_projection(contrast)
            with patch("numpy.linalg.pinv", side_effect=AssertionError(
                "Prepared projections must not repeat the decomposition"
            )):
                for data in (beta, -beta):
                    obs = _constrain_est_identity(data, contrast, projections)
                    npt.assert_allclose(obs, np.maximum(data, 0), atol=1e-14)

        # With no constraints, every coefficient vector is already feasible.
        contrast = np.empty((0, 2))
        projections = _prep_trend_projection(contrast)
        npt.assert_array_equal(
            _constrain_est_identity(beta, contrast, projections), beta
        )

        # Large systems retain the existing SLSQP fallback.
        contrast = np.tile(np.eye(2), (6, 1))
        self.assertIsNone(_prep_trend_projection(contrast))
        npt.assert_allclose(
            _constrain_est_identity(beta, contrast), np.maximum(beta, 0), atol=1e-8
        )

    @patch("skbio.stats.composition._ancombc._dunn_global")
    def test_mdfdr_dunnett(self, mock_dunn_global):
        mock_dunn_global.return_value = pd.DataFrame(
            {"reject": [True, False, False]}
        )
        W = np.array([[4.0, 0.0], [3.0, 2.0], [1.0, 2.0]])

        obs_pval, obs_qval = _mdfdr_dunnett(
            W=W,
            dof=10.0,
            fwer_ctrl="holm",
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

    def test_mdfdr_fractional_families(self):
        # Five features, two screened, three comparisons: effective size = 7.5.
        selected = np.array([True, True, False, False, False])
        W = np.tile(t.isf(np.array([0.01, 0.02, 0.03]) / 2, 10), (5, 1))
        expected = {
            "holm": [0.075, 0.13, 0.165],
            "bonf": [0.075, 0.15, 0.225],
            "bh": [0.075] * 3,
            "by": [0.1944642857142857] * 3,
        }
        for method, exp in expected.items():
            with patch("skbio.stats.composition._ancombc._global_test",
                       return_value=(None, None, None, selected)):
                pair = _mdfdr_pairwise(
                    W, 10., method, None, None, None, None, 0.05)
            with patch("skbio.stats.composition._ancombc._dunn_global",
                       return_value=pd.DataFrame({"reject": selected})):
                dunn = _mdfdr_dunnett(
                    W, 10., method, 1, 0.05, np.random.default_rng(0))
            for pval, qval in (pair, dunn):
                npt.assert_allclose(qval[selected], np.tile(exp, (2, 1)))
                npt.assert_array_equal(qval[~selected], 1.)
                npt.assert_array_equal(pval[~selected], 1.)

    def test_mdfdr_inflated_families(self):
        from statsmodels.stats.multitest import multipletests

        W = np.array([[4., 2.], [3., 1.], [1., 2.]])
        for selected in ([True, False, False], [True, True, False],
                         [True, True, True], [False, False, False]):
            selected = np.array(selected)
            for method, sm_method in (("holm", "holm"), ("bh", "fdr_bh"),
                                      ("bonf", "bonferroni"), ("by", "fdr_by")):
                with patch("skbio.stats.composition._ancombc._global_test",
                           return_value=(None, None, None, selected)):
                    pair = _mdfdr_pairwise(
                        W, 10., method, None, None, None, None, 0.05)
                with patch("skbio.stats.composition._ancombc._dunn_global",
                           return_value=pd.DataFrame({"reject": selected})):
                    dunn = _mdfdr_dunnett(
                        W, 10., method, 1, 0.05, np.random.default_rng(0))
                for pval, qval in (pair, dunn):
                    exp = np.ones_like(W)
                    if selected.any():
                        n_tests = W.shape[1] * W.shape[0] // selected.sum()
                        for i, col in enumerate(pval):
                            padded = np.pad(col, (0, n_tests - col.size),
                                            constant_values=1.)
                            exp[i] = multipletests(padded, method=sm_method)[1][:2]
                    npt.assert_allclose(qval, exp, rtol=1e-14, atol=0)
                    npt.assert_array_equal(pval[~selected], 1.)


class StrucZeroTests(TestCase):
    def test_struc_zero(self):
        samples = [f'S{i}' for i in range(1, 7)]
        features = [f'F{i}' for i in range(1, 8)]
        data = np.array(
            [[ 2,  1,  4,  7,  0,  0,  1],
             [ 1,  0,  0,  6,  5,  0, 10],
             [ 3,  2,  2,  9,  6,  0,  1],
             [ 0, 12,  1,  0,  0,  3,  2],
             [ 2,  8, 27,  0,  0,  7,  3],
             [10,  9,  0,  0,  4,  4,  3]])
        table = pd.DataFrame(data, index=samples, columns=features)
        grouping = ["well"] * 3 + ["sick"] * 3
        meta = pd.Series(grouping, index=samples, name="status").to_frame()
        obs = struc_zero(table, meta, "status")
        exp = pd.DataFrame(np.array(
            [[0, 0, 0, 1, 0, 0, 0],
             [0, 0, 0, 0, 0, 1, 0]], dtype=bool).T,
            index=features, columns=["sick", "well"])  # sorted alphabetically
        pdt.assert_frame_equal(obs, exp)

    def test_struc_zero_pseq_sub(self):
        """Test on the HITChip Atlas dataset."""
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        cats = ["lean", "overweight", "obese"]
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # Groups are sorted alphabetically in the result.
        exp = pd.DataFrame(False, index=table.columns, columns=sorted(cats))
        obs = struc_zero(table, meta, "bmi")
        pdt.assert_frame_equal(obs, exp)

        # Use negative lower bound to detect structural zeros in the "overweight" group.
        exp.loc["Cyanobacteria", "overweight"] = True
        obs = struc_zero(table, meta, "bmi", neg_lb=True)
        pdt.assert_frame_equal(obs, exp)


if __name__ == "__main__":
    main()
