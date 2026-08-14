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
from skbio.stats.composition._ancombc2 import (
    _estimate_params,
    _estimate_bias_em,
    _sample_fractions,
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
        samples = [f'S{i}' for i in range(1, 7)]
        features = [f'F{i}' for i in range(1, 8)]
        self.data1 = np.array(
            [[ 2,  1,  4,  7,  0,  0,  1],
             [ 1,  0,  0,  6,  5,  1, 10],
             [ 3,  2,  2,  9,  6,  0,  1],
             [ 0, 12,  1,  2,  0,  3,  2],
             [ 2,  8, 27,  0,  0,  7,  3],
             [10,  9,  0,  0,  4,  4,  3]])
        self.table1 = pd.DataFrame(data, index=samples, columns=features)
        grouping = ["well"] * 3 + ["sick"] * 3
        self.meta1 = pd.Series(grouping, index=samples, name="status").to_frame()

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

    def test_estimate_params(self):
        data = np.log1p(self.table.to_numpy())
        dmat = dmatrix("grouping", self.grouping.to_frame())
        obs = _estimate_params(data, dmat)

        exp_var_hat = [[0.00112214, 0.14814216],
                       [0.00031849, 0.00931848],
                       [0.00748883, 0.01449789],
                       [0.00023415, 0.01428563],
                       [0.00420568, 0.01320567],
                       [0.00023415, 0.00516622],
                       [0.00023415, 0.00498806]]

        exp_beta = [[ 3.11935683,  3.10585971,  2.46951019, 2.39789527,  2.47828263,
                      2.39789527,  2.39789527],
                    [-1.26579628, -0.62095306, -0.01593022, 0.02900379, -0.08038735,
                     -0.20204527, -0.03177006]]

        exp_theta = [0.12293081, 0.10931507, -0.23224588, -0.03514176, 0.00627999,
                     0.02886177]

        # NOTE: atol is set because occassionally slightly different results will be
        # generated during the CI workflow. Although SciPy optimizers should be
        # deterministic, this happens in some cases. The initial estimation of
        # parameters is usually precise, but the subsequent iterative optimization is
        # prone to this problem.
        for o, e in zip(obs[0], exp_var_hat):
            npt.assert_allclose(o, e, atol=1e-5)

        for o, e in zip(obs[1], exp_beta):
            npt.assert_allclose(o, e, atol=1e-5)

        for o, e in zip(obs[2], exp_theta):
            npt.assert_allclose(o, e, atol=1e-5)

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
        data = np.log1p(self.table.to_numpy())
        dmat = dmatrix("grouping", self.grouping.to_frame())
        var_hat, beta, _, _ = _estimate_params(data, dmat)

        obs_0 = _estimate_bias_em(beta[0], var_hat[:, 0], max_iter=100)
        obs_1 = _estimate_bias_em(beta[1], var_hat[:, 1], max_iter=100)
        exp_0 = np.array([2.40007051, 2.4000710, 5.809086e-05])
        exp_1 = np.array([-0.08410937, -0.0847577, 1.395714e-03])

        npt.assert_allclose(obs_0, exp_0, atol=1e-5)
        npt.assert_allclose(obs_1, exp_1, atol=1e-5)

    def test_sample_bias(self):
        data = np.log1p(self.table.to_numpy())
        dmat = dmatrix("grouping", self.grouping.to_frame())
        var_hat, beta, _, _ = _estimate_params(data, dmat)
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
        data = np.log1p(self.table.to_numpy())
        dmat = dmatrix("grouping", self.grouping.to_frame())
        var_hat, beta, _, _ = _estimate_params(data, dmat)
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

    def test_post_hoc_methods_recalculate(self):
        res = ancombc(self.table + 1, self.grouping.to_frame(), "grouping")

        first = res.global_test("grouping")
        second = res.global_test("grouping")
        self.assertIsNot(first, second)
        pdt.assert_frame_equal(first, second)

    def test_post_hoc_methods_inherit_fit_settings(self):
        res = ancombc(
            self.table + 1,
            self.grouping.to_frame(),
            "grouping",
            alpha=0.1,
            p_adjust="bh",
        )

        inherited = res.global_test("grouping")
        explicit = res.global_test("grouping", alpha=0.1, p_adjust="bh")
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
        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc(table + 1, meta, formula="age + region + bmi")
        obs = res.res
        exp = pd.read_table(get_data_path("pseq_sub_ancombc_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp, atol=1e-3)

        # global test
        obs = res.global_test("bmi")
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
        res = ancombc2(self.table, self.grouping.to_frame(), "grouping")
        self.assertEqual(res.method, "ANCOM-BC2")
        self.assertIsInstance(res._dmat, DesignMatrix)
        self.assertEqual(res._dmat.design_info.column_names, ["Intercept", "grouping[T.treatment]"])
        self.assertEqual(res.global_test("grouping").shape[0], self.table.shape[1])
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

    def test_ancombc2_pseq_sub(self):
        cats = ["lean", "overweight", "obese"]
        table = pd.read_csv(get_data_path("pseq_sub_feature_table.csv"), index_col=0)
        meta = pd.read_csv(get_data_path("pseq_sub_meta_data.csv"), index_col=0)
        meta["bmi"] = pd.Categorical(meta["bmi"], categories=cats)

        # core test
        res = ancombc2(table, meta, formula="age + region + bmi")
        obs = res.res
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_main.tsv"), index_col=(0, 1))
        exp["Signif"] = exp["Signif"].astype("boolean")
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # global test
        obs = res.global_test("bmi")
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_global.tsv"), index_col=0)
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # pairwise test
        obs = res.pairwise_test("bmi")
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_pair.tsv"), index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # dunnett test
        obs = res.dunnett_test("bmi", seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_dunn.tsv"), index_col=(0, 1))
        pdt.assert_frame_equal(obs, exp.iloc[:, :-2], atol=1e-3)

        # trend test
        obs = res.trend_test("bmi", seed=123)
        exp = pd.read_table(get_data_path("pseq_sub_ancombc2_trend.tsv"), index_col=0)
        pdt.assert_frame_equal(obs[["W", "Signif"]], exp[["W", "Signif"]], atol=1e-3)
        # NOTE: Trend test is highly stochastic, therefore we cannot directly compare
        # p- and q-values. See its documentation.

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

    #     var_hat, beta, _, vcov_hat = _estimate_params(feature_table, dmat)

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
