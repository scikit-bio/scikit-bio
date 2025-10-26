# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
from patsy import dmatrix

from skbio.util import get_data_path
from skbio.stats.composition._ancombc import (
    _estimate_params,
    _estimate_bias_em,
    _sample_fractions,
    _calc_statistics,
    _init_bias_params,
    _global_test,
    struc_zero,
    ancombc,
)


class AncombcTests(TestCase):
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

        npt.assert_allclose(obs_0, exp_0, atol=1e-2)
        npt.assert_allclose(obs_1, exp_1, atol=1e-2)

    def test_sample_bias(self):
        data = np.log1p(self.table.to_numpy())
        dmat = dmatrix("grouping", self.grouping.to_frame())
        var_hat, beta, _, _ = _estimate_params(data, dmat)
        bias = np.empty((2, 3))
        for i in range(2):
            res = _estimate_bias_em(beta[i], var_hat[:, i], max_iter=1)
            bias[i] = res
        delta_em = bias[:, 0]
        beta_hat = beta.T - delta_em

        obs = _sample_fractions(data, dmat, beta_hat)
        exp = np.array(
            [2.43809627, 2.42448053, 2.08291958, 2.36465192, 2.40607366, 2.42865545]
        )

        npt.assert_allclose(obs, exp, atol=1e-2)

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

        obs = _calc_statistics(beta_hat, var_hat)
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

        exp_p = np.array([[2.36547636e-102, 2.14901260e-003],
                          [0.00000000e+000, 2.76164163e-008],
                          [4.20463527e-001, 5.68306514e-001],
                          [9.01266346e-001, 3.41748381e-001],
                          [2.26166816e-001, 9.70561497e-001],
                          [9.01266346e-001, 1.02343585e-001],
                          [9.01266346e-001, 4.54205953e-001]])

        exp_q = np.array([[1.41928582e-101, 1.28940756e-002],
                          [0.00000000e+000, 1.93314914e-007],
                          [1.00000000e+000, 1.00000000e+000],
                          [1.00000000e+000, 1.00000000e+000],
                          [1.00000000e+000, 1.00000000e+000],
                          [1.00000000e+000, 5.11717926e-001],
                          [1.00000000e+000, 1.00000000e+000]])

        for o, e in zip(obs[0], exp_se_hat):
            npt.assert_allclose(o, e, atol=1e-2)

        for o, e in zip(obs[1], exp_W):
            npt.assert_allclose(o, e, atol=1e-2)

        for o, e in zip(obs[2], exp_p):
            npt.assert_allclose(o, e, atol=1e-2)

        for o, e in zip(obs[3], exp_q):
            npt.assert_allclose(o, e, atol=1e-2)

    def test_ancombc_fail_alpha(self):
        with self.assertRaises(ValueError):
            ancombc(self.table + 1, self.grouping.to_frame(), "grouping", alpha=-1)
        with self.assertRaises(ValueError):
            ancombc(self.table + 1, self.grouping.to_frame(), "grouping", alpha=1.1)

    def test_struc_zero(self):
        # This test returns all False results (i.e., none of the features have
        # structural zeros). Please see the doctest for an example that generates both
        # False and True results. Also, the original (un-subsampled) HITChip Atlas
        # dataset should produce some True results.
        table = pd.read_csv(
            get_data_path("pseq_feature_table_subset.csv.gz"), index_col=0
        )
        features = table.columns

        meta_data = pd.read_csv(
            get_data_path("pseq_meta_data_subset.csv.gz"), index_col=0
        )
        meta_data = meta_data.dropna(axis=1, how="any")
        categories = ["obese", "overweight", "lean"]
        meta_data["bmi"] = pd.Categorical(
            meta_data["bmi"], categories=["obese", "overweight", "lean"]
        )

        obs = struc_zero(table, meta_data, "bmi", neg_lb=False)
        exp = np.zeros((len(features), len(categories)), dtype=bool)

        # note: groups are sorted alphabetically
        exp = pd.DataFrame(exp, index=features, columns=["lean", "obese", "overweight"])
        pdt.assert_frame_equal(obs, exp)

        obs = struc_zero(table, meta_data, "bmi", neg_lb=True)
        pdt.assert_frame_equal(obs, exp)

    def test_global_test(self):
        table = pd.read_csv(
            get_data_path("pseq_feature_table_subset.csv.gz"), index_col=0
        )
        meta_data = pd.read_csv(
            get_data_path("pseq_meta_data_subset.csv.gz"), index_col=0
        )
        meta_data = meta_data.dropna(axis=1, how="any")
        meta_data["bmi"] = pd.Categorical(
            meta_data["bmi"], categories=["obese", "overweight", "lean"]
        )
        feature_table = np.log1p(table.to_numpy())
        dmat = dmatrix("age + region + bmi", meta_data)
        covars = dmat.design_info.column_names
        n_covars = len(covars)

        var_hat, beta, _, vcov_hat = _estimate_params(feature_table, dmat)

        bias = np.empty((n_covars, 3))
        for i in range(n_covars):
            bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], tol=1e-5, max_iter=100)
        delta_em = bias[:, 0]

        beta_hat = beta.T - delta_em

        obs = _global_test(dmat, "bmi", beta_hat, vcov_hat, 0.05, "holm")[-1]
        exp = np.array([False,  True, False, False,  True, False, False,  True,  True,
                        False, False, False, False, False, False,  True, False, False,
                        False, False,  True])
        npt.assert_array_equal(obs, exp)

    def test_ancombc(self):
        # ancom-bc results of test dataset
        res = ancombc(self.table + 1, self.grouping.to_frame(), "grouping")
        obs = res["Signif"].to_numpy()

        # expected differential abundance of intercept and grouping
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

        res = ancombc(self.table.to_numpy() + 1, self.grouping.to_frame(), "grouping")
        obs = res["Signif"].to_numpy()
        npt.assert_array_equal(obs, exp)

        # Load the HITChip Atlas dataset.
        # This dataset is adopted from the official ANCOM-BC tutorial:
        #   https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/
        #   inst/doc/ANCOMBC.html
        # The original dataset was described in:
        #   Lahti, Leo, et al. "Tipping elements in the human intestinal ecosystem."
        #   Nature communications 5.1 (2014): 4344.
        # A subset of the dataset with aggregated features is used in testing for
        # simplicity. We followed the ANCOM-BC tutorial to preprocess the data and
        # aggregate features at the family level. The metadata was filtered to retain
        # attributes of interests, including a continuous covariate of age, and two
        # categorical covariates of region and bmi according to the formula used in
        # the ANCOM-BC tutorial.
        table = pd.read_csv(
            get_data_path("pseq_feature_table_subset.csv.gz"), index_col=0
        )
        meta_data = pd.read_csv(
            get_data_path("pseq_meta_data_subset.csv.gz"), index_col=0
        )
        meta_data = meta_data.dropna(axis=1, how="any")
        meta_data["bmi"] = pd.Categorical(
            meta_data["bmi"], categories=["obese", "overweight", "lean"]
        )

        # run ancom-bc reimplemented in python for the HITChip Atlas dataset
        res = ancombc(table + 1, meta_data, "age + region + bmi")

        # format multi-index dataframe
        obs = res["Signif"].unstack()
        obs.columns.name = None
        obs.index.name = "taxon"
        obs = obs.rename(columns={"Intercept": "(Intercept)"})
        for c in obs.columns:
            obs = obs.rename(columns={c: c.replace("[T.", "").replace("]", "")})

        # load ancom-bc results generated by the R package ANCOMBC
        exp = pd.read_csv(
            get_data_path("pseq_subset_out_res_diff_abn.csv"), index_col="taxon"
        ).drop("Unnamed: 0", axis=1)

        similarity = exp.eq(obs).sum().sum() / exp.size
        npt.assert_equal(similarity, 1.0)

    def test_ancombc_global(self):
        table = pd.read_csv(
            get_data_path("pseq_feature_table_subset.csv.gz"), index_col=0
        )
        meta_data = pd.read_csv(
            get_data_path("pseq_meta_data_subset.csv.gz"), index_col=0
        )
        meta_data = meta_data.dropna(axis=1, how="any")
        meta_data["bmi"] = pd.Categorical(
            meta_data["bmi"], categories=["obese", "overweight", "lean"]
        )

        # run ancom-bc reimplemented in python for the HITChip Atlas dataset
        res = ancombc(table + 1, meta_data, "age + region + bmi", "bmi")
        obs = res[1]["Signif"].to_numpy()
        exp = np.array([False,  True, False, False,  True, False, False,  True,  True,
                        False, False, False, False, False, False,  True, False, False,
                        False, False,  True])
        npt.assert_array_equal(obs, exp)

        meta_data["bmi"] = 1
        with self.assertRaises(ValueError):
            ancombc(table + 1, meta_data, "age + region + bmi", "bmi")


if __name__ == "__main__":
    main()
