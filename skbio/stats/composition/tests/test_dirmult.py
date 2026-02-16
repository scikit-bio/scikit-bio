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

from skbio.stats.composition._dirmult import (
    dirmult_ttest, dirmult_lme,
)


class DirMultTTestTests(TestCase):
    def setUp(self):
        np.random.seed(0)
        # Create sample data for testing
        self.data = {
            'feature1': [5, 8, 12, 15, 20],
            'feature2': [3, 6, 9, 12, 15],
            'feature3': [10, 15, 20, 25, 30],
        }
        self.table = pd.DataFrame(self.data)
        self.grouping = pd.Series(['Group1', 'Group1', 'Group2', 'Group2', 'Group2'])
        self.treatment = 'Group2'
        self.reference = 'Group1'

        d = 50
        n = 200
        self.depth = depth = 1000
        p1 = np.random.lognormal(0, 1, size=d) * 10
        p2 = np.random.lognormal(0.01, 1, size=d) * 10
        self.p1, self.p2 = p1 / p1.sum(), p2 / p2.sum()
        self.data2 = np.vstack((
            np.random.multinomial(depth, self.p1, size=n),
            np.random.multinomial(depth, self.p2, size=n)))
        self.table2 = pd.DataFrame(self.data2)
        self.grouping2 = pd.Series(['Group1'] * n + ['Group2'] * n)

    def test_dirmult_ttest_demo(self):
        # The same example as in doctest
        data = np.array([
            [ 20, 110, 100, 101, 100, 103, 104],
            [ 33, 110, 120, 100, 101, 100, 102],
            [ 12, 110, 100, 110, 100,  50,  90],
            [202, 201,   9,  10,  10,  11,  11],
            [200, 202,  10,  10,  13,  10,  10],
            [203, 201,  14,  10,  10,  13,  12],
        ])
        samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
        features = ["b1", "b2", "b3", "b4", "b5", "b6", "b7"]
        table = pd.DataFrame(data, index=samples, columns=features)
        labels = ["treatment", "treatment", "treatment",
                  "placebo", "placebo", "placebo"]
        grouping = pd.Series(labels, index=samples)

        obs = dirmult_ttest(table, grouping, 'treatment', 'placebo', seed=0)
        self.assertTupleEqual(obs.shape, (7, 7))
        self.assertListEqual(obs.index.to_list(), features)
        exp = {
            "T-statistic": [-17.179, -16.873,  6.943,  6.523,  6.654,  3.84,   7.601],
            "Log2(FC)":    [ -4.992,  -2.534,  1.628,  1.707,  1.528,  1.182,  1.48 ],
            "CI(2.5)":     [ -7.884,  -3.595, -1.048, -0.467, -1.037, -0.703, -0.601],
            "CI(97.5)":    [ -2.293,  -1.462,  4.751,  4.165,  3.978,  3.556,  4.044],
            "pvalue":      [  0.003,   0.001,  0.021,  0.013,  0.019,  0.045,  0.017],
            "qvalue":      [  0.02 ,   0.007,  0.068,  0.066,  0.068,  0.068,  0.068],
        }
        for key, value in exp.items():
            npt.assert_array_equal(obs[key].to_numpy().round(3), np.array(value))
        exp = np.array([True, True, False, False, False, False, False])
        npt.assert_array_equal(obs["Signif"].to_numpy(), exp)

    def test_dirmult_ttest_toy(self):
        p1 = np.array([5, 6, 7])
        p2 = np.array([4, 7, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 1000
        n = 100
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        grouping = pd.Series(['Group1'] * n + ['Group2'] * n)

        exp_lfc = np.log2([4/5, 7/6, 7/7])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        res = dirmult_ttest(table, grouping, self.treatment, self.reference)

        npt.assert_array_less(exp_lfc, res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], exp_lfc)

    def test_dirmult_ttest_toy_depth(self):
        p1 = np.array([5, 6, 7, 8, 9, 4])
        p2 = np.array([4, 7, 7, 6, 5, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 100
        n = 100
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        grouping = pd.Series(['Group1'] * n + ['Group2'] * n)
        exp_lfc = np.log2([4/5, 7/6, 7/7, 6/8, 5/9, 7/4])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates
        res_100 = dirmult_ttest(table, grouping, self.treatment, self.reference)

        # increase sequencing depth by 100 fold
        depth = 10000
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        res_10000 = dirmult_ttest(table, grouping, self.treatment, self.reference)

        # when the sequencing depth increases, the confidence intervals
        # should also shrink
        npt.assert_array_less(res_100['CI(2.5)'], res_10000['CI(2.5)'])
        npt.assert_array_less(res_10000['CI(97.5)'], res_100['CI(97.5)'])

    def test_dirmult_ttest_output(self):
        exp_lfc = np.log2(self.p2 / self.p1)
        exp_lfc = exp_lfc - exp_lfc.mean()
        res = dirmult_ttest(self.table2, self.grouping2,
                            self.treatment, self.reference)

        npt.assert_array_less(res['Log2(FC)'], res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], res['Log2(FC)'])

        # a couple of things that complicate the tests
        # first, there is going to be a little bit of a fudge factor due
        # to the pseudocount, so we will define it via log2(0.5)
        eps = np.abs(np.log2(0.5))

        # second, the confidence interval is expected to be inaccurate
        # for (1/20) of the tests. So we should double check to
        # see if the confidence intervals were able to capture
        # 95% of the log-fold changes correctly
        self.assertGreater(np.mean(res['CI(2.5)'] - eps < exp_lfc), 0.95)
        self.assertGreater(np.mean(res['CI(97.5)'] + eps > exp_lfc), 0.95)

    def test_dirmult_ttest_valid_input(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape[1], 7)  # Expected number of columns
        pdt.assert_index_equal(result.index,
                               pd.Index(['feature1', 'feature2', 'feature3']))

    def test_dirmult_ttest_array_input(self):
        result = dirmult_ttest(self.table.to_numpy(), self.grouping, self.treatment,
                               self.reference)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIsInstance(result.index, pd.RangeIndex)

    def test_dirmult_ttest_no_group(self):
        result = dirmult_ttest(self.table, self.grouping)
        self.assertIsInstance(result, pd.DataFrame)
        result = dirmult_ttest(self.table, self.grouping, treatment=self.treatment)
        self.assertIsInstance(result, pd.DataFrame)
        result = dirmult_ttest(self.table, self.grouping, reference=self.treatment)
        self.assertIsInstance(result, pd.DataFrame)

    def test_dirmult_ttest_no_pseudocount(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference, pseudocount=None)
        self.assertIsInstance(result, pd.DataFrame)

    def test_dirmult_ttest_no_p_adjust(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference, p_adjust=None)
        pdt.assert_series_equal(result['pvalue'], result['qvalue'], check_names=False)

    def test_dirmult_ttest_invalid_table_type(self):
        with self.assertRaises(TypeError):
            dirmult_ttest("invalid_table", self.grouping, self.treatment,
                          self.reference)

    def test_dirmult_ttest_negative_values_in_table(self):
        self.table.iloc[0, 0] = -5  # Modify a value to be negative
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_missing_values_in_grouping(self):
        self.grouping[1] = np.nan  # Introduce a missing value in grouping
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_missing_values_in_table(self):
        self.table.iloc[2, 1] = np.nan  # Introduce a missing value in the table
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_inconsistent_indexes(self):
        self.table.index = ['a', 'b', 'c', 'd', 'e']  # Change table index
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

class DirMultLMETests(TestCase):
    # NOTE: `dirmult_lme` performs numerical optimization, which might (though rarely)
    # generate slightly different results on different platforms. The following tests
    # have specific numbers commented out, just to be safe in the CI workflow. But one
    # may check the accuracy of results locally by restoring the commented code.
    def setUp(self):
        np.random.seed(0)
        index = ["subject1", "subject2", "subject3", "subject4", "subject5", "subject6"]
        columns = ["feature1", "feature2", "feature3", "feature4"]

        # create sample data for testing
        self.table = pd.DataFrame(
            [[20, 110, 100, 101],
             [33, 110, 120, 100],
             [12, 110, 100, 110],
             [202, 201, 9, 10],
             [200, 202, 10, 10],
             [203, 201, 14, 10]],
            index=index,
            columns=columns)

        self.metadata = pd.DataFrame(
            {"Covar1": [1,1,2,2,3,3],
             "Covar2": [1,1,1,1,2,2],
             "Covar3": [1,2,1,2,1,2]},
            index=index)

    def test_dirmult_lme_demo(self):
        # a regular analysis
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak")
        exp = """
  FeatureID Covariate  Reps  Log2(FC)   CI(2.5)  CI(97.5)    pvalue    qvalue  Signif
0  feature1    Covar2     1  2.708376 -0.818699  6.235451  0.132319  0.247129   False
1  feature1    Covar3     1  1.696770 -1.224053  4.617594  0.254876  0.444790   False
2  feature2    Covar2     1  0.956017 -0.736692  2.648726  0.268312  0.464632   False
3  feature2    Covar3     1  0.325451 -0.657935  1.308836  0.516565  0.766291   False
4  feature3    Covar2     1 -1.990268 -4.362246  0.381711  0.100061  0.190110   False
5  feature3    Covar3     1 -0.812892 -2.910615  1.284830  0.447548  0.694797   False
6  feature4    Covar2     1 -1.674125 -3.885791  0.537540  0.137915  0.256810   False
7  feature4    Covar3     1 -1.209329 -3.204871  0.786213  0.234925  0.414660   False
""".strip("\n")
        # self.assertEqual(str(res), exp)
        self.assertIsInstance(res, pd.DataFrame)
        self.assertTupleEqual(res.shape, (8, 9))
        pdt.assert_index_equal(res.columns, pd.Index([
            "FeatureID", "Covariate", "Reps", "Log2(FC)", "CI(2.5)", "CI(97.5)",
            "pvalue", "qvalue", "Signif"]))
        self.assertListEqual(res["FeatureID"].tolist(), [
            "feature1", "feature1", "feature2", "feature2", "feature3", "feature3",
            "feature4", "feature4"])
        self.assertListEqual(res["Covariate"].tolist(), [
            "Covar2", "Covar3", "Covar2", "Covar3", "Covar2", "Covar3", "Covar2",
            "Covar3"])
        self.assertTrue((res["Reps"] == 1).all())
        # npt.assert_array_equal(res["Log2(FC)"].round(5), np.array([
        #     2.70838, 1.69677, 0.95602, 0.32545, -1.99027, -0.81289, -1.67413,
        #     -1.20933]))
        # npt.assert_array_equal(res["CI(2.5)"].round(5), np.array([
        #     -0.8187, -1.22405, -0.73669, -0.65793, -4.36225, -2.91061, -3.88579,
        #     -3.20487]))
        # npt.assert_array_equal(res["CI(97.5)"].round(5), np.array([
        #     6.23545, 4.61759, 2.64873, 1.30884, 0.38171, 1.28483, 0.53754, 0.78621]))
        # npt.assert_array_equal(res["pvalue"].round(5), np.array([
        #     0.13232, 0.25488, 0.26831, 0.51657, 0.10006, 0.44755, 0.13792, 0.23492]))
        # npt.assert_array_equal(res["qvalue"].round(5), np.array([
        #     0.24713, 0.44479, 0.46463, 0.76629, 0.19011, 0.6948 , 0.25681, 0.41466]))
        self.assertTrue((res["Signif"] == False).all())

        # confirm that 2.5% < fold-change < 97.5%
        npt.assert_array_less(res["Log2(FC)"], res["CI(97.5)"])
        npt.assert_array_less(res["CI(2.5)"], res["Log2(FC)"])

    def test_dirmult_lme_array_input(self):
        res = dirmult_lme(
            table=self.table.to_numpy(), metadata=self.metadata,
            formula="Covar2 + Covar3", grouping="Covar1", draws=1, seed=0,
            p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_alt_grouping(self):
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping=self.metadata["Covar1"].to_numpy(), draws=1, seed=0,
            p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_pseudocount(self):
        res = dirmult_lme(
            table=self.table + 0.5, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak", pseudocount=None)
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_intercept(self):
        # "-1" at the end of formula suppresses intercept
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3 - 1",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.988, 0.998, 0.129, 0.761, 0.766, 0.988, 0.971, 0.965]))

    def test_dirmult_lme_re_formula(self):
        # "1": random effect only in intercept (the default scenario)
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak", re_formula="1")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

        # add a random slope for Covar2
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak",
            re_formula="1 + Covar2", fit_method="bfgs")
        self.assertIsInstance(res, pd.DataFrame)
        # NOTE: This test function is prone to generate different results on different
        # platforms. If the following test is to be executed, be prepared that lower
        # precision is needed to get the test pass.
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.274, 0.463, 0.645, 0.766, 0.251, 0.703, 0.400, 0.415]))

    def test_dirmult_lme_vc_formula(self):
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak",
            vc_formula={"Covar2": "0 + C(Covar2)"})
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_p_adjust_and_reml(self):
        result = dirmult_lme(
            self.table, self.metadata, "Covar2 + Covar3", "Covar1",
            draws=1, seed=0, p_adjust=None)
        pdt.assert_series_equal(result['pvalue'], result['qvalue'], check_names=False)

        fit_kwargs = {"reml": False}
        res_ml = dirmult_lme(
            self.table, self.metadata, "Covar2 + Covar3", "Covar1",
            draws=1, seed=0, p_adjust=None, fit_kwargs=fit_kwargs)

        with self.assertRaises(AssertionError):
            npt.assert_allclose(result['CI(97.5)'], res_ml['CI(97.5)'])
            npt.assert_allclose(result['CI(2.5)'], res_ml['CI(2.5)'])
            npt.assert_allclose(result['pvalue'], res_ml['pvalue'])

    def test_dirmult_lme_fit_warnings(self):
        # Convergence warnings are frequently raised during model fitting.
        from statsmodels.tools.sm_exceptions import ConvergenceWarning

        with self.assertWarns(ConvergenceWarning):
            dirmult_lme(table=self.table, metadata=self.metadata,
                        formula="Covar2 + Covar3", grouping="Covar1",
                        draws=1, seed=0, p_adjust="sidak",
                        fit_warnings=True)

    def test_dirmult_lme_fail_all(self):
        # Supply a non-existent optimization method to make it fail.
        msg = "LME fit failed for all features in all replicates."
        with self.assertRaises(ValueError) as cm:
            dirmult_lme(table=self.table, metadata=self.metadata,
                        formula="Covar2 + Covar3", grouping="Covar1",
                        draws=1, seed=0, p_adjust="sidak",
                        fit_method="not_a_method")
        self.assertEqual(str(cm.exception), msg)

    # def test_dirmult_lme_fail_some(self):
    #     # With the BFGS method, LME model fitting will not converge on two features.
    #     # Output will be NaN.
    #     msg = "LME fit failed for 2 features in all replicates, reporting NaNs."
    #     with self.assertWarns(UserWarning) as cm:
    #         res = dirmult_lme(table=self.table, metadata=self.metadata,
    #                           formula="Covar2 + Covar3", grouping="Covar1",
    #                           draws=1, seed=0, p_adjust="sidak",
    #                           fit_method="bfgs", fit_converge=True)
    #     self.assertEqual(str(cm.warning), msg)
    #     self.assertTrue(res.query(
    #         "FeatureID == ['feature1', 'feature3']")["Log2(FC)"].isnull().all())
    #     self.assertTrue(res.query(
    #         "FeatureID == ['feature2', 'feature4']")["Log2(FC)"].notnull().all())

    def test_dirmult_lme_invalid_table_type(self):
        with self.assertRaises(TypeError):
            dirmult_lme("not a table", self.metadata, "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_invalid_metadata_type(self):
        with self.assertRaises(TypeError):
            dirmult_lme(self.table, "not metadata", "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_invalid_grouping(self):
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar2 + Covar3", "hello")

    def test_dirmult_lme_inconsistent_indexes(self):
        # change table index
        self.table.index = ["a", "b", "c", "d", "e", "f"]
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_single_feature(self):
        # only keep the first column of the table
        self.table = self.table.iloc[:, :1]
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar1", "Covar1")

    def test_dirmult_lme_toy_data(self):
        # simulate a dataset of 20 samples by 3 features, with 10 repeated measurements
        # (covar1) and 2 sample groups (covar2)
        p1 = np.array([5, 6, 7])
        p2 = np.array([4, 7, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 1000
        n = 10
        index_range = range(1, n * 2 + 1)
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = ["feature1", "feature2", "feature3"]
        table.index = [f"subject{i}" for i in index_range]

        exp_lfc = np.log2([4/5, 7/6, 7/7])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        metadata = pd.DataFrame({
            "covar1": [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10],
            "covar2": [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]},
            index=[f"subject{i}" for i in index_range])

        res = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        npt.assert_array_equal(res["Log2(FC)"].round(5), [-0.28051, 0.32118, -0.04067])
        npt.assert_array_equal(res["CI(2.5)"].round(5), [-0.41639, 0.18745, -0.16542])
        npt.assert_array_equal(res["CI(97.5)"].round(5), [-0.14359, 0.46473, 0.07706])

        # confirm expected fold change is within confidence interval
        npt.assert_array_less(exp_lfc, res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], exp_lfc)

    def test_dirmult_lme_toy_data_depth(self):
        p1 = np.array([5, 6, 7, 8, 9, 4])
        p2 = np.array([4, 7, 7, 6, 5, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 100
        n = 10
        index_range = range(1, n * 2 + 1)
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = [
            "feature1", "feature2", "feature3", "feature4", "feature5", "feature6"]
        table.index = [f"subject{i}" for i in index_range]

        metadata = pd.DataFrame({
            "covar1": [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10],
            "covar2": [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]},
            index=[f"subject{i}" for i in index_range])

        exp_lfc = np.log2([4/5, 7/6, 7/7, 6/8, 5/9, 7/4])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        res_100 = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        # increase sequencing depth by 100 fold
        depth = 10000
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = [
            "feature1", "feature2", "feature3", "feature4", "feature5", "feature6"]
        table.index = [f"subject{i}" for i in index_range]
        metadata.index = [f"subject{i}" for i in index_range]
        res_10000 = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        # when the sequencing depth increases, the confidence intervals
        # should also shrink
        npt.assert_array_less(res_100['CI(2.5)'], res_10000['CI(2.5)'])
        npt.assert_array_less(res_10000['CI(97.5)'], res_100['CI(97.5)'])


if __name__ == "__main__":
    main()
