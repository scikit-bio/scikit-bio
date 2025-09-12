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
from patsy import dmatrix

from skbio.stats._ancombc import _est_params, _bias_em, _correct_coefficients


class AncombcTests(TestCase):
    def setUp(self):
        self.table = pd.DataFrame([[12, 11, 10, 10, 10, 10, 10],
                                   [9,  11, 12, 10, 10, 10, 10],
                                   [1,  11, 10, 11, 10, 5,  9],
                                   [22, 21, 9,  10, 10, 10, 10],
                                   [20, 22, 10, 10, 13, 10, 10],
                                   [23, 21, 14, 10, 10, 10, 10]],
                                  index=['s1', 's2', 's3', 's4', 's5', 's6'],
                                  columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6',
                                           'b7'])
        self.grouping = pd.Series(['treatment', 'treatment', 'treatment',
                                   'placebo', 'placebo', 'placebo'],
                                  index=['s1', 's2', 's3', 's4', 's5', 's6'])

    def test_est_params(self):
        data = np.log1p(self.table)
        meta = self.grouping.to_frame()
        meta.columns = ["X0"]
        formula = "X0"
        dmat = dmatrix(formula, data=meta, return_type="dataframe")
        obs = _est_params(data.to_numpy(), dmat.to_numpy())

        exp_var_hat = [
            [0.00112214, 0.14814216],
            [0.00031849, 0.00931848],
            [0.00748883, 0.01449789],
            [0.00023415, 0.01428563],
            [0.00420568, 0.01320567],
            [0.00023415, 0.00516622],
            [0.00023415, 0.00498806]
        ]

        exp_beta_hat = [
            [ 3.11935683,  3.10585971,  2.46951019,  2.39789527,  2.47828263, 2.39789527,  2.39789527],
            [-1.26579628, -0.62095306, -0.01593022,  0.02900379, -0.08038735, -0.20204527, -0.03177006]
        ]

        exp_theta = [ 0.12293081,  0.10931507, -0.23224588, -0.03514176,  0.00627999, 0.02886177]

        for o, e in zip(obs[0], exp_var_hat):
            npt.assert_array_equal(o.round(8), e)

        for o, e in zip(obs[1], exp_beta_hat):
            npt.assert_array_equal(o.round(8), e)

        for o, e in zip(obs[2], exp_theta):
            npt.assert_array_equal(o.round(8), e)

    def test_bias_em(self):
        data = np.log1p(self.table)
        meta = self.grouping.to_frame()
        meta.columns = ["X0"]
        formula = "X0"
        dmat = dmatrix(formula, data=meta, return_type="dataframe")
        var_hat, beta_hat, _ = _est_params(data.to_numpy(), dmat.to_numpy())

        obs = _bias_em(beta_hat[0], var_hat[:, 0], max_iter=1)
        exp_bias = np.array([2.39979368e+00,  2.53583125e+00,  6.12616582e-05])

        npt.assert_array_equal(np.array(obs).round(8), exp_bias.round(8))

    def test_correct_coefficients(self):
        data = np.log1p(self.table)
        meta = self.grouping.to_frame()
        meta.columns = ["X0"]
        formula = "X0"
        dmat = dmatrix(formula, data=meta, return_type="dataframe")

        var_hat, beta_hat, _ = _est_params(data.to_numpy(), dmat.to_numpy())
        bias = np.empty((2, 3))
        for i in range(2):
            res = _bias_em(beta_hat[i], var_hat[:, i], max_iter=1)
            bias[i] = res
        delta_em = bias[:, 0]
        obs = _correct_coefficients(beta_hat, delta_em)
        exp_beta = np.array([[ 0.71956315, -1.18116807],
                             [ 0.70606603, -0.53632484],
                             [ 0.06971651,  0.068698  ],
                             [-0.0018984 ,  0.11363201],
                             [ 0.07848895,  0.00424087],
                             [-0.0018984 , -0.11741705],
                             [-0.0018984 ,  0.05285816]])

        for o, e in zip(obs, exp_beta):
            npt.assert_array_equal(o.round(8), e)

if __name__ == '__main__':
    main()
