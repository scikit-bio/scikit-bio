# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import main, TestCase

import numpy as np
import numpy.testing as npt

from skbio.diversity.beta.tests.test_unifrac_base import StatsTests
from skbio.diversity.beta import unweighted_unifrac_fast, weighted_unifrac_fast
from skbio.diversity.beta._fast_unifrac import (_boundary_case, make_pdist,
                                                unifrac, w_unifrac,
                                                _branch_correct)


class FastUnifracTests(StatsTests, TestCase):
    _method = {'unweighted_unifrac': unweighted_unifrac_fast,
               'weighted_unifrac': weighted_unifrac_fast}

    def test_make_pdist_unweighted(self):
        f, counts, length = make_pdist(self.b1, np.asarray(self.oids1),
                                       self.t1, metric=unifrac)
        exp = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        obs = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                obs[i, j] = f(counts[i], counts[j])
                exp[i, j] = self.unweighted_unifrac(self.b1[i], self.b1[j],
                                                    self.oids1, self.t1)
        npt.assert_almost_equal(obs, exp)

    def test_make_pdist_weighted(self):
        f, counts, length = make_pdist(self.b1, np.asarray(self.oids1),
                                       self.t1, metric=w_unifrac)
        exp = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        obs = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                obs[i, j] = f(counts[i], counts[j])
                exp[i, j] = self.weighted_unifrac(self.b1[i], self.b1[j],
                                                  self.oids1, self.t1)
        npt.assert_almost_equal(obs, exp)

    def test_make_pdist_weighted_normalized(self):
        f, counts, length = make_pdist(self.b1, np.asarray(self.oids1),
                                       self.t1, metric=w_unifrac,
                                       normalized=True)
        exp = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        obs = np.zeros((len(self.b1), len(self.b1)), dtype=float)
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                obs[i, j] = f(counts[i], counts[j])
                exp[i, j] = self.weighted_unifrac(self.b1[i], self.b1[j],
                                                  self.oids1, self.t1,
                                                  normalized=True)
        npt.assert_almost_equal(obs, exp)

    def test_boundary_case(self):
        self.assertEqual(_boundary_case(100, 1000), None)
        self.assertEqual(_boundary_case(100, 1000, normalized=True), None)
        self.assertEqual(_boundary_case(100, 1000, unweighted=False), None)
        self.assertEqual(_boundary_case(100, 1000, normalized=True,
                                        unweighted=False), None)

        self.assertEqual(_boundary_case(0, 1), 1.0)
        self.assertEqual(_boundary_case(1, 0), 1.0)
        self.assertEqual(_boundary_case(1, 0, normalized=True), 1.0)
        self.assertEqual(_boundary_case(1, 0, normalized=True,
                                        unweighted=False), 1.0)
        self.assertEqual(_boundary_case(1, 0, normalized=False,
                                        unweighted=False), None)

        self.assertEqual(_boundary_case(0, 0), 0.0)
        self.assertEqual(_boundary_case(0, 0, normalized=True), 0.0)
        self.assertEqual(_boundary_case(0, 0, normalized=True,
                                        unweighted=False), 0.0)
        self.assertEqual(_boundary_case(0, 0, normalized=False,
                                        unweighted=False), 0.0)

    def test_branch_correct(self):
        # for ((a:1, b:2)c:3,(d:4,e:5)f:6)root;"
        tip_ds = np.array([4, 5, 10, 11, 0, 0, 0])[:, np.newaxis]
        u_counts = np.array([1, 1, 0, 0, 2, 0, 2])
        v_counts = np.array([0, 2, 1, 0, 2, 1, 3])
        u_sum = 2  # counts at the tips
        v_sum = 3
        exp = np.array([2.0,
                        5.0 * (.5 + (2.0/3.0)),
                        10.0 * (1.0 / 3.0),
                        0.0]).sum()
        obs = _branch_correct(tip_ds, u_counts, v_counts, u_sum, v_sum)
        self.assertEqual(obs, exp)

    def test_unifrac(self):
        # adapted from PyCogent unit tests
        m = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, 0],
                      [0, 1, 1], [1, 1, 1], [0, 1, 1], [1, 1, 1]])
        # lengths from ((a:1,b:2):4,(c:3,(d:1,e:1):2):3)
        bl = np.array([1, 2, 1, 1, 3, 2, 4, 3, 0], dtype=float)
        self.assertEqual(unifrac(bl, m[:,0], m[:,1]), 10/16.0)
        self.assertEqual(unifrac(bl, m[:,0], m[:,2]), 8/13.0)
        self.assertEqual(unifrac(bl, m[:,1], m[:,2]), 8/17.0)

    def test_w_unifrac(self):
        # lengths from ((a:1,b:2):4,(c:3,(d:1,e:1):2):3)
        bl = np.array([1, 2, 1, 1, 3, 2, 4, 3, 0], dtype=float)

        # adapted from PyCogent unit tests
        m = np.array([[1, 0, 1],  # a
                      [1, 1, 0],  # b
                      [0, 1, 0],  # d
                      [0, 0, 1],  # e
                      [0, 1, 0],  # c
                      [0, 1, 1],  # parent of (d, e)
                      [2, 1, 1],  # parent of a, b
                      [0, 2, 1],  # parent of c (d, e)
                      [2, 3, 2]]) # root

        # sum just the counts at the tips
        m0s = m[:5, 0].sum()
        m1s = m[:5, 1].sum()
        m2s = m[:5, 2].sum()

        # scores computed by educational implementation
        self.assertAlmostEqual(w_unifrac(bl, m[:, 0], m[:, 1], m0s, m1s), 7.5)
        self.assertAlmostEqual(w_unifrac(bl, m[:, 0], m[:, 2], m0s, m2s), 6.0)
        self.assertAlmostEqual(w_unifrac(bl, m[:, 1], m[:, 2], m1s, m2s), 4.5)


if __name__ == '__main__':
    main()
