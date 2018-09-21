# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.stats import subsample_counts
from skbio.diversity.alpha import lladser_pe, lladser_ci
from skbio.diversity.alpha._lladser import (
    _expand_counts, _lladser_point_estimates,
    _get_interval_for_r_new_otus, _lladser_ci_series, _lladser_ci_from_r)


def create_fake_observation():
    """Create a subsample with defined property"""

    # Create a subsample of a larger sample such that we can compute
    # the expected probability of the unseen portion.
    # This is used in the tests of lladser_pe and lladser_ci
    counts = np.ones(1001, dtype='int64')
    counts[0] = 9000
    total = counts.sum()

    fake_obs = subsample_counts(counts, 1000)
    exp_p = 1 - sum([x/total for (x, y) in zip(counts, fake_obs) if y > 0])

    return fake_obs, exp_p


class LladserTests(unittest.TestCase):
    def test_lladser_pe(self):
        """lladser_pe returns point estimates within the expected variance"""

        obs = lladser_pe([3], r=4)
        self.assertTrue(np.isnan(obs))

        np.random.seed(123456789)
        fake_obs, exp_p = create_fake_observation()
        reps = 100
        sum = 0
        for i in range(reps):
            sum += lladser_pe(fake_obs, r=30)
        obs = sum / reps

        # Estimator has variance of (1-p)^2/(r-2),
        # which for r=30 and p~=0.9 is 0.0289
        npt.assert_almost_equal(obs, exp_p, decimal=2)

    def test_lladser_ci_nan(self):
        """lladser_ci returns nan if sample is too short to make an estimate"""
        obs = lladser_ci([3], r=4)
        self.assertTrue(len(obs) == 2 and
                        np.isnan(obs[0]) and
                        np.isnan(obs[1]))

    def test_lladser_ci(self):
        """lladser_ci estimate using defaults contains p with 95% prob"""
        np.random.seed(12345678)
        reps = 100
        sum = 0
        for i in range(reps):
            fake_obs, exp_p = create_fake_observation()
            (low, high) = lladser_ci(fake_obs, r=10)
            if (low <= exp_p <= high):
                sum += 1

        self.assertTrue(sum/reps >= 0.95)

    def test_lladser_ci_f3(self):
        """lladser_ci estimate using f=3 contains p with 95% prob"""

        # Test different values of f=3 and r=14, which lie exactly on the
        # 95% interval line. For 100 reps using simple cumulative binomial
        # probs we expect to have more than 5 misses of the interval in 38%
        # of all test runs. To make this test pass reliable we thus have to
        # set a defined seed
        np.random.seed(12345678)
        reps = 100
        sum = 0
        for i in range(reps):
            # re-create the obs for every estimate, such that they are truly
            # independent events
            fake_obs, exp_p = create_fake_observation()
            (low, high) = lladser_ci(fake_obs, r=14, f=3)
            if (low <= exp_p <= high):
                sum += 1

        self.assertTrue(sum/reps >= 0.95)

    def test_expand_counts(self):
        arr = np.array([2, 0, 1, 2])
        npt.assert_array_equal(_expand_counts(arr), np.array([0, 0, 2, 3, 3]))

    def test_lladser_point_estimates(self):
        s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5, 3]
        r = 3
        observed = list(_lladser_point_estimates(s, r))
        self.assertEqual(len(observed), 3)

        for k in range(3):
            x = observed[k]
            t = x[2]
            self.assertEqual(x[0], (r - 1) / t)

        # Estimator has variance of (1-p)^2/(r-2),
        # which for r=7 and p=0.5 is 0.05
        seq = "WBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBW"
        reps = 1000
        sum = 0
        for i in range(reps):
            p, _, _ = list(_lladser_point_estimates(seq, r=7))[0]
            sum += p
        self.assertTrue(0.45 < sum / reps and sum / reps < 0.55)

    def test_lladser_point_estimates_invalid_r(self):
        with self.assertRaises(ValueError):
            list(_lladser_point_estimates([5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5, 3],
                                          2))

    def test_get_interval_for_r_new_otus(self):
        s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5]
        expected = [(3, set([5]), 4, 0),
                    (4, set([5, 1]), 6, 1),
                    (4, set([5, 1, 2]), 9, 4)]
        for x, y in zip(_get_interval_for_r_new_otus(s, 2), expected):
            self.assertEqual(x, y)

        s = [5, 5, 5, 5, 5]
        # never saw new one
        self.assertEqual(list(_get_interval_for_r_new_otus(s, 2)), [])

    def test_lladser_ci_series_exact(self):
        # have seen RWB
        urn_1 = 'RWBWWBWRRWRYWRPPZ'
        results = list(_lladser_ci_series(urn_1, r=4))
        self.assertEqual(len(results), 3)

    def test_lladser_ci_series_random(self):
        seq = "WBWBWBWBWBWB"
        observations = []
        alpha = 0.95
        reps = 1000
        for i in range(reps):
            obs = list(_lladser_ci_series(seq, r=4, alpha=alpha))[0]
            observations.append(obs)
        tps = list(filter(lambda a_b: a_b[0] < 0.5 and 0.5 < a_b[1],
                          observations))
        self.assertTrue(len(tps) >= alpha * reps)  # 100%-95%

    def test_lladser_ci_from_r(self):
        f = 10
        t = 10
        r = 4
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f)
        npt.assert_almost_equal(obs_low, 0.0806026244)
        npt.assert_almost_equal(obs_high, 0.806026244)

        r = 20
        t = 100
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f)
        npt.assert_almost_equal(obs_low, 0.02787923964)
        npt.assert_almost_equal(obs_high, 0.2787923964)

        # make sure we test with each possible alpha
        alpha = 0.99
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        npt.assert_almost_equal(obs_low, 0.03184536992)
        npt.assert_almost_equal(obs_high, 0.3184536992)

        alpha = 0.9
        r = 3
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        npt.assert_almost_equal(obs_low, 0.005635941995)
        npt.assert_almost_equal(obs_high, 0.05635941995)

        # test other ci_types
        ci_type = 'ULCU'
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                               ci_type=ci_type)
        npt.assert_almost_equal(obs_low, 0.01095834700)
        npt.assert_almost_equal(obs_high, 0.1095834700)

        alpha = 0.95
        t = 10
        ci_type = 'U'
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                               ci_type=ci_type)
        npt.assert_almost_equal(obs_low, 0)
        npt.assert_almost_equal(obs_high, 0.6295793622)

        ci_type = 'L'
        obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                               ci_type=ci_type)
        npt.assert_almost_equal(obs_low, 0.0817691447)
        npt.assert_almost_equal(obs_high, 1)

    def test_lladser_ci_from_r_invalid_input(self):
        # unsupported alpha for ci_type='U'
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=3, t=10, f=10, alpha=0.90, ci_type='U')

        # unsupported r for ci_type='U'
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=42, t=10, f=10, alpha=0.95, ci_type='U')

        # unsupported alpha for ci_type='L'
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=3, t=10, f=10, alpha=0.90, ci_type='L')

        # unsupported r for ci_type='L'
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=50, t=10, f=10, alpha=0.95, ci_type='L')

        # unknown ci_type
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=4, t=10, f=10, alpha=0.95, ci_type='brofist')

        # requesting CI for not precomputed values
        with self.assertRaises(ValueError):
            _lladser_ci_from_r(r=500, t=10, f=10)


if __name__ == '__main__':
    unittest.main()
