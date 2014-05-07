#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
from nose.tools import (assert_equal, assert_almost_equal, assert_raises,
                        assert_true)

from skbio.maths.diversity.alpha.lladser import (
    _expand_counts, _lladser_point_estimates, _get_interval_for_r_new_species,
    _lladser_ci_series, _lladser_ci_from_r)


def test_expand_counts():
    arr = np.array([2, 0, 1, 2])
    npt.assert_array_equal(_expand_counts(arr), np.array([0, 0, 2, 3, 3]))


def test_lladser_point_estimates():
    s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5, 3]
    r = 3
    observed = list(_lladser_point_estimates(s, r))
    assert_equal(len(observed), 3)

    for k in (range(3)):
        x = observed[k]
        t = x[2]
        assert_equal(x[0], (r - 1) / t)

    # Estimator has variance of (1-p)^2/(r-2),
    # which for r=7 and p=0.5 is 0.05
    seq = "WBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBW"
    reps = 1000
    sum = 0
    for i in range(reps):
        p, _, _ = list(_lladser_point_estimates(seq, r=7))[0]
        sum += p
    assert_true(0.45 < sum / reps and sum / reps < 0.55)


def test_get_interval_for_r_new_species():
    s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5]
    expected = [(3, set([5]), 4, 0),
                (4, set([5, 1]), 6, 1),
                (4, set([5, 1, 2]), 9, 4)]
    for x, y in zip(_get_interval_for_r_new_species(s, 2), expected):
        assert_equal(x, y)

    s = [5, 5, 5, 5, 5]
    # never saw new one
    assert_equal(list(_get_interval_for_r_new_species(s, 2)), [])


def test_lladser_ci_series_exact():
    # have seen RWB
    urn_1 = 'RWBWWBWRRWRYWRPPZ'
    results = list(_lladser_ci_series(urn_1, r=4))
    assert_equal(len(results), 3)


def test_lladser_ci_series_random():
    seq = "WBWBWBWBWBWB"
    observations = []
    alpha = 0.95
    reps = 1000
    for i in range(reps):
        obs = list(_lladser_ci_series(seq, r=4, alpha=alpha))[0]
        observations.append(obs)
    tps = list(filter(lambda a_b: a_b[0] < 0.5 and 0.5 < a_b[1], observations))
    assert_true(len(tps) >= alpha * reps)  # 100%-95%


def test_lladser_ci_from_r():
    f = 10
    t = 10
    r = 4
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f)
    assert_almost_equal(obs_low, 0.0806026244)
    assert_almost_equal(obs_high, 0.806026244)

    r = 20
    t = 100
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f)
    assert_almost_equal(obs_low, 0.02787923964)
    assert_almost_equal(obs_high, 0.2787923964)

    # make sure we test with each possible alpha
    alpha = 0.99
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
    assert_almost_equal(obs_low, 0.03184536992)
    assert_almost_equal(obs_high, 0.3184536992)

    alpha = 0.9
    r = 3
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
    assert_almost_equal(obs_low, 0.005635941995)
    assert_almost_equal(obs_high, 0.05635941995)

    # test other ci_types
    ci_type = 'ULCU'
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                           ci_type=ci_type)
    assert_almost_equal(obs_low, 0.01095834700)
    assert_almost_equal(obs_high, 0.1095834700)

    alpha = 0.95
    t = 10
    ci_type = 'U'
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                           ci_type=ci_type)
    assert_almost_equal(obs_low, 0)
    assert_almost_equal(obs_high, 0.6295793622)

    ci_type = 'L'
    obs_low, obs_high = _lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha,
                                           ci_type=ci_type)
    assert_almost_equal(obs_low, 0.0817691447)
    assert_almost_equal(obs_high, 1)

    # Requesting CI for not precomputed values raises error
    with assert_raises(ValueError):
        _lladser_ci_from_r(r=500, t=t, f=f, alpha=alpha, ci_type=ci_type)


if __name__ == '__main__':
    import nose
    nose.runmodule()
