from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from numpy.testing import assert_allclose
from nose.tools import assert_almost_equal, assert_equal

from skbio.math.stats.evolve import hommola_cospeciation
from skbio.math.stats.evolve.hommola import _get_dist, _gen_lists


def test_hommola_cospeciation_sig():

    np.random.seed(1)

    hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
        8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
    pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
        8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
    matrix = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
        0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])

    obs_p, obs_r, obs_perm_stats = hommola_cospeciation(
        hdist, pdist, matrix, 9)
    exp_p = .1
    exp_r = 0.83170965463247915
    exp_perm_stats = np.array([-0.14928122, 0.26299538, -0.21125858,
                               0.24143838, 0.61557855, -0.24258293,
                               0.09885203, 0.02858, 0.42742399])
    assert_almost_equal(obs_p, exp_p)
    assert_almost_equal(obs_r, exp_r)

    assert_allclose(obs_perm_stats, exp_perm_stats)


def test_hommola_cospeciation_no_sig():

    np.random.seed(1)

    hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
        8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
    pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
        8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
    # this matrix was picked because it will generate an r value that's less
    # than a standard deviation away from the mean of the normal distribution
    # of r vals
    randomized_matrix = np.array(
        [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0], [1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1]])

    obs_p, obs_r, obs_perm_stats = hommola_cospeciation(
        hdist, pdist, randomized_matrix, 9)
    exp_p = .6
    exp_r = -0.013679391379114569
    exp_perm_stats = np.array([-0.22216543, -0.14836061, -0.04434843,
                               0.1478281, -0.29105645, 0.56395839, 0.47304992,
                               0.79125657, 0.06804138])
    assert_almost_equal(obs_p, exp_p)
    assert_almost_equal(obs_r, exp_r)
    assert_allclose(obs_perm_stats, exp_perm_stats)


def test_get_dist():
    labels = np.array([0, 1, 1, 2, 3])
    k_labels, t_labels = _gen_lists(labels)
    dists = np.array([[0, 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7], [3, 4, 7, 0]])
    index = np.array([2, 3, 1, 0])

    expected_vec = np.array([7, 7, 5, 6, 0, 4, 3, 4, 3, 2])
    actual_vec = _get_dist(k_labels, t_labels, dists, index)

    assert_allclose(actual_vec, expected_vec)

def test_gen_lists():
     matrix = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
        0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])
     exp_pars_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4])
     exp_pars_t_labels = np.array([1, 2, 3, 4, 4, 2, 3, 4, 4, 3, 4, 4, 4, 4, 4])
     exp_hosts_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3])
     exp_hosts_t_labels = np.array([1, 2, 3, 3, 4, 2, 3, 3, 4, 3, 3, 4, 3, 4, 4])
     pars, hosts = np.nonzero(matrix)
     obs_pars_k_labels, obs_pars_t_labels = _gen_lists(pars)
     obs_hosts_k_labels, obs_hosts_t_labels = _gen_lists(hosts)
     assert_allclose(exp_pars_k_labels, obs_pars_k_labels)
     assert_allclose(exp_pars_t_labels, obs_pars_t_labels)
     assert_allclose(exp_hosts_k_labels, obs_hosts_k_labels)
     assert_allclose(exp_hosts_t_labels, obs_hosts_t_labels)

if __name__ == '__main__':
    import nose
    nose.runmodule()
