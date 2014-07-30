from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from random import seed

import numpy as np
from numpy.testing import assert_allclose
from nose.tools import assert_almost_equal, assert_raises, assert_equal

from skbio.math.stats.evolve import hommola_cospeciation
from skbio.math.stats.evolve.hommola import _get_dist

def test_hommola_cospeciation_sig():
  
  seed(1)

  hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
                8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
  pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
                8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
  matrix = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
                 0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])

  obs_p, obs_r, obs_perm_stats = hommola_cospeciation(hdist, pdist, matrix, 9)
  exp_p = .1
  exp_r = 0.83170965463247915
  exp_perm_stats = np.array([-0.09826424,  0.63310285, -0.17556286,  0.73840484,  0.,-0.02529993,  0.47279242,  0.10184894,  0.42742399])
  assert_almost_equal(obs_p, exp_p)
  assert_almost_equal(obs_r, exp_r)

  #print(obs_perm_stats)
  #print(exp_perm_stats)

  assert_allclose(obs_perm_stats, exp_perm_stats)

def test_hommola_cospeciation_no_sig():

  seed(1)

  hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
                8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
  pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
                8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
  # this matrix was picked because it will generate an r value that's less than
  # a standard deviation away from the mean of the normal distribution of
  # r vals
  randomized_matrix = np.array(
      [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])

  obs_p, obs_r, obs_perm_stats = hommola_cospeciation(hdist, pdist, randomized_matrix, 9)
  exp_p = .5
  exp_r = -0.013679391379114569
  exp_perm_stats = np.array([-0.15686395, -0.14836061,  0.27165268, -0.21887026, -0.37597226,
        0.28203804, -0.42305706,  0.56395839,  0.25098232])
  assert_almost_equal(obs_p, exp_p)
  assert_almost_equal(obs_r, exp_r)
  assert_allclose(obs_perm_stats, exp_perm_stats)

def test_get_dist():
    labels = [0, 1, 1, 2, 3]
    dists = np.array([[0, 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7], [3, 4, 7, 0]])
    index = [2, 3, 1, 0]

    expected_vec = [7, 7, 5, 6, 0, 4, 3, 4, 3, 2]
    actual_vec = _get_dist(labels, dists, index)

    assert_equal(actual_vec, expected_vec)

if __name__ == '__main__':
    import nose
    nose.runmodule()

