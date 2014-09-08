from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix
from skbio.diversity.beta import pw_distances, pw_distances_from_table


class HelperBiomTable(object):
    """An object that looks like a BIOM table, for use in testing

    This allows us to test passing BIOM-like objects, without having to
    depend on the biom-format project (since this would ultimately be a
    circular dependency).
    """

    def __init__(self, data, observation_ids, sample_ids):
        self._data = data.T
        self.observation_ids = observation_ids
        self.sample_ids = sample_ids

    def data(self, sample_id):
        i = self.sample_ids.index(sample_id)
        return self._data[i]


class BaseTests(TestCase):
    def setUp(self):
        self.t1 = [[1, 5],
                   [2, 3],
                   [0, 1]]
        self.ids1 = list('ABC')

        self.t2 = [[23, 64, 14, 0, 0, 3, 1],
                   [0, 3, 35, 42, 0, 12, 1],
                   [0, 5, 5, 0, 40, 40, 0],
                   [44, 35, 9, 0, 1, 0, 0],
                   [0, 2, 8, 0, 35, 45, 1],
                   [0, 0, 25, 35, 0, 19, 0]]
        self.ids2 = list('ABCDEF')

        # In the future, if necessary, it should be possible to just replace
        # HelperBiomTable with Table in the following lines to test with the
        # biom.table.Table object directly (i.e., this constructor
        # interface aligns with the biom.table.Table constructor
        # interface).
        self.table1 = HelperBiomTable(
            np.array(self.t1).T, observation_ids=range(2),
            sample_ids=self.ids1)
        self.table2 = HelperBiomTable(
            np.array(self.t2).T, observation_ids=range(7),
            sample_ids=self.ids2)

    def test_pw_distances_invalid_input(self):
        # number of ids doesn't match the number of samples
        self.assertRaises(ValueError, pw_distances, self.t1, list('AB'),
                          'euclidean')

    def test_pw_distances_euclidean(self):
        actual_dm = pw_distances(self.t1, self.ids1, 'euclidean')
        self.assertEqual(actual_dm.shape, (3, 3))
        npt.assert_almost_equal(actual_dm['A', 'A'], 0.0)
        npt.assert_almost_equal(actual_dm['B', 'B'], 0.0)
        npt.assert_almost_equal(actual_dm['C', 'C'], 0.0)
        npt.assert_almost_equal(actual_dm['A', 'B'], 2.23606798)
        npt.assert_almost_equal(actual_dm['B', 'A'], 2.23606798)
        npt.assert_almost_equal(actual_dm['A', 'C'], 4.12310563)
        npt.assert_almost_equal(actual_dm['C', 'A'], 4.12310563)
        npt.assert_almost_equal(actual_dm['B', 'C'], 2.82842712)
        npt.assert_almost_equal(actual_dm['C', 'B'], 2.82842712)

        actual_dm = pw_distances(self.t2, self.ids2, 'euclidean')
        expected_data = [
            [0., 80.8455317, 84.0297566, 36.3042697, 86.0116271, 78.9176786],
            [80.8455317, 0., 71.0844568, 74.4714710, 69.3397433, 14.422205],
            [84.0297566, 71.0844568, 0., 77.2851861, 8.3066238, 60.7536007],
            [36.3042697, 74.4714710, 77.2851861, 0., 78.7908624, 70.7389567],
            [86.0116271, 69.3397433, 8.3066238, 78.7908624, 0., 58.4807660],
            [78.9176786, 14.422205, 60.7536007, 70.7389567, 58.4807660, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.ids2)
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_pw_distances_braycurtis(self):
        actual_dm = pw_distances(self.t1, self.ids1, 'braycurtis')
        self.assertEqual(actual_dm.shape, (3, 3))
        npt.assert_almost_equal(actual_dm['A', 'A'], 0.0)
        npt.assert_almost_equal(actual_dm['B', 'B'], 0.0)
        npt.assert_almost_equal(actual_dm['C', 'C'], 0.0)
        npt.assert_almost_equal(actual_dm['A', 'B'], 0.27272727)
        npt.assert_almost_equal(actual_dm['B', 'A'], 0.27272727)
        npt.assert_almost_equal(actual_dm['A', 'C'], 0.71428571)
        npt.assert_almost_equal(actual_dm['C', 'A'], 0.71428571)
        npt.assert_almost_equal(actual_dm['B', 'C'], 0.66666667)
        npt.assert_almost_equal(actual_dm['C', 'B'], 0.66666667)

        actual_dm = pw_distances(self.t2, self.ids2, 'braycurtis')
        expected_data = [
            [0., 0.78787879, 0.86666667, 0.30927835, 0.85714286, 0.81521739],
            [0.78787879, 0., 0.78142077, 0.86813187, 0.75, 0.1627907],
            [0.86666667, 0.78142077, 0., 0.87709497, 0.09392265, 0.71597633],
            [0.30927835, 0.86813187, 0.87709497, 0., 0.87777778, 0.89285714],
            [0.85714286, 0.75, 0.09392265, 0.87777778, 0., 0.68235294],
            [0.81521739, 0.1627907, 0.71597633, 0.89285714, 0.68235294, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.ids2)
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_pw_distances_from_table_euclidean(self):
        # results are equal when passed as Table or matrix
        m_dm = pw_distances(self.t1, self.ids1, 'euclidean')
        t_dm = npt.assert_warns(
            UserWarning, pw_distances_from_table, self.table1, 'euclidean')
        for id1 in self.ids1:
            for id2 in self.ids1:
                npt.assert_almost_equal(m_dm[id1, id2], t_dm[id1, id2])

        m_dm = pw_distances(self.t2, self.ids2, 'euclidean')
        t_dm = npt.assert_warns(
            UserWarning, pw_distances_from_table, self.table2, 'euclidean')
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(m_dm[id1, id2], t_dm[id1, id2])

    def test_pw_distances_from_table_braycurtis(self):
        # results are equal when passed as Table or matrix
        m_dm = pw_distances(self.t1, self.ids1, 'braycurtis')
        t_dm = npt.assert_warns(
            UserWarning, pw_distances_from_table, self.table1, 'braycurtis')
        for id1 in self.ids1:
            for id2 in self.ids1:
                npt.assert_almost_equal(m_dm[id1, id2], t_dm[id1, id2])

        m_dm = pw_distances(self.t2, self.ids2, 'braycurtis')
        t_dm = npt.assert_warns(
            UserWarning, pw_distances_from_table, self.table2, 'braycurtis')
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(m_dm[id1, id2], t_dm[id1, id2])
