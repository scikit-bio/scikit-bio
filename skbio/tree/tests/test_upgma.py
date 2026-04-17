import io
from unittest import TestCase, main

import numpy as np
from skbio import DistanceMatrix, TreeNode

from skbio.tree._upgma import (
    upgma
)


class UpgmaTests(TestCase):

    def setUp(self):
        data = [[0,  5,  9,  9,  8],
                [5,  0, 10, 10,  9],
                [9, 10,  0,  8,  7],
                [9, 10,  8,  0,  3],
                [8,  9,  7,  3,  0]]
        ids = list('abcde')
        self.dm = DistanceMatrix(data, ids)
        self.dm_condensed = DistanceMatrix(self.dm, condensed=True)
        # Note that upgma and wpgma have different edge
        # estimations, even with small examples.
        self.expected_str_upgma = (
            "((a:2.5,b:2.5):2.083333333333333,(c:3.75,(d:1.5,e:1.5):2.25):0.833333333333333);")
        self.expected_str_wpgma = (
            "((a:2.5,b:2.5):2.125,(c:3.75,(d:1.5,e:1.5):2.25):0.875);")
        self.expected_TreeNode_upgma = TreeNode.read(
            io.StringIO(self.expected_str_upgma))
        self.expected_TreeNode_wpgma = TreeNode.read(
            io.StringIO(self.expected_str_wpgma))

    def test_upgma_dm(self):
        actual_TreeNode = upgma(self.dm)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected_TreeNode_upgma), 0.0, places=10)

    def test_upgma_dm_condensed(self):
        actual_TreeNode = upgma(self.dm_condensed)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected_TreeNode_upgma), 0.0, places=10)

    def test_wpgma_dm(self):
        actual_TreeNode = upgma(self.dm, weighted=True)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected_TreeNode_wpgma), 0.0, places=10)

    def test_upgma_fixed_point_scale(self):
        scale = 0.5
        float_data = self.dm.data.astype(np.float64)
        fixed_data = np.rint(float_data / scale).astype(np.int16)

        dm_float = DistanceMatrix(float_data, self.dm.ids)
        dm_fixed = DistanceMatrix(fixed_data, self.dm.ids, scale=scale)
        dm_fixed_condensed = DistanceMatrix(
            fixed_data, self.dm.ids, condensed=True, scale=scale
        )

        exp = upgma(dm_float)
        obs = upgma(dm_fixed)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0, places=10)
        obs = upgma(dm_fixed_condensed)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0, places=10)

        exp_w = upgma(dm_float, weighted=True)
        obs_w = upgma(dm_fixed, weighted=True)
        self.assertAlmostEqual(obs_w.compare_cophenet(exp_w), 0.0, places=10)
