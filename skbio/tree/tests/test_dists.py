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

from skbio import TreeNode
from skbio.tree._dists import rf_dists


class DistsTests(TestCase):

    def test_rf_dists(self):
        """Calculate Robinson-Foulds distances among trees."""
        # all trees have identical taxa
        nwks = ["(((a,b),c),d,e);",
                "((a,(b,c)),d,e);",
                "((a,b),(c,d),e);",
                "(a,b,(c,(d,e)));"]
        trees = [TreeNode.read([x]) for x in nwks]

        obs = rf_dists(trees)
        self.assertSequenceEqual(obs.ids, list("0123"))
        exp = np.array([[0, 2, 2, 0],
                        [2, 0, 4, 2],
                        [2, 4, 0, 2],
                        [0, 2, 2, 0]], dtype=float)
        npt.assert_array_almost_equal(obs.data, exp)

        # with IDs
        ids = list("ABCD")
        obs = rf_dists(trees, ids)
        self.assertSequenceEqual(obs.ids, ids)
        npt.assert_array_almost_equal(obs.data, exp)

        # rooted
        obs = rf_dists(trees, rooted=True).data
        exp = np.array([[0, 2, 2, 4],
                        [2, 0, 4, 4],
                        [2, 4, 0, 4],
                        [4, 4, 4, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # proportion
        obs = rf_dists(trees, proportion=True).data
        exp = np.array([[0, .5, .5, 0],
                        [.5, 0, 1, .5],
                        [.5, 1, 0, .5],
                        [0, .5, .5, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # not all taxa are shared
        nwks = ["(((a,b),c),((d,e),f),g);",
                "((a,(b,c)),(d,(e,f)),h);",
                "(((a,b),(c,e)),(f,h),i);",
                "((a,((b,c),f)),g,(h,e));"]
        trees = [TreeNode.read([x]) for x in nwks]
        obs = rf_dists(trees).data
        exp = np.array([[0, 2, 2, 4],
                        [2, 0, 4, 2],
                        [2, 4, 0, 4],
                        [4, 2, 4, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # pairwise sharing
        obs = rf_dists(trees, shared_by_all=False).data
        exp = np.array([[0, 4, 2, 6],
                        [4, 0, 6, 4],
                        [2, 6, 0, 6],
                        [6, 4, 6, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # no shared taxon among trees
        nwks = ["((a,b),(c,d));",
                "((e,f),(g,h));"]
        trees = [TreeNode.read([x]) for x in nwks]
        with self.assertRaises(ValueError):
            rf_dists(trees)


if __name__ == "__main__":
    main()
