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
from skbio.tree._compare import rf_dists, wrf_dists, path_dists


class CompareTests(TestCase):

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

        # with unique taxa
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

        # duplicate IDs
        with self.assertRaises(ValueError):
            rf_dists(trees, list("xxxx"))

        # ID number doesn't match
        with self.assertRaises(ValueError):
            rf_dists(trees, list("abc"))

        # no shared taxon
        nwks = ["((a,b),(c,d));",
                "((e,f),(g,h));"]
        trees = [TreeNode.read([x]) for x in nwks]
        obs = rf_dists(trees).data
        self.assertFalse(np.any(obs))

        # only one tree
        obs = rf_dists([TreeNode.read(["((a,b),(c,d));"])]).data
        npt.assert_array_equal(obs, np.zeros((1, 1)))

        # no tree
        with self.assertRaises(ValueError):
            rf_dists([])

    def test_wrf_dists(self):
        """Calculate weighted Robinson-Foulds distances among trees."""
        # all trees have identical taxa
        nwks = ["((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);",
                "((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);",
                "((a:1,c:6):2,(b:3,(d:2,e:3):1):2,f:7);"]
        trees = [TreeNode.read([x]) for x in nwks]

        obs = wrf_dists(trees)
        self.assertSequenceEqual(obs.ids, list("012"))
        exp = np.array([[0, 16, 15],
                        [16, 0, 27],
                        [15, 27, 0]], dtype=float)
        npt.assert_array_almost_equal(obs.data, exp)

        # with IDs
        ids = list("ABC")
        obs = wrf_dists(trees, ids)
        self.assertSequenceEqual(obs.ids, ids)
        npt.assert_array_almost_equal(obs.data, exp)

        # rooted
        obs = wrf_dists(trees, rooted=True).data
        exp = np.array([[0, 18, 15],
                        [18, 0, 27],
                        [15, 27, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # Euclidean distances
        obs = wrf_dists(trees, metric="euclidean").data
        exp = np.array([[0, 6.164414, 5],
                        [6.164414, 0, 9.21954446],
                        [5, 9.21954446, 0]])
        npt.assert_array_almost_equal(obs, exp)

        # unit correlation distances
        obs = wrf_dists(trees, metric="unitcorr").data
        exp = np.array([[0, 0.15822237, 0.13211014],
                        [0.15822237, 0, 0.32552865],
                        [0.13211014, 0.32552865, 0]])
        npt.assert_array_almost_equal(obs, exp)

        # invalide distance metric
        with self.assertRaises(ValueError):
            wrf_dists(trees, metric=100)
        with self.assertRaises(AttributeError):
            wrf_dists(trees, metric="hello")

        # with unique taxa
        nwks = ["((a:1,b:2):1,(c:3,d:4):2);",
                "((a:1,(b:2,d:4):1):2,e:5);",
                "(b:2,((c:3,d:4):1,e:5):2);"]
        trees = [TreeNode.read([x]) for x in nwks]
        obs = wrf_dists(trees).data
        exp = np.array([[0, 3, 0],
                        [3, 0, 3],
                        [0, 3, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # pairwise sharing
        obs = wrf_dists(trees, shared_by_all=False).data
        exp = np.array([[0, 4, 0],
                        [4, 0, 6],
                        [0, 6, 0]], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

    def test_path_dists(self):
        """Calculate path-length distances among trees."""
        # all trees have identical taxa
        nwks = ["((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);",
                "((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);",
                "((a:1,c:6):2,(b:3,(d:2,e:3):1):2,f:7);"]
        trees = [TreeNode.read([x]) for x in nwks]

        obs = path_dists(trees)
        self.assertSequenceEqual(obs.ids, list("012"))
        exp = np.array([13.7113092, 11.87434209, 19.5192213])
        npt.assert_array_almost_equal(obs.condensed_form(), exp)

        # with IDs
        ids = list("ABC")
        obs = path_dists(trees, ids)
        self.assertSequenceEqual(obs.ids, ids)
        npt.assert_array_almost_equal(obs.condensed_form(), exp)

        # city block distances
        obs = path_dists(trees, metric="cityblock").condensed_form()
        exp = np.array([48, 37, 61])
        npt.assert_array_almost_equal(obs, exp)

        # unit correlation distances
        obs = path_dists(trees, metric="unitcorr").condensed_form()
        exp = np.array([0.14130903, 0.2805951, 0.48826796])
        npt.assert_array_almost_equal(obs, exp)

        # random sampling
        obs = path_dists(trees, sample=3, shuffler=list.sort).condensed_form()
        exp = np.array([4.24264069, 7.87400787, 9.2736185])
        npt.assert_array_almost_equal(obs, exp)

        # too many samples
        with self.assertRaises(ValueError):
            path_dists(trees, sample=10)

        # with unique taxa
        nwks = ["((a:1,b:2):1,(c:3,d:4):2,e:5);",
                "((a:1,(b:2,d:4):1):2,e:5,f:6);",
                "(a:1,b:2,((c:3,d:4):1,e:5):2);"]
        trees = [TreeNode.read([x]) for x in nwks]
        obs = path_dists(trees).condensed_form()
        exp = np.array([4.47213595, 1.73205081, 4.35889894])
        npt.assert_array_almost_equal(obs, exp)

        # pairwise sharing
        obs = path_dists(trees, shared_by_all=False).condensed_form()
        exp = np.array([4.47213595, 2.0, 4.35889894])
        npt.assert_array_almost_equal(obs, exp)


if __name__ == "__main__":
    main()
