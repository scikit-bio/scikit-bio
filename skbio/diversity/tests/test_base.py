# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
from io import StringIO

import numpy as np
import numpy.testing as npt

from skbio.diversity._base import (
    _validate_counts_vector, _validate_counts_vectors,
    _validate_otu_ids_and_tree)
from skbio import TreeNode
from skbio.tree import DuplicateNodeError, MissingNodeError


class BaseTests(TestCase):

    def test_validate_counts_vector(self):
        # python list
        obs = _validate_counts_vector([0, 2, 1, 3])
        npt.assert_array_equal(obs, np.array([0, 2, 1, 3]))
        self.assertEqual(obs.dtype, int)

        # numpy array (no copy made)
        data = np.array([0, 2, 1, 3])
        obs = _validate_counts_vector(data)
        npt.assert_array_equal(obs, data)
        self.assertEqual(obs.dtype, int)
        self.assertTrue(obs is data)

        # single element
        obs = _validate_counts_vector([42])
        npt.assert_array_equal(obs, np.array([42]))
        self.assertEqual(obs.dtype, int)
        self.assertEqual(obs.shape, (1,))

        # suppress casting to int
        obs = _validate_counts_vector([42.2, 42.1, 0], suppress_cast=True)
        npt.assert_array_equal(obs, np.array([42.2, 42.1, 0]))
        self.assertEqual(obs.dtype, float)

        # all zeros
        obs = _validate_counts_vector([0, 0, 0])
        npt.assert_array_equal(obs, np.array([0, 0, 0]))
        self.assertEqual(obs.dtype, int)

        # all zeros (single value)
        obs = _validate_counts_vector([0])
        npt.assert_array_equal(obs, np.array([0]))
        self.assertEqual(obs.dtype, int)

    def test_validate_counts_vector_invalid_input(self):
        # wrong dtype
        with self.assertRaises(TypeError):
            _validate_counts_vector([0, 2, 1.2, 3])

        # wrong number of dimensions (2-D)
        with self.assertRaises(ValueError):
            _validate_counts_vector([[0, 2, 1, 3], [4, 5, 6, 7]])

        # wrong number of dimensions (scalar)
        with self.assertRaises(ValueError):
            _validate_counts_vector(1)

        # negative values
        with self.assertRaises(ValueError):
            _validate_counts_vector([0, 0, 2, -1, 3])

    def test_validate_counts_vectors(self):
        # basic valid input (n=2)
        obs_u, obs_v = _validate_counts_vectors([0, 1, 1, 0, 2],
                                                [0, 0, 2, 1, 3])
        npt.assert_array_equal(obs_u, np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(obs_v, np.array([0, 0, 2, 1, 3]))

        # basic valid input (n=3)
        actual = _validate_counts_vectors([0, 1, 1, 0, 2],
                                          [0, 0, 2, 1, 3],
                                          [1, 1, 1, 1, 1])
        npt.assert_array_equal(actual[0], np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(actual[1], np.array([0, 0, 2, 1, 3]))
        npt.assert_array_equal(actual[2], np.array([1, 1, 1, 1, 1]))

        # empty counts vectors
        obs_u, obs_v = _validate_counts_vectors(np.array([], dtype=int),
                                                np.array([], dtype=int))
        npt.assert_array_equal(obs_u, np.array([]))
        npt.assert_array_equal(obs_v, np.array([]))

    def test_validate_counts_vectors_suppress_cast(self):
        # suppress_cast is passed through to _validate_counts_vector
        obs_u, obs_v = _validate_counts_vectors(
            [42.2, 42.1, 0], [42.2, 42.1, 1.0], suppress_cast=True)
        npt.assert_array_equal(obs_u, np.array([42.2, 42.1, 0]))
        npt.assert_array_equal(obs_v, np.array([42.2, 42.1, 1.0]))
        self.assertEqual(obs_u.dtype, float)
        self.assertEqual(obs_v.dtype, float)
        with self.assertRaises(TypeError):
            _validate_counts_vectors([0.0], [1], suppress_cast=False)

    def test_validate_counts_vectors_invalid_input(self):
        # checks that are caught by the calls to _validate_counts_vector
        with self.assertRaises(ValueError):
            _validate_counts_vectors([0, 1, 1, 0, 2], [0, 0, 2, -1, 3])
        with self.assertRaises(ValueError):
            _validate_counts_vectors([0, 0, 2, -1, 3], [0, 1, 1, 0, 2])

        # len of vectors not equal
        u_counts = [1, 2]
        v_counts = [1, 1, 1]
        self.assertRaises(ValueError, _validate_counts_vectors, u_counts,
                          v_counts)
        u_counts = [1, 2, 3]
        v_counts = [1, 1]
        self.assertRaises(ValueError, _validate_counts_vectors, u_counts,
                          v_counts)

    def test_validate_otu_ids_and_tree(self):
        # basic valid input
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all tips observed
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # no tips observed
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = []
        otu_ids = []
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all counts zero
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [0, 0, 0, 0, 0]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # single node tree
        t = TreeNode.read(StringIO(u'root;'))
        counts = []
        otu_ids = []
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

    def test_validate_otu_ids_and_tree_invalid_input(self):
        # tree has duplicated tip ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, _validate_otu_ids_and_tree,
                          counts, otu_ids, t)

        # unrooted tree as input
        t = TreeNode.read(StringIO(u'((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   u'OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO(u'((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO(u'(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.25,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU32']
        self.assertRaises(MissingNodeError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)


if __name__ == '__main__':
    main()
