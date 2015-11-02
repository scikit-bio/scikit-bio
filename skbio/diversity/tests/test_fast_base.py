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

from skbio import TreeNode
from skbio.diversity._fast_base import index_tree, counts_and_index


class FastBaseTests(TestCase):
    def setUp(self):
        self.t = TreeNode.read(StringIO(u"((a:1, b:2)c:3)root;"))

    def test_index_tree(self):
        indexed = index_tree(self.t)
        npt.assert_equal(indexed['length'], np.array([1, 2, 3, 0],
                                                     dtype=float))

    def test_counts_and_index(self):
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        count_array, indexed = counts_and_index(counts, np.array(['a', 'b']),
                                                self.t, None)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts)

    def test_counts_and_length_with_index(self):
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        indexed = index_tree(self.t)
        count_array, indexed = counts_and_index(counts, np.array(['a', 'b']),
                                                self.t, indexed)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts)


if __name__ == '__main__':
    main()
