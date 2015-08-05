# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from six import StringIO

from skbio import TreeNode
from skbio.diversity.alpha._stats import phylogenetic_diversity

class StatsTests(TestCase):

    def setUp(self):
        self.b1 = np.array(
           [[1, 3, 0, 1, 0],
            [0, 2, 0, 4, 4],
            [0, 0, 6, 2, 1],
            [0, 0, 1, 1, 1]])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1,6)]
        self.t1 = TreeNode.read(StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:'
            '1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))

    def test_phylogenetic_diversity_empty(self):
        actual = phylogenetic_diversity(
            [0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)

    def test_phylogenetic_diversity_empty(self):
        actual = phylogenetic_diversity(
            [1, 1, 1, 1, 1], self.oids1, self.t1)
        expected = sum([n.length for n in self.t1.traverse()
                        if n.length is not None])
        self.assertAlmostEqual(actual, expected)

        actual = phylogenetic_diversity(
            [1, 2, 3, 4, 5], self.oids1, self.t1)
        expected = sum([n.length for n in self.t1.traverse()
                        if n.length is not None])
        self.assertAlmostEqual(actual, expected)

    def test_phylogenetic_diversity(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # phylogenetic diversity implementation
        actual = phylogenetic_diversity(
            self.b1[0], self.oids1, self.t1)
        expected = 4.5
        self.assertAlmostEqual(actual, expected)
        actual = phylogenetic_diversity(
            self.b1[1], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)
        actual = phylogenetic_diversity(
            self.b1[2], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)
        actual = phylogenetic_diversity(
            self.b1[3], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)

if __name__ == "__main__":
    main()
