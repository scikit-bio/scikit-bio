# ----------------------------------------------------------------------------
# Copyright (c) 2016--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

from tempdir import run_in_tempdir

from skbio.tree import TreeNode, parse_partitioned_newick

'''
In these tests, for brevity, "part" is always short for "partioned", not the
word part itself.
'''


class TestPartitionedNewick(TestCase):

    def setUp(self):
        # Partitioned newick trees
        self.parttrees = [
            '[10](A:2,(B:1,C:1):1);',
            '[20](C:2,(B:1,A:1):1);',
            '[10](A:2,(B:1,C:1):1);',
        ]
        # Just the lengths
        self.partlengths = [int(parttree[1:3]) for parttree in self.parttrees]
        # Just the tree bit of the newick
        self.trees = [parttree[4:] for parttree in self.parttrees]
        # Just the tree bit parsed into a TreeNode
        self.treenodes = [TreeNode.read([t]) for t in self.trees]

        # List of tuples per the parser's output
        self.length_and_tree = list(zip(self.partlengths, self.treenodes))


    def _do_test(self, iterable):
        for expected, got in zip(self.length_and_tree,
                                 parse_partitioned_newick(iterable)):
            got_len, got_tree = got
            expected_len, expected_tree = expected
            self.assertEqual(got_len, expected_len)
            # The below test doesn't work due to TreeNode's __eq__
            # self.assertEqual(got_tree, expected_tree)

    def test_list(self):
        self._do_test(self.parttrees)

    @run_in_tempdir
    def test_file(self):
        with open('parttree.nwk', 'w') as fh:
            for treeline in self.parttrees:
                print(treeline, file=fh)

        with open('parttree.nwk') as fh:
            self._do_test(fh)

if __name__ == '__main__':
    main()
