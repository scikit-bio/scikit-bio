# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np
from skbio import TreeNode

from skbio.tree._util import (
    _ordered, _pair_lca,
    _subtree_root, _parent, _sibling,
    _ancestors, _subtree, _move_subtree,
    _move_node, _array_to_tree)
from skbio.tree._cutils import num_dist_cy


class UtilTests(TestCase):

    def setUp(self):
        self.ordered_array = [3, 4, 1, 11, 25, 26, 12, 5, 13, 14, 6, 2, 0]
        self.tree_array = np.array([[1, 5, 13, 14, 0], [1, 2, 3, 4, 0]])
        self.expected_array = np.array([[1, 5, 27, 28, 0], [1, 2, 3, 4, 0]])
        self.expected_str = ("((b,(c,(d,e))))a;")
        self.expected_TreeNode = TreeNode.read(
                io.StringIO(self.expected_str))

    def test_ordered(self):
        expected_array = np.array([[1, 5, 13, 14, 6, 2, 0], [1, 2, 3, 4, 0, 0, 0]])
        self.assertEqual(_ordered(self.tree_array), expected_array)

    def test_pair_lca(self):
        # Values are node locations
        data_1 = [7, 14]
        data_2 = [12, 13]
        data_3 = [12, 14]
        self.assertEqual(_pair_lca(data_1[0], data_1[1]), 0)
        self.assertEqual(_pair_lca(data_2[0], data_2[1]), 2)
        self.assertEqual(_pair_lca(data_3[0], data_3[1]), 2)
    
    def test_num_dist(self):
        # checking ancestral line
        data_1 = [1, 8]
        # checking non-ancestral line
        data_2 = [1, 2]
        # checking distance to self
        data_3 = [1, 1]
        self.assertEqual(num_dist_cy(data_1[0], data_1[1]), 2)
        self.assertEqual(num_dist_cy(data_1[0], data_1[1]), -1)
        self.assertEqual(num_dist_cy(data_1[0], data_1[1]), 0)
    
    def test_subtree_root(self):
        data = [7, 17, 18]
        self.assertEqual(_subtree_root(data), 3)
        
    def test_parent(self):
        data_1 = 7
        # root is only node without a parent
        # and a value of -1 is returned
        data_2 = 0
        self.assertEqual(_parent(data_1), 3)
        self.assertEqual(_parent(data_2), -1)
        
    def test_sibling(self):
        # odd value location
        data_1 = 7
        # even value location
        data_2 = 4
        # root is only node without a sibling
        # and a value of -1 is returned
        data_3 = 0
        self.assertEqual(_sibling(data_1), 8)
        self.assertEqual(_sibling(data_2), 3)
        self.assertEqual(_sibling(data_3), -1)
        
    def test_ancestors(self):
        data_node_1 = 12
        # only root as ancestor
        data_node_2 = 1
        # root has no ancestor
        # and returns empty list
        data_node_3 = 0
        self.assertEqual(_ancestors(data_node_1, self.ordered_array)[0], 5)
        self.assertEqual(_ancestors(data_node_1, self.ordered_array)[1], 2)
        self.assertEqual(_ancestors(data_node_1, self.ordered_array)[2], 0)
        self.assertEqual(_ancestors(data_node_2, self.ordered_array)[0], 0)
        self.assertEqual(_ancestors(data_node_3, self.ordered_array).size, 0)
        
    def test_subtree(self):
        # subtree is an internal node
        data_node_1 = 12
        # subtree is a leaf
        data_node_2 = 26
        self.assertEqual(_subtree(data_node_1, self.ordered_array)[0], 25)
        self.assertEqual(_subtree(data_node_1, self.ordered_array)[1], 26)
        self.assertEqual(_subtree(data_node_1, self.ordered_array)[2], 12)
        self.assertEqual(_subtree(data_node_2, self.ordered_array)[0], 26)
        
    def test_move_subtree(self):
        # subtree of node at position 6
        data_subtree = [13, 14]
        # old and new lca's
        data_lca = [6, 13]
        actual_tree = self.tree_array
        _move_subtree(actual_tree, data_subtree, data_lca[0], data_lca[1])
        self.assertEqual(actual_tree, self.expected_array)
        
    def test_move_node(self):
        # move to sibling position
        data_1 = [7, 3, 4]
        # move down the tree
        data_2 = [7, 3, 1]
        # move up the tree
        data_3 = [7, 3, 15]
        self.assertEqual(_move_node(data_1[0], data_1[1], data_1[2]), 9)
        self.assertEqual(_move_node(data_2[0], data_2[1], data_2[2]), 3)
        self.assertEqual(_move_node(data_3[0], data_3[1], data_3[2]), 31)
        
    def test_array_to_tree(self):
        data_taxa = ['a', 'b', 'c', 'd', 'e', 'f']
        actual_tree = _array_to_tree(data_taxa, self.tree_array)
        self.assertEqual(actual_tree, self.expected_TreeNode)


if __name__ == "__main__":
    main()
