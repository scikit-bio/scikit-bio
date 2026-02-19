# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd
from skbio.tree import TreeNode
from skbio.tree._dataframe import (dataframe_to_treenode, treenode_to_dataframe,
                                   tip_to_root_conversion)


class TestTreeNodeDataFrameConversion(unittest.TestCase):
    """Test conversion between TreeNode and DataFrame representations."""
    
    def setUp(self):
        """Set up test data with a simple tree."""
        # Create a simple test tree in memory
        #          /c
        #          │
        #      /b──┤ 
        #      │   │
        # ──a──┤   \d
        #      │
        #      \e
        node_a = TreeNode('a')
        node_b = TreeNode('b')
        node_c = TreeNode('c') 
        node_d = TreeNode('d')
        node_e = TreeNode('e')
        node_a.append(node_b)
        node_b.append(node_c)
        node_b.append(node_d)
        node_a.append(node_e)
        self.test_tree = node_a
        
    def test_treenode_to_dataframe_conversion(self):
        """Test converting TreeNode to DataFrame."""
        df = treenode_to_dataframe(self.test_tree)
        
        # Check that we have the expected structure
        expected_columns = ['parent', 'node', 'name', 'length', 'is_tip']
        self.assertEqual(list(df.columns), expected_columns)
        
        # Check that we have the right number of nodes
        self.assertEqual(len(df), 5)  # a, b, c, d, e
        
        # Check that node names are preserved
        names = set(df['name'].tolist())
        expected_names = {'a', 'b', 'c', 'd', 'e'}
        self.assertEqual(names, expected_names)
        
    def test_dataframe_to_treenode_conversion(self):
        """Test converting DataFrame back to TreeNode."""
        # Convert to dataframe and back
        df = treenode_to_dataframe(self.test_tree)
        reconstructed_tree = dataframe_to_treenode(df)
        
        # Check basic properties are preserved
        original_tip_count = self.test_tree.count(tips=True)
        reconstructed_tip_count = reconstructed_tree.count(tips=True)
        self.assertEqual(original_tip_count, reconstructed_tip_count)
        
        # Check that tip names are preserved
        original_tips = {tip.name for tip in self.test_tree.tips()}
        reconstructed_tips = {tip.name for tip in reconstructed_tree.tips()}
        self.assertEqual(original_tips, reconstructed_tips)
        
    def test_roundtrip_conversion(self):
        """Test that tree -> dataframe -> tree preserves structure."""
        df = treenode_to_dataframe(self.test_tree)
        reconstructed = dataframe_to_treenode(df)
        
        # The reconstructed tree should have the same topology
        # (RF distance of 0 indicates identical topology)
        rf_distance = self.test_tree.compare_rfd(reconstructed)
        self.assertEqual(rf_distance, 0.0, 
                        f"Expected identical topology but RF distance was {rf_distance}")


class TestTipToRootConversion(unittest.TestCase):
    """Test tip-to-root conversion functionality."""
    
    def setUp(self):
        """Set up test data with a larger tree."""
        # Create a more complex test tree
        # ((c,d)b,e)a
        node_a = TreeNode('a')
        node_b = TreeNode('b') 
        node_c = TreeNode('c')
        node_d = TreeNode('d')
        node_e = TreeNode('e')
        node_a.append(node_b)
        node_a.append(node_e)
        node_b.append(node_c)
        node_b.append(node_d)
        self.tree_root = node_a
        
    def test_tip_to_root_with_all_tips(self):
        """Test tip-to-root conversion when all tips are included."""
        df = treenode_to_dataframe(self.tree_root)
        all_tips = [tip.name for tip in self.tree_root.tips()]
        
        reconstructed = tip_to_root_conversion(df, all_tips)
        
        # Should be identical to original
        rf_distance = self.tree_root.compare_rfd(reconstructed)
        self.assertEqual(rf_distance, 0.0)
        
    def test_tip_to_root_subset(self):
        """Test tip-to-root conversion with a subset of tips.""" 
        df = treenode_to_dataframe(self.tree_root)
        subset_tips = ['c', 'd']  # Just these two tips
        
        reconstructed = tip_to_root_conversion(df, subset_tips)
        
        # Check that requested tips are present
        reconstructed_tip_names = {tip.name for tip in reconstructed.tips()}
        for tip in subset_tips:
            self.assertIn(tip, reconstructed_tip_names)
            
        # Check that we have the right number of nodes (should include ancestors)
        self.assertGreaterEqual(len(list(reconstructed.traverse())), len(subset_tips))


if __name__ == '__main__':
    unittest.main()