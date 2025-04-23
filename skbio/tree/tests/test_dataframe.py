import skbio
from skbio import TreeNode
from unittest import TestCase
from skbio.tree._dataframe import treenode_to_dataframe,dataframe_to_treenode,tip_to_root_conversion
import pandas as pd

import cProfile
import re

class unit_test_dataframe_to_treenode_85_otsu_tree(TestCase):

    def setUp(self):
        node_a = TreeNode('a')
        node_b = TreeNode('b')
        node_c = TreeNode('c')
        node_d = TreeNode('d')
        node_e = TreeNode('e')
        node_a.append(node_b)
        node_b.append(node_c)
        node_b.append(node_d)
        node_a.append(node_e)
        self.TreeRoot = node_a
        #
        #          /c
        #          |
        #      /b--|
        #      |   |
        # --a---|   \d
        #      |
        #      \e

        list_df = [{'parent': -1, 'node': 1, 'name': 'a', 'length': None},
                   {'parent': 1, 'node': 2, 'name': 'b', 'length': None},
                   {'parent': 2, 'node': 3, 'name': 'c', 'length': None},
                   {'parent': 2, 'node': 4, 'name': 'd', 'length': None},
                   {'parent': 1, 'node': 5, 'name': 'e', 'length': None}
                   ]

        self.expected_tree = skbio.TreeNode.read('enter file path', format='newick')
        observed_tree = treenode_to_dataframe(self.expected_tree)
        self.observed_tree = dataframe_to_treenode(observed_tree)

    def test_count_tips(self):
        real_tree_tips = self.expected_tree.count(tips=True)
        test_tree_tips = self.observed_tree.count(tips=True)
        self.assertEqual(test_tree_tips, real_tree_tips)

    def test_root(self):
        real_root = self.expected_tree.root().name
        test_tree_root = self.observed_tree.root().name
        self.assertEqual(test_tree_root, real_root)

    def test_tip(self):
        real_tips = {n.name for n in self.expected_tree.tips()}
        # test_tree_tips = list(self.test_tree.tips())
        est_tree_tips = {n.name for n in self.observed_tree.tips()}
        self.assertEqual(real_tips, est_tree_tips)

        # print(self.test_tree.ascii_art())

    def test_rfd(self):
        # self.expected_tree.assign_ids()
        # print(type(self.expected_tree))
        self.assertEqual(self.expected_tree.compare_rfd(self.observed_tree), 0.0)
        # print(self.expected_tree)

    def test_tip_distance(self):
        self.expected_tree.assign_ids()
        self.assertEqual(self.expected_tree.compare_tip_distances(self.observed_tree, sample=1000), 0.0)

    def test_paths(self):
        observed_tips = {tip.name for tip in self.observed_tree.tips()}
        # expeted_tips = {expeted_tips for expeted_tips in self.expected_tree.tips()}
        # observed_tips = {observed_tips for observed_tips in self.observed_tree.tips()}
        expected_tips = {tip.name for tip in self.expected_tree.tips()}
        if expected_tips != observed_tips:
            self.assertEqual(0, 1)
        for exp_tips in expected_tips:
            e_tips = self.expected_tree.find(exp_tips)
            o_tips = self.observed_tree.find(exp_tips)
            ancestor_exp = [node.name for node in e_tips.ancestors()]
            ancestor_obs = [node.name for node in o_tips.ancestors()]
            if ancestor_exp != ancestor_obs:
                self.assertEqual(0, 1)
        self.assertEqual(1, 1)

class unit_test_dataframe_to_treenode_simple_tree(TestCase):
    def setUp(self):
        node_a = TreeNode('a')
        node_b = TreeNode(name = 'b',length =2)
        node_c = TreeNode('c')
        node_d = TreeNode(name = 'd',length =1 )
        node_e = TreeNode('e')
        node_a.append(node_b)
        node_b.append(node_c)
        node_b.append(node_d)
        node_a.append(node_e)
        self.TreeRoot = node_a
        #
        #          /c
        #          |
        #      /b--|
        #      |   |
        # --a---|   \d
        #      |
        #      \e

        list_df = [{'parent': -1, 'node': 1, 'name': 'a', 'length': None},
                   {'parent': 1, 'node': 2, 'name': 'b', 'length': 2},
                   {'parent': 2, 'node': 3, 'name': 'c', 'length': None},
                   {'parent': 2, 'node': 4, 'name': 'd', 'length': 1},
                   {'parent': 1, 'node': 5, 'name': 'e', 'length': None}
                   ]
        observed_tree = pd.DataFrame(list_df)
        self.observed_tree = dataframe_to_treenode(observed_tree)

    def test_count_tips(self):
        real_tree_tips = self.TreeRoot.count(tips=True)
        test_tree_tips = self.observed_tree.count(tips=True)
        self.assertEqual(test_tree_tips, real_tree_tips)

    def test_root(self):
        real_root = self.TreeRoot.root().name
        test_tree_root = self.observed_tree.root().name
        self.assertEqual(test_tree_root, real_root)

    def test_tip(self):
        real_tips = {n.name for n in self.TreeRoot.tips()}
        # test_tree_tips = list(self.test_tree.tips())
        est_tree_tips = {n.name for n in self.observed_tree.tips()}
        self.assertEqual(real_tips, est_tree_tips)

        # print(self.test_tree.ascii_art())

    def test_rfd(self):
        # self.expected_tree.assign_ids()
        # print(type(self.expected_tree))
        self.assertEqual(self.TreeRoot.compare_rfd(self.observed_tree), 0.0)
        # print(self.expected_tree)

    def test_tip_distance(self):
        self.TreeRoot.assign_ids()
        self.assertEqual(self.TreeRoot.compare_tip_distances(self.observed_tree), 0.0) # comes out nan if no length is provided

    def test_paths(self):
        observed_tips = {tip.name for tip in self.observed_tree.tips()}
        # expeted_tips = {expeted_tips for expeted_tips in self.expected_tree.tips()}
        # observed_tips = {observed_tips for observed_tips in self.observed_tree.tips()}
        expected_tips = {tip.name for tip in self.TreeRoot.tips()}
        if expected_tips != observed_tips:
            self.assertEqual(0, 1)
        for exp_tips in expected_tips:
            e_tips = self.TreeRoot.find(exp_tips)
            o_tips = self.observed_tree.find(exp_tips)
            ancestor_exp = [node.name for node in e_tips.ancestors()]
            ancestor_obs = [node.name for node in o_tips.ancestors()]
            if ancestor_exp != ancestor_obs:
                self.assertEqual(0, 1)
        self.assertEqual(1, 1)



class unit_test_treenode_to_dataframe(TestCase):

    def setUp(self):
        node_a = TreeNode('a')
        self.test_tree = node_a
        node_b = TreeNode('b')
        node_c = TreeNode('c')
        node_d = TreeNode('d')
        node_e = TreeNode('e')
        node_a.append(node_b)
        node_b.append(node_c)
        node_b.append(node_d)
        node_a.append(node_e)
        self.TreeRoot = node_a
        #
        #          /c
        #          |
        #      /b--|
        #      |   |
        # --a---|   \d
        #      |
        #      \e

        node_a.assign_ids()
        self.data_frame = treenode_to_dataframe(node_a)

        self.parent_list = []
        self.node_list = []
        for node in node_a.traverse(include_self=True):
            self.parent_list.append(node.parent.id if node.parent is not None else -1)
            self.node_list.append(node.id)

    def test_parent(self):
        test_parent_list = [x for x in self.data_frame['parent']]
        self.assertEqual(test_parent_list, self.parent_list)

    def test_node(self):
        test_node_list = [y for y in self.data_frame['node']]
        self.assertEqual(self.node_list, test_node_list)

    def test_parent_size(self):
        parent_len = len(self.data_frame['parent'])
        self.assertEqual(parent_len, 5)

    def test_node_size(self):
        node_len = len(self.data_frame['node'])
        self.assertEqual(node_len, 5)


class TestTipToRootConversion(TestCase):
    def setUp(self):
        # Creating a simple test tree
        node_a = TreeNode(name='a')
        node_b = TreeNode(name='b')
        node_c = TreeNode(name='c')
        node_d = TreeNode(name='d')
        node_a.append(node_b)
        node_a.append(node_c)
        node_b.append(node_d)

        #      /b----d
        # a----|
        #      \c
        self.tree_root = node_a  # Store the original tree
        self.tip_name = "d"

        # Convert TreeNode to DataFrame before calling the function
        self.df_tree = treenode_to_dataframe(self.tree_root)

        # Call the function with the correct DataFrame format
        self.reconstructed_tree = tip_to_root_conversion(self.df_tree, [self.tip_name])

    def test_root_is_correct(self):
        #Test that the reconstructed tree has the correct root.
        self.assertEqual(self.reconstructed_tree.root().name, self.tree_root.root().name)

    def test_tip_is_present(self):
        #Test that the expected tip is present in the reconstructed tree.
        tip_names = {tip.name for tip in self.reconstructed_tree.tips()}
        self.assertIn(self.tip_name, tip_names)

    def test_node_count(self):
        #Test that the number of nodes in the reconstructed tree is as expected.
        expected_node_count = 3  # Nodes a, b, d
        actual_node_count = self.reconstructed_tree.count()
        self.assertEqual(actual_node_count, expected_node_count)

    def test_all_nodes_present_in_reconstructed_tree(self):
        tips_to_keep = ['c', 'd']
        reconstructed = tip_to_root_conversion(self.df_tree, tips_to_keep)

        # Create expected tree by shearing original tree without pruning
        expected = self.tree_root.shear(tips_to_keep, prune=False)

        # Collect all node names from both trees (internal + tips)
        expected_node_names = {n.name for n in expected.traverse() if n.name is not None}
        reconstructed_node_names = {n.name for n in reconstructed.traverse() if n.name is not None}

        # Assert all expected nodes are present
        self.assertEqual(expected_node_names, reconstructed_node_names)

    def test_rfd_of_reconstructed_tree(self):
        # Tips to reconstruct from
        tips_to_keep = ['c', 'd']

        # Reconstruct from DataFrame
        reconstructed = tip_to_root_conversion(self.df_tree, tips_to_keep)

        # Generate expected tree via shear without pruning
        expected = self.tree_root.shear(tips_to_keep, prune=False)

        # Compare RF distance (should be 0.0 if topologies are identical)
        rf_distance = expected.compare_rfd(reconstructed)
        self.assertEqual(rf_distance, 0.0, f"Expected RF distance of 0.0 but got {rf_distance}")


class TestTipToRootConversion_85OTUs(TestCase):
    def setUp(self):
        # Load tree from file
        self.tree_path = 'enter file path'
        self.tree_root = TreeNode.read(self.tree_path, format='newick')
        self.tips_to_keep = ['give names of tip as string']

        # Convert tree to dataframe
        self.df_tree = treenode_to_dataframe(self.tree_root)

        # Reconstruct tree from the selected tips
        self.reconstructed_tree = tip_to_root_conversion(self.df_tree, self.tips_to_keep)

        # Expected tree using shear with prune=False
        self.expected_tree = self.tree_root.shear(self.tips_to_keep, prune=False)

    def test_all_nodes_present(self):
        #Ensure all nodes in expected tree are present in reconstructed tree.
        expected_names = {n.name for n in self.expected_tree.traverse() if n.name is not None}
        reconstructed_names = {n.name for n in self.reconstructed_tree.traverse() if n.name is not None}
        self.assertEqual(expected_names, reconstructed_names)

    def test_rf_distance(self):
        #RF distance between expected and reconstructed should be 0.0 (topologically identical).
        rf_distance = self.expected_tree.compare_rfd(self.reconstructed_tree)
        self.assertEqual(rf_distance, 0.0, f"Expected RF distance of 0.0 but got {rf_distance}")

    def test_tips_present(self):
        #Check that the tips we requested are present in the reconstructed tree.
        reconstructed_tips = {tip.name for tip in self.reconstructed_tree.tips()}
        for tip in self.tips_to_keep:
            self.assertIn(tip, reconstructed_tips)

    def test_branch_length_preservation(self):
        tips_to_keep = ['enter tip name as list']
        reconstructed = tip_to_root_conversion(self.df_tree, tips_to_keep)

        # Compare branch lengths for nodes with matching names
        for node in reconstructed.traverse():
            if node.name:
                orig_node = self.tree_root.find(node.name)
                self.assertEqual(node.length, orig_node.length,
                                 f"Branch length mismatch at node {node.name}")

    def test_ancestor_paths_match(self):
        tips = ['enter tip name as list']
        reconstructed = tip_to_root_conversion(self.df_tree, tips)

        for tip in tips:
            # Find the ancestor paths for the original tree
            original_path = [n.name for n in self.tree_root.find(tip).ancestors()]
            # Find the ancestor paths for the reconstructed tree
            new_path = [n.name for n in reconstructed.find(tip).ancestors()]

            # Assert that both ancestor paths match
            self.assertEqual(original_path, new_path,
                             f"Ancestor paths for tip {tip} do not match: original path {original_path}, new path {new_path}")

    def test_full_tree_reconstruction(self):
        # Get all tips from the original tree
        all_tips = [tip.name for tip in self.tree_root.tips()]

        # Reconstruct the tree using all the tips
        reconstructed = tip_to_root_conversion(self.df_tree, all_tips)

        # Calculate the RF distance between the original and reconstructed trees
        rf_distance = self.tree_root.compare_rfd(reconstructed)

        # Assert that the RF distance is 0, indicating the trees are identical
        self.assertEqual(rf_distance, 0.0, f"RF distance is not 0.0; it is {rf_distance}")