# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np

from skbio import TreeNode
from skbio.tree import majority_rule
from skbio.tree._majority_rule import (_walk_clades, _filter_clades,
                                       _build_trees)


class MajorityRuleTests(TestCase):
    def test_majority_rule(self):
        trees = [
            TreeNode.read(
                io.StringIO("(A,(B,(H,(D,(J,(((G,E),(F,I)),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(D,((J,H),(((G,E),(F,I)),C)))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(D,(H,(J,(((G,E),(F,I)),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,(G,((F,I),((J,(H,D)),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,(G,((F,I),(((J,H),D),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,((F,I),(G,((J,(H,D)),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,((F,I),(G,(((J,H),D),C))))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,((G,(F,I)),((J,(H,D)),C)))));")),
            TreeNode.read(
                io.StringIO("(A,(B,(E,((G,(F,I)),(((J,H),D),C)))));"))]

        exp = TreeNode.read(io.StringIO("(((E,(G,(F,I),(C,(D,J,H)))),B),A);"))
        obs = majority_rule(trees)
        self.assertEqual(exp.compare_subsets(obs[0]), 0.0)
        self.assertEqual(len(obs), 1)

        tree = obs[0]
        exp_supports = sorted([9.0, 9.0, 9.0, 6.0, 6.0, 6.0])
        obs_supports = sorted([n.support for n in tree.non_tips()])
        self.assertEqual(obs_supports, exp_supports)

        obs = majority_rule(trees, weights=np.ones(len(trees)) * 2)
        self.assertEqual(exp.compare_subsets(obs[0]), 0.0)
        self.assertEqual(len(obs), 1)

        tree = obs[0]
        exp_supports = sorted([18.0, 18.0, 12.0, 18.0, 12.0, 12.0])
        obs_supports = sorted([n.support for n in tree.non_tips()])

        with self.assertRaises(ValueError):
            majority_rule(trees, weights=[1, 2])

    def test_majority_rule_multiple_trees(self):
        trees = [
            TreeNode.read(io.StringIO("((a,b),(c,d),(e,f));")),
            TreeNode.read(io.StringIO("(a,(c,d),b,(e,f));")),
            TreeNode.read(io.StringIO("((c,d),(e,f),b);")),
            TreeNode.read(io.StringIO("(a,(c,d),(e,f));"))]

        trees = majority_rule(trees)
        self.assertEqual(len(trees), 4)

        for tree in trees:
            self.assertIs(type(tree), TreeNode)

        exp = set([
                  frozenset(['a']),
                  frozenset(['b']),
                  frozenset([None, 'c', 'd']),
                  frozenset([None, 'e', 'f'])])

        obs = set([frozenset([n.name for n in t.traverse()]) for t in trees])
        self.assertEqual(obs, exp)

    def test_majority_rule_tree_node_class(self):
        class TreeNodeSubclass(TreeNode):
            pass

        trees = [
            TreeNode.read(io.StringIO("((a,b),(c,d),(e,f));")),
            TreeNode.read(io.StringIO("(a,(c,d),b,(e,f));")),
            TreeNode.read(io.StringIO("((c,d),(e,f),b);")),
            TreeNode.read(io.StringIO("(a,(c,d),(e,f));"))]

        trees = majority_rule(trees, tree_node_class=TreeNodeSubclass)
        self.assertEqual(len(trees), 4)

        for tree in trees:
            self.assertIs(type(tree), TreeNodeSubclass)

        exp = set([
                  frozenset(['a']),
                  frozenset(['b']),
                  frozenset([None, 'c', 'd']),
                  frozenset([None, 'e', 'f'])])

        obs = set([frozenset([n.name for n in t.traverse()]) for t in trees])
        self.assertEqual(obs, exp)

    def test_walk_clades(self):
        trees = [TreeNode.read(io.StringIO("((A,B),(D,E));")),
                 TreeNode.read(io.StringIO("((A,B),(D,(E,X)));"))]
        exp_clades = [
            (frozenset(['A']), 2.0),
            (frozenset(['B']), 2.0),
            (frozenset(['A', 'B']), 2.0),
            (frozenset(['D', 'E']), 1.0),
            (frozenset(['D', 'E', 'A', 'B']), 1.0),
            (frozenset(['D']), 2.0),
            (frozenset(['E']), 2.0),
            (frozenset(['X']), 1.0),
            (frozenset(['E', 'X']), 1.0),
            (frozenset(['D', 'E', 'X']), 1.0),
            (frozenset(['A', 'B', 'D', 'E', 'X']), 1.0)]

        exp_lengths_nolength = {
            frozenset(['A']): None,
            frozenset(['B']): None,
            frozenset(['A', 'B']): None,
            frozenset(['D', 'E']): None,
            frozenset(['D', 'E', 'A', 'B']): None,
            frozenset(['D']): None,
            frozenset(['E']): None,
            frozenset(['X']): None,
            frozenset(['E', 'X']): None,
            frozenset(['D', 'E', 'X']): None,
            frozenset(['A', 'B', 'D', 'E', 'X']): None}

        exp_lengths = {
            frozenset(['A']): 2.0,
            frozenset(['B']): 2.0,
            frozenset(['A', 'B']): 2.0,
            frozenset(['D', 'E']): 1.0,
            frozenset(['D', 'E', 'A', 'B']): 1.0,
            frozenset(['D']): 2.0,
            frozenset(['E']): 2.0,
            frozenset(['X']): 1.0,
            frozenset(['E', 'X']): 1.0,
            frozenset(['D', 'E', 'X']): 1.0,
            frozenset(['A', 'B', 'D', 'E', 'X']): 1.0}

        obs_clades, obs_lengths = _walk_clades(trees, np.ones(len(trees)))
        self.assertEqual(set(obs_clades), set(exp_clades))
        self.assertEqual(obs_lengths, exp_lengths_nolength)

        for t in trees:
            for n in t.traverse(include_self=True):
                n.length = 2.0

        obs_clades, obs_lengths = _walk_clades(trees, np.ones(len(trees)))

        self.assertEqual(set(obs_clades), set(exp_clades))
        self.assertEqual(obs_lengths, exp_lengths)

    def test_filter_clades(self):
        clade_counts = [(frozenset(['A', 'B']), 8),
                        (frozenset(['A', 'C']), 7),
                        (frozenset(['A']), 6),
                        (frozenset(['B']), 5)]
        obs = _filter_clades(clade_counts, 2)
        exp = {frozenset(['A', 'B']): 8,
               frozenset(['A']): 6,
               frozenset(['B']): 5}
        self.assertEqual(obs, exp)

        clade_counts = [(frozenset(['A']), 8),
                        (frozenset(['B']), 7),
                        (frozenset(['C']), 7),
                        (frozenset(['A', 'B']), 6),
                        (frozenset(['A', 'B', 'C']), 5),
                        (frozenset(['D']), 2)]
        obs = _filter_clades(clade_counts, 4)
        exp = {frozenset(['A']): 8,
               frozenset(['B']): 7,
               frozenset(['C']): 7,
               frozenset(['A', 'B']): 6,
               frozenset(['A', 'B', 'C']): 5}
        self.assertEqual(obs, exp)

    def test_build_trees(self):
        clade_counts = {frozenset(['A', 'B']): 6,
                        frozenset(['A']): 7,
                        frozenset(['B']): 8}
        edge_lengths = {frozenset(['A', 'B']): 1,
                        frozenset(['A']): 2,
                        frozenset(['B']): 3}
        tree = _build_trees(clade_counts, edge_lengths, 'foo', TreeNode)[0]
        self.assertEqual(tree.foo, 6)
        tree_foos = set([c.foo for c in tree.children])
        tree_lens = set([c.length for c in tree.children])
        self.assertEqual(tree_foos, set([7, 8]))
        self.assertEqual(tree_lens, set([2, 3]))


if __name__ == '__main__':
    main()
