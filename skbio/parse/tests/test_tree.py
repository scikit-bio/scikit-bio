#!/usr/bin/env python
"""Unit tests for tree parsers"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.parse.tree import DndParser, DndTokenizer, RecordError
from bipy.core.tree import TreeNode
from bipy.util.unit_test import TestCase, main

__credits__ = ["Rob Knight", "Peter Maxwell", "Daniel McDonald"]

sample = """
(
(
xyz:0.28124,
(
def:0.24498,
mno:0.03627)
:0.17710)
:0.04870,

abc:0.05925,
(
ghi:0.06914,
jkl:0.13776)
:0.09853);
"""

node_data_sample = """
(
(
xyz:0.28124,
(
def:0.24498,
mno:0.03627)
'A':0.17710)
B:0.04870,

abc:0.05925,
(
ghi:0.06914,
jkl:0.13776)
C:0.09853);
"""

minimal = "();"
no_names = "((,),(,));"
missing_tip_name = "((a,b),(c,));"

empty = '();'
single = '(abc:3);'
double = '(abc:3, def:4);'
onenest = '(abc:3, (def:4, ghi:5):6 );'
nodedata = '(abc:3, (def:4, ghi:5)jkl:6 );'

class DndTokenizerTests(TestCase):
    """Tests of the DndTokenizer factory function."""

    def test_gdata(self):
        """DndTokenizer should work as expected on real data"""
        exp = \
        ['(', '(', 'xyz', ':', '0.28124',',', '(', 'def', ':', '0.24498',\
        ',', 'mno', ':', '0.03627', ')', ':', '0.17710', ')', ':', '0.04870', \
        ',', 'abc', ':', '0.05925', ',', '(', 'ghi', ':', '0.06914', ',', \
        'jkl', ':', '0.13776', ')', ':', '0.09853', ')', ';']
        #split it up for debugging on an item-by-item basis
        obs = list(DndTokenizer(sample))
        self.assertEqual(len(obs), len(exp))
        for i, j in zip(obs, exp):
            self.assertEqual(i, j)
        #try it all in one go
        self.assertEqual(list(DndTokenizer(sample)), exp)

    def test_nonames(self):
        """DndTokenizer should work as expected on trees with no names"""
        exp = ['(','(',',',')',',','(',',',')',')',';']
        obs = list(DndTokenizer(no_names))
        self.assertEqual(obs, exp)

    def test_missing_tip_name(self):
        """DndTokenizer should work as expected on trees with a missing name"""
        exp = ['(','(','a',',','b',')',',','(','c',',',')',')',';']
        obs = list(DndTokenizer(missing_tip_name))
        self.assertEqual(obs, exp)

    def test_minimal(self):
        """DndTokenizer should work as expected a minimal tree without names"""
        exp = ['(',')',';']
        obs = list(DndTokenizer(minimal))
        self.assertEqual(obs, exp)

class DndParserTests(TestCase):
    """Tests of the DndParser factory function."""
    def test_nonames(self):
        """DndParser should produce the correct tree when there are no names"""
        obs = DndParser(no_names)
        exp = TreeNode()
        exp.append(TreeNode())
        exp.append(TreeNode())
        exp.Children[0].append(TreeNode())
        exp.Children[0].append(TreeNode())
        exp.Children[1].append(TreeNode())
        exp.Children[1].append(TreeNode())
        self.assertEqual(str(obs), str(exp))

    def test_minimal(self):
        """DndParser should produce the correct minimal tree"""
        obs = DndParser(minimal)
        exp = TreeNode()
        exp.append(TreeNode())
        self.assertEqual(str(obs), str(exp))

    def test_missing_tip_name(self):
        """DndParser should produce the correct tree when missing a name"""
        obs = DndParser(missing_tip_name)
        exp = TreeNode()
        exp.append(TreeNode())
        exp.append(TreeNode())
        exp.Children[0].append(TreeNode(Name='a'))
        exp.Children[0].append(TreeNode(Name='b'))
        exp.Children[1].append(TreeNode(Name='c'))
        exp.Children[1].append(TreeNode())
        self.assertEqual(str(obs), str(exp))

    def test_gsingle(self):
        """DndParser should produce a single-child TreeNode on minimal data"""
        t = DndParser(single)
        self.assertEqual(len(t), 1)
        child = t[0]
        self.assertEqual(child.Name, 'abc')
        self.assertEqual(child.Length, 3)
        self.assertEqual(str(t), '(abc:3.0);')

    def test_gdouble(self):
        """DndParser should produce a double-child TreeNode from data"""
        t = DndParser(double)
        self.assertEqual(len(t), 2)
        self.assertEqual(str(t), '(abc:3.0,def:4.0);')

    def test_gonenest(self):
        """DndParser should work correctly with nested data"""
        t = DndParser(onenest)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  #first child is terminal
        self.assertEqual(len(t[1]), 2)  #second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0):6.0);')

    def test_gnodedata(self):
        """DndParser should assign Name to internal nodes correctly"""
        t = DndParser(nodedata)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  #first child is terminal
        self.assertEqual(len(t[1]), 2)  #second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0)jkl:6.0);')
        info_dict = {}
        for node in t.traverse():
            info_dict[node.Name] = node.Length
        self.assertEqual(info_dict['abc'], 3.0)
        self.assertEqual(info_dict['def'], 4.0)
        self.assertEqual(info_dict['ghi'], 5.0)
        self.assertEqual(info_dict['jkl'], 6.0)

    def test_data(self):
        """DndParser should work as expected on real data"""
        t = DndParser(sample)
        self.assertEqual(str(t), '((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);')
        tdata = DndParser(node_data_sample, unescape_name=True)
        self.assertEqual(str(tdata), "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)C:0.09853);")

    def test_gbad(self):
        """DndParser should fail if parens unbalanced"""
        left = '((abc:3)'
        right = '(abc:3))'
        self.assertRaises(RecordError, DndParser, left)
        self.assertRaises(RecordError, DndParser, right)

    def test_DndParser(self):
        """DndParser tests"""
        t_str = "(A_a,(B:1.0,C),'D_e':0.5)E;"
        tree_unesc = DndParser(t_str, TreeNode, unescape_name=True)
        tree_esc = DndParser(t_str, TreeNode, unescape_name=False)

        self.assertEqual(tree_unesc.Name, 'E')
        self.assertEqual(tree_unesc.Children[0].Name, 'A a')
        self.assertEqual(tree_unesc.Children[1].Children[0].Name, 'B')
        self.assertEqual(tree_unesc.Children[1].Children[0].Length, 1.0)
        self.assertEqual(tree_unesc.Children[1].Children[1].Name, 'C')
        self.assertEqual(tree_unesc.Children[2].Name, 'D_e')
        self.assertEqual(tree_unesc.Children[2].Length, 0.5)

        self.assertEqual(tree_esc.Name, 'E')
        self.assertEqual(tree_esc.Children[0].Name, 'A_a')
        self.assertEqual(tree_esc.Children[1].Children[0].Name, 'B')
        self.assertEqual(tree_esc.Children[1].Children[0].Length, 1.0)
        self.assertEqual(tree_esc.Children[1].Children[1].Name, 'C')
        self.assertEqual(tree_esc.Children[2].Name, "'D_e'")
        self.assertEqual(tree_esc.Children[2].Length, 0.5)

        reload_test = tree_esc.to_newick(with_distances=True, \
                                         escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True), \
                         tree_esc.to_newick(with_distances=True))
        reload_test = tree_unesc.to_newick(with_distances=True, \
                                           escape_name=False)
        obs = DndParser(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True), \
                         tree_unesc.to_newick(with_distances=True))

class TreeNodeTests(TestCase):
    """Check that TreeNode works the way I think"""
    def test_gops(self):
        """Basic TreeNode operations should work as expected"""
        p = TreeNode()
        self.assertEqual(str(p), ';')
        p.Name = 'abc'
        self.assertEqual(str(p), 'abc;')
        p.Length = 3
        self.assertEqual(str(p), 'abc:3;')   #don't suppress branch from root
        q = TreeNode()
        p.append(q)
        self.assertEqual(str(p), '()abc:3;')
        r = TreeNode()
        q.append(r)
        self.assertEqual(str(p), '(())abc:3;')
        r.Name = 'xyz'
        self.assertEqual(str(p), '((xyz))abc:3;')
        q.Length = 2
        self.assertEqual(str(p), '((xyz):2)abc:3;')

if __name__ == '__main__':
    main()
