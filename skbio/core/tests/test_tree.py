#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import numpy.testing as nptest
from unittest import TestCase, main
from skbio.core.tree import TreeNode, _dnd_tokenizer
from skbio.core.exception import NoLengthError, TreeError, RecordError
from skbio.maths.stats.test import correlation_t

class TreeTests(TestCase):
    def setUp(self):
        """Prep the self"""
        self.simple_t = TreeNode.from_newick("((a,b)i1,(c,d)i2)root;")
        nodes = dict([(x, TreeNode(x)) for x in 'abcdefgh'])
        nodes['a'].append(nodes['b'])
        nodes['b'].append(nodes['c'])
        nodes['c'].append(nodes['d'])
        nodes['c'].append(nodes['e'])
        nodes['c'].append(nodes['f'])
        nodes['f'].append(nodes['g'])
        nodes['a'].append(nodes['h'])
        self.TreeNode = nodes
        self.TreeRoot = nodes['a']

    def test_copy(self):
        """copy a tree"""
        self.simple_t.Children[0].Length = 1.2
        self.simple_t.Children[1].Children[0].Length = 0.5
        cp = self.simple_t.copy()
        gen = zip(cp.traverse(include_self=True),
                  self.simple_t.traverse(include_self=True))

        for a,b in gen:
            self.assertIsNot(a, b)
            self.assertEqual(a.Name, b.Name)
            self.assertEqual(a.Length, b.Length)

    def test_append(self):
        """Append a node to a tree"""
        second_tree = TreeNode.from_newick("(x,y)z;")
        self.simple_t.append(second_tree)

        self.assertEqual(self.simple_t.Children[0].Name, 'i1')
        self.assertEqual(self.simple_t.Children[1].Name, 'i2')
        self.assertEqual(self.simple_t.Children[2].Name, 'z')
        self.assertEqual(len(self.simple_t.Children), 3)
        self.assertEqual(self.simple_t.Children[2].Children[0].Name, 'x')
        self.assertEqual(self.simple_t.Children[2].Children[1].Name, 'y')
        self.assertEqual(second_tree.Parent, self.simple_t)

    def test_extend(self):
        """Extend a few nodes"""
        second_tree = TreeNode.from_newick("(x1,y1)z1;")
        third_tree = TreeNode.from_newick("(x2,y2)z2;")
        self.simple_t.extend([second_tree, third_tree])

        self.assertEqual(self.simple_t.Children[0].Name, 'i1')
        self.assertEqual(self.simple_t.Children[1].Name, 'i2')
        self.assertEqual(self.simple_t.Children[2].Name, 'z1')
        self.assertEqual(self.simple_t.Children[3].Name, 'z2')
        self.assertEqual(len(self.simple_t.Children), 4)
        self.assertEqual(self.simple_t.Children[2].Children[0].Name, 'x1')
        self.assertEqual(self.simple_t.Children[2].Children[1].Name, 'y1')
        self.assertEqual(self.simple_t.Children[3].Children[0].Name, 'x2')
        self.assertEqual(self.simple_t.Children[3].Children[1].Name, 'y2')
        self.assertIs(second_tree.Parent, self.simple_t)
        self.assertIs(third_tree.Parent, self.simple_t)

    def test_extend_empty(self):
        """Extend on the empty case should work"""
        self.simple_t.extend([])
        self.assertEqual(self.simple_t.Children[0].Name, 'i1')
        self.assertEqual(self.simple_t.Children[1].Name, 'i2')
        self.assertEqual(len(self.simple_t.Children), 2)

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

    def test_pop(self):
        """Pop off a node"""
        second_tree = TreeNode.from_newick("(x1,y1)z1;")
        third_tree = TreeNode.from_newick("(x2,y2)z2;")
        self.simple_t.extend([second_tree, third_tree])

        i1 = self.simple_t.pop(0)
        z2 = self.simple_t.pop()

        self.assertEqual(i1.Name, 'i1')
        self.assertEqual(z2.Name, 'z2')
        self.assertEqual(i1.Children[0].Name, 'a')
        self.assertEqual(i1.Children[1].Name, 'b')
        self.assertEqual(z2.Children[0].Name, 'x2')
        self.assertEqual(z2.Children[1].Name, 'y2')

        self.assertEqual(self.simple_t.Children[0].Name, 'i2')
        self.assertEqual(self.simple_t.Children[1].Name, 'z1')
        self.assertEqual(len(self.simple_t.Children), 2)

    def test_remove(self):
        """Remove nodes"""
        self.assertTrue(self.simple_t.remove(self.simple_t.Children[0]))
        self.assertEqual(len(self.simple_t.Children), 1)
        n = TreeNode()
        self.assertFalse(self.simple_t.remove(n))

    def test_adopt(self):
        """Adopt a node!"""
        n1 = TreeNode(Name='n1')
        n2 = TreeNode(Name='n2')
        n3 = TreeNode(Name='n3')

        self.simple_t._adopt(n1)
        self.simple_t.Children[-1]._adopt(n2)
        n2._adopt(n3)

        # adopt doesn't update .Children
        self.assertEqual(len(self.simple_t.Children), 2)

        self.assertIs(n1.Parent, self.simple_t)
        self.assertIs(n2.Parent, self.simple_t.Children[-1])
        self.assertIs(n3.Parent, n2)

    def test_remove_node(self):
        """Remove a node by index"""
        n = self.simple_t._remove_node(-1)
        self.assertEqual(n.Parent, None)
        self.assertEqual(len(self.simple_t.Children), 1)
        self.assertEqual(len(n.Children), 2)
        self.assertNotIn(n, self.simple_t.Children)

    def test_prune(self):
        """Collapse single descendent nodes"""
        # check the identity case
        cp = self.simple_t.copy()
        self.simple_t.prune()

        gen = zip(cp.traverse(include_self=True),
                  self.simple_t.traverse(include_self=True))

        for a,b in gen:
            self.assertIsNot(a, b)
            self.assertEqual(a.Name, b.Name)
            self.assertEqual(a.Length, b.Length)

        # create a single descendent by removing tip 'a'
        n = self.simple_t.Children[0]
        n.remove(n.Children[0])
        self.simple_t.prune()

        self.assertEqual(len(self.simple_t.Children), 2)
        self.assertEqual(self.simple_t.Children[0].Name, 'i2')
        self.assertEqual(self.simple_t.Children[1].Name, 'b')

    def test_subset(self):
        """subset should return set of leaves that descends from node"""
        t = self.simple_t
        self.assertEqual(t.subset(), frozenset('abcd'))
        c = t.Children[0]
        self.assertEqual(c.subset(), frozenset('ab'))
        leaf = c.Children[1]
        self.assertEqual(leaf.subset(), frozenset(''))

    def test_subsets(self):
        """subsets should return all subsets descending from a set"""
        t = self.simple_t
        self.assertEqual(t.subsets(), frozenset(
                    [frozenset('ab'), frozenset('cd')]))


    def test_is_tip(self):
        """see if we're a tip or not"""
        self.assertFalse(self.simple_t.is_tip())
        self.assertFalse(self.simple_t.Children[0].is_tip())
        self.assertTrue(self.simple_t.Children[0].Children[0].is_tip())

    def test_is_root(self):
        """see if we're at the root or not"""
        self.assertTrue(self.simple_t.is_root())
        self.assertFalse(self.simple_t.Children[0].is_root())
        self.assertFalse(self.simple_t.Children[0].Children[0].is_root())

    def test_root(self):
        """Get the root!"""
        root = self.simple_t
        self.assertIs(root, self.simple_t.root())
        self.assertIs(root, self.simple_t.Children[0].root())
        self.assertIs(root, self.simple_t.Children[1].Children[1].root())

    def test_ancestors(self):
        """Get all the ancestors"""
        exp = ['i1','root']
        obs = self.simple_t.Children[0].Children[0].ancestors()
        self.assertEqual([o.Name for o in obs], exp)

        exp = ['root']
        obs = self.simple_t.Children[0].ancestors()
        self.assertEqual([o.Name for o in obs], exp)

        exp = []
        obs = self.simple_t.ancestors()
        self.assertEqual([o.Name for o in obs], exp)

    def test_siblings(self):
        """Get the siblings"""
        exp = []
        obs = self.simple_t.siblings()
        self.assertEqual(obs, exp)

        exp = ['i2']
        obs = self.simple_t.Children[0].siblings()
        self.assertEqual([o.Name for o in obs], exp)

        exp = ['c']
        obs = self.simple_t.Children[1].Children[1].siblings()
        self.assertEqual([o.Name for o in obs], exp)

        self.simple_t.append(TreeNode(Name="foo"))
        self.simple_t.append(TreeNode(Name="bar"))
        exp = ['i1','foo','bar']
        obs = self.simple_t.Children[1].siblings()
        self.assertEqual([o.Name for o in obs], exp)

    def test_ascii_art(self):
        """Make some ascii trees"""
        # unlabeled internal node
        tr = TreeNode.from_newick("(B:0.2,(C:0.3,D:0.4):0.6)F;")
        obs = tr.ascii_art(show_internal=True, compact=False)
        exp = """          /-B\n-F-------|\n         |          /-C\n          \\--------|\n                    \\-D"""
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=True, compact=True)
        exp = """-F------- /-B\n          \-------- /-C\n                    \-D"""
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=False, compact=False)
        exp = """          /-B\n---------|\n         |          /-C\n          \\--------|\n                    \\-D"""
        self.assertEqual(obs, exp)

    def test_accumulate_to_ancestor(self):
        """Get the distance from a node to its ancestor"""
        t = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;")
        a = t.find('a')
        exp_to_root = 0.1 + 0.3
        obs_to_root = a._accumulate_to_ancestor(t)
        self.assertEqual(obs_to_root, exp_to_root)

    def test_distance(self):
        """Get the distance between two nodes"""
        t = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;")
        tips = sorted([n for n in t.tips()], key=lambda x: x.Name)

        nptest.assert_almost_equal(tips[0].distance(tips[0]), 0.0)
        nptest.assert_almost_equal(tips[0].distance(tips[1]), 0.3)
        nptest.assert_almost_equal(tips[0].distance(tips[2]), 1.3)
        with self.assertRaises(NoLengthError):
            _ = tips[0].distance(tips[3])

        nptest.assert_almost_equal(tips[1].distance(tips[0]), 0.3)
        nptest.assert_almost_equal(tips[1].distance(tips[1]), 0.0)
        nptest.assert_almost_equal(tips[1].distance(tips[2]), 1.4)
        with self.assertRaises(NoLengthError):
            _ = tips[1].distance(tips[3])

        self.assertEqual(tips[2].distance(tips[0]), 1.3)
        self.assertEqual(tips[2].distance(tips[1]), 1.4)
        self.assertEqual(tips[2].distance(tips[2]), 0.0)
        with self.assertRaises(NoLengthError):
            _ = tips[2].distance(tips[3])

    def test_lowest_common_ancestor(self):
        """TreeNode lowestCommonAncestor should return LCA for set of tips"""
        t1 = TreeNode.from_newick("((a,(b,c)d)e,f,(g,h)i)j;")
        t2 = t1.copy()
        t3 = t1.copy()
        t4 = t1.copy()
        input1 = ['a'] # return self
        input2 = ['a','b'] # return e
        input3 = ['b','c'] # return d
        input4 = ['a','h','g'] # return j
        exp1 = t1.find('a')
        exp2 = t2.find('e')
        exp3 = t3.find('d')
        exp4 = t4
        obs1 = t1.lowest_common_ancestor(input1)
        obs2 = t2.lowest_common_ancestor(input2)
        obs3 = t3.lowest_common_ancestor(input3)
        obs4 = t4.lowest_common_ancestor(input4)
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)
        self.assertEqual(obs4, exp4)

        # verify multiple calls work
        t_mul = t1.copy()
        exp_1 = t_mul.find('d')
        exp_2 = t_mul.find('i')
        obs_1 = t_mul.lowest_common_ancestor(['b','c'])
        obs_2 = t_mul.lowest_common_ancestor(['g','h'])
        self.assertEqual(obs_1, exp_1)
        self.assertEqual(obs_2, exp_2)

    def test_get_max_distance(self):
        """get_max_distance should get max tip distance across tree"""
        tree = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;")
        dist, nodes = tree.get_max_distance()
        nptest.assert_almost_equal(dist, 1.6)
        self.assertEqual(sorted([n.Name for n in nodes]), ['b','e'])

    def test_set_max_distance(self):
        """set_max_distance sets MaxDistTips across tree"""
        tree = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;")
        tree._set_max_distance()
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 1.6)
        self.assertEqual(sorted([tip_a[1].Name, tip_b[1].Name]), ['b','e'])

    def test_compare_tip_distances(self):
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        t2 = TreeNode.from_newick('(((H:1,G:1,O:1):2,R:3):1,X:4);')
        obs = t.compare_tip_distances(t2)
        #note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = correlation_t(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1-r)/2)

    def test_compare_tip_distances_sample(self):
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        t2 = TreeNode.from_newick('(((H:1,G:1,O:1):2,R:3):1,X:4);')
        obs = t.compare_tip_distances(t2, sample=3, shuffle_f=sorted)
        #note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = correlation_t(m1.flat, m2.flat)[0]
        self.assertEqual(obs, (1-r)/2)

        # 4 common taxa, still picking H, G, R
        s = '((H:1,G:1):2,(R:0.5,M:0.7,Q:5):3);'
        t = TreeNode.from_newick(s, TreeNode)
        s3 = '(((H:1,G:1,O:1):2,R:3,Q:10):1,X:4);'
        t3 = TreeNode.from_newick(s3, TreeNode)
        obs = t.compare_tip_distances(t3, sample=3, shuffle_f=sorted)

    def test_tip_tip_distances_endpoints(self):
        """Test getting specifc tip distances  with tipToTipDistances"""
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        nodes = [t.find('H'), t.find('G'), t.find('M')]
        names = ['H', 'G', 'M']
        exp = np.array([[0, 2.0, 6.7], [2.0, 0, 6.7], [6.7, 6.7, 0.0]])
        exp_order = nodes

        obs, obs_order = t.tip_tip_distances(endpoints=names)
        nptest.assert_almost_equal(obs, exp)
        self.assertEqual(obs_order, exp_order)

        obs, obs_order = t.tip_tip_distances(endpoints=nodes)
        nptest.assert_almost_equal(obs, exp)
        self.assertEqual(obs_order, exp_order)

    def test_neighbors(self):
        """Get neighbors of a node"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        exp = t.Children
        obs = t.neighbors()
        self.assertEqual(obs, exp)

        exp = t.Children[0].Children + [t]
        obs = t.Children[0].neighbors()
        self.assertEqual(obs, exp)

        exp = [t.Children[0].Children[0]] + [t]
        obs = t.Children[0].neighbors(ignore=t.Children[0].Children[1])
        self.assertEqual(obs, exp)

        exp = [t.Children[0]]
        obs = t.Children[0].Children[0].neighbors()
        self.assertEqual(obs, exp)

    def test_has_children(self):
        """Test if has children"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        self.assertTrue(t.has_children())
        self.assertTrue(t.Children[0].has_children())
        self.assertTrue(t.Children[1].has_children())
        self.assertFalse(t.Children[0].Children[0].has_children())
        self.assertFalse(t.Children[0].Children[1].has_children())
        self.assertFalse(t.Children[1].Children[0].has_children())
        self.assertFalse(t.Children[1].Children[1].has_children())

    def test_index_tree(self):
        """index_tree should produce correct index and node map"""
        #test for first tree: contains singleton outgroup
        t1 = TreeNode.from_newick('(((a,b),c),(d,e))')
        t2 = TreeNode.from_newick('(((a,b),(c,d)),(e,f))')
        t3 = TreeNode.from_newick('(((a,b,c),(d)),(e,f))')

        id_1, child_1 = t1.index_tree()
        nodes_1 = [n._leaf_index for n in t1.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_1, [0, 1, 2, 3, 6, 4, 5, 7, 8])
        self.assertEqual(child_1, [(2, 0, 1), (6, 2, 3), (7, 4, 5), (8, 6, 7)])

        #test for second tree: strictly bifurcating
        id_2, child_2 = t2.index_tree()
        nodes_2 = [n._leaf_index for n in t2.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_2, [0, 1, 4, 2, 3, 5, 8, 6, 7, 9, 10])
        self.assertEqual(child_2, [(4, 0, 1), (5, 2, 3), (8, 4, 5), (9, 6, 7),
                                   (10, 8, 9)])

        #test for third tree: contains trifurcation and single-child parent
        id_3, child_3 = t3.index_tree()
        nodes_3 = [n._leaf_index for n in t3.traverse(self_before=False,
                   self_after=True)]

        self.assertEqual(nodes_3, [0, 1, 2, 4, 3, 5, 8, 6, 7, 9, 10])
        self.assertEqual(child_3, [(4, 0, 2), (5, 3, 3), (8, 4, 5), (9, 6, 7),
                                   (10, 8, 9)])

    def test_root_at(self):
        """Form a new root"""
        t = TreeNode.from_newick("(((a,b)c,(d,e)f)g,h)i;")
        with self.assertRaises(TreeError):
            _ = t.root_at(t.find('h'))

        exp = "(a,b,((d,e)f,(h)g)c)root;"
        rooted = t.root_at('c')
        obs = str(rooted)
        self.assertEqual(obs, exp)

    def test_root_at_midpoint(self):
        """Root at the midpoint"""
        nodes, tree = self.TreeNode, self.TreeRoot
        tree1 = tree.copy()

        for n in tree1.traverse():
            n.Length = 1

        result = tree1.root_at_midpoint()
        self.assertEqual(result.distance(result.find('e')), 1.5)
        self.assertEqual(result.distance(result.find('g')), 2.5)
        exp_dist, exp_order = tree1.tip_tip_distances()
        obs_dist, obs_order = result.tip_tip_distances()
        nptest.assert_almost_equal(obs_dist, exp_dist)
        self.assertEqual([n.Name for n in obs_order],
                         [n.Name for n in exp_order])

    def test_compare_subsets(self):
        """compare_subsets should return the fraction of shared subsets"""
        t = TreeNode.from_newick('((H,G),(R,M));')
        t2 = TreeNode.from_newick('(((H,G),R),M);')
        t4 = TreeNode.from_newick('(((H,G),(O,R)),X);')

        result = t.compare_subsets(t)
        self.assertEqual(result, 0)

        result = t2.compare_subsets(t2)
        self.assertEqual(result, 0)

        result = t.compare_subsets(t2)
        self.assertEqual(result, 0.5)

        result = t.compare_subsets(t4)
        self.assertEqual(result, 1-2./5)

        result = t.compare_subsets(t4, exclude_absent_taxa=True)
        self.assertEqual(result, 1-2./3)

        result = t.compare_subsets(self.TreeRoot, exclude_absent_taxa=True)
        self.assertEqual(result, 1)

        result = t.compare_subsets(self.TreeRoot)
        self.assertEqual(result, 1)

    def test_assign_ids(self):
        """Assign IDs to the tree"""
        t1 = TreeNode.from_newick("(((a,b),c),(e,f),(g));")
        t2 = TreeNode.from_newick("(((a,b),c),(e,f),(g));")
        t3 = TreeNode.from_newick("((g),(e,f),(c,(a,b)));")
        t1_copy = t1.copy()

        t1.assign_ids()
        t2.assign_ids()
        t3.assign_ids()
        t1_copy.assign_ids()

        self.assertEqual([(n.Name, n.Id) for n in t1.traverse()],
                         [(n.Name, n.Id) for n in t2.traverse()])
        self.assertEqual([(n.Name, n.Id) for n in t1.traverse()],
                         [(n.Name, n.Id) for n in t1_copy.traverse()])
        self.assertNotEqual([(n.Name, n.Id) for n in t1.traverse()],
                            [(n.Name, n.Id) for n in t3.traverse()])


    def test_unrooted_deepcopy(self):
        """Do an unrooted_copy"""
        t = TreeNode.from_newick("((a,(b,c)d)e,(f,g)h)i;")
        exp = "(b,c,(a,((f,g)h)e)d)root;"
        obs = t.find('d').unrooted_deepcopy()
        self.assertEqual(str(obs), exp)

        t_ids = {id(n) for n in t.traverse()}
        obs_ids = {id(n) for n in obs.traverse()}

        self.assertEqual(t_ids.intersection(obs_ids), set())

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
        obs = list(_dnd_tokenizer(sample))
        self.assertEqual(len(obs), len(exp))
        for i, j in zip(obs, exp):
            self.assertEqual(i, j)
        #try it all in one go
        self.assertEqual(list(_dnd_tokenizer(sample)), exp)

    def test_nonames(self):
        """DndTokenizer should work as expected on trees with no names"""
        exp = ['(','(',',',')',',','(',',',')',')',';']
        obs = list(_dnd_tokenizer(no_names))
        self.assertEqual(obs, exp)

    def test_missing_tip_name(self):
        """DndTokenizer should work as expected on trees with a missing name"""
        exp = ['(','(','a',',','b',')',',','(','c',',',')',')',';']
        obs = list(_dnd_tokenizer(missing_tip_name))
        self.assertEqual(obs, exp)

    def test_minimal(self):
        """DndTokenizer should work as expected a minimal tree without names"""
        exp = ['(',')',';']
        obs = list(_dnd_tokenizer(minimal))
        self.assertEqual(obs, exp)

class DndParserTests(TestCase):
    """Tests of the DndParser factory function."""
    def test_nonames(self):
        """DndParser should produce the correct tree when there are no names"""
        obs = TreeNode.from_newick(no_names)
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
        obs = TreeNode.from_newick(minimal)
        exp = TreeNode()
        exp.append(TreeNode())
        self.assertEqual(str(obs), str(exp))

    def test_missing_tip_name(self):
        """DndParser should produce the correct tree when missing a name"""
        obs = TreeNode.from_newick(missing_tip_name)
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
        t = TreeNode.from_newick(single)
        self.assertEqual(len(t), 1)
        child = t[0]
        self.assertEqual(child.Name, 'abc')
        self.assertEqual(child.Length, 3)
        self.assertEqual(str(t), '(abc:3.0);')

    def test_gdouble(self):
        """DndParser should produce a double-child TreeNode from data"""
        t = TreeNode.from_newick(double)
        self.assertEqual(len(t), 2)
        self.assertEqual(str(t), '(abc:3.0,def:4.0);')

    def test_gonenest(self):
        """DndParser should work correctly with nested data"""
        t = TreeNode.from_newick(onenest)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  #first child is terminal
        self.assertEqual(len(t[1]), 2)  #second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0):6.0);')

    def test_gnodedata(self):
        """DndParser should assign Name to internal nodes correctly"""
        t = TreeNode.from_newick(nodedata)
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
        t = TreeNode.from_newick(sample)
        self.assertEqual(str(t), '((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);')
        tdata = TreeNode.from_newick(node_data_sample, unescape_name=True)
        self.assertEqual(str(tdata), "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)C:0.09853);")

    def test_gbad(self):
        """DndParser should fail if parens unbalanced"""
        left = '((abc:3)'
        right = '(abc:3))'
        self.assertRaises(RecordError, TreeNode.from_newick, left)
        self.assertRaises(RecordError, TreeNode.from_newick, right)

    def test_DndParser(self):
        """DndParser tests"""
        t_str = "(A_a,(B:1.0,C),'D_e':0.5)E;"
        tree_unesc = TreeNode.from_newick(t_str, unescape_name=True)
        tree_esc = TreeNode.from_newick(t_str, unescape_name=False)

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
        obs = TreeNode.from_newick(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True), \
                         tree_esc.to_newick(with_distances=True))
        reload_test = tree_unesc.to_newick(with_distances=True, \
                                           escape_name=False)
        obs = TreeNode.from_newick(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True), \
                         tree_unesc.to_newick(with_distances=True))

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

if __name__ == '__main__':
    main()
