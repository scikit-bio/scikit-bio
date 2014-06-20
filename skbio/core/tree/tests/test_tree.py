#!/usr/bin/env python

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
import numpy.testing as nptest
from scipy.stats import pearsonr
from future.utils.six import StringIO

from skbio.core.tree import TreeNode
from skbio.core.tree.tree import _dnd_tokenizer
from skbio.core.distance import DistanceMatrix
from skbio.core.exception import (DuplicateNodeError, NoLengthError, TreeError,
                                  RecordError, MissingNodeError, NoParentError)


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

    def test_count(self):
        """Get node counts"""
        exp = 7
        obs = self.simple_t.count()
        self.assertEqual(obs, exp)

        exp = 4
        obs = self.simple_t.count(tips=True)
        self.assertEqual(obs, exp)

    def test_copy(self):
        """copy a tree"""
        self.simple_t.children[0].length = 1.2
        self.simple_t.children[1].children[0].length = 0.5
        cp = self.simple_t.copy()
        gen = zip(cp.traverse(include_self=True),
                  self.simple_t.traverse(include_self=True))

        for a, b in gen:
            self.assertIsNot(a, b)
            self.assertEqual(a.name, b.name)
            self.assertEqual(a.length, b.length)

    def test_append(self):
        """Append a node to a tree"""
        second_tree = TreeNode.from_newick("(x,y)z;")
        self.simple_t.append(second_tree)

        self.assertEqual(self.simple_t.children[0].name, 'i1')
        self.assertEqual(self.simple_t.children[1].name, 'i2')
        self.assertEqual(self.simple_t.children[2].name, 'z')
        self.assertEqual(len(self.simple_t.children), 3)
        self.assertEqual(self.simple_t.children[2].children[0].name, 'x')
        self.assertEqual(self.simple_t.children[2].children[1].name, 'y')
        self.assertEqual(second_tree.parent, self.simple_t)

    def test_extend(self):
        """Extend a few nodes"""
        second_tree = TreeNode.from_newick("(x1,y1)z1;")
        third_tree = TreeNode.from_newick("(x2,y2)z2;")
        self.simple_t.extend([second_tree, third_tree])

        self.assertEqual(self.simple_t.children[0].name, 'i1')
        self.assertEqual(self.simple_t.children[1].name, 'i2')
        self.assertEqual(self.simple_t.children[2].name, 'z1')
        self.assertEqual(self.simple_t.children[3].name, 'z2')
        self.assertEqual(len(self.simple_t.children), 4)
        self.assertEqual(self.simple_t.children[2].children[0].name, 'x1')
        self.assertEqual(self.simple_t.children[2].children[1].name, 'y1')
        self.assertEqual(self.simple_t.children[3].children[0].name, 'x2')
        self.assertEqual(self.simple_t.children[3].children[1].name, 'y2')
        self.assertIs(second_tree.parent, self.simple_t)
        self.assertIs(third_tree.parent, self.simple_t)

    def test_extend_empty(self):
        """Extend on the empty case should work"""
        self.simple_t.extend([])
        self.assertEqual(self.simple_t.children[0].name, 'i1')
        self.assertEqual(self.simple_t.children[1].name, 'i2')
        self.assertEqual(len(self.simple_t.children), 2)

    def test_iter(self):
        """iter wraps children"""
        exp = ['i1', 'i2']
        obs = [n.name for n in self.simple_t]
        self.assertEqual(obs, exp)

    def test_gops(self):
        """Basic TreeNode operations should work as expected"""
        p = TreeNode()
        self.assertEqual(str(p), ';')
        p.name = 'abc'
        self.assertEqual(str(p), 'abc;')
        p.length = 3
        self.assertEqual(str(p), 'abc:3;')  # don't suppress branch from root
        q = TreeNode()
        p.append(q)
        self.assertEqual(str(p), '()abc:3;')
        r = TreeNode()
        q.append(r)
        self.assertEqual(str(p), '(())abc:3;')
        r.name = 'xyz'
        self.assertEqual(str(p), '((xyz))abc:3;')
        q.length = 2
        self.assertEqual(str(p), '((xyz):2)abc:3;')

    def test_pop(self):
        """Pop off a node"""
        second_tree = TreeNode.from_newick("(x1,y1)z1;")
        third_tree = TreeNode.from_newick("(x2,y2)z2;")
        self.simple_t.extend([second_tree, third_tree])

        i1 = self.simple_t.pop(0)
        z2 = self.simple_t.pop()

        self.assertEqual(i1.name, 'i1')
        self.assertEqual(z2.name, 'z2')
        self.assertEqual(i1.children[0].name, 'a')
        self.assertEqual(i1.children[1].name, 'b')
        self.assertEqual(z2.children[0].name, 'x2')
        self.assertEqual(z2.children[1].name, 'y2')

        self.assertEqual(self.simple_t.children[0].name, 'i2')
        self.assertEqual(self.simple_t.children[1].name, 'z1')
        self.assertEqual(len(self.simple_t.children), 2)

    def test_remove(self):
        """Remove nodes"""
        self.assertTrue(self.simple_t.remove(self.simple_t.children[0]))
        self.assertEqual(len(self.simple_t.children), 1)
        n = TreeNode()
        self.assertFalse(self.simple_t.remove(n))

    def test_remove_deleted(self):
        """Remove nodes by function"""
        f = lambda node: node.name in ['b', 'd']
        self.simple_t.remove_deleted(f)
        exp = "((a)i1,(c)i2)root;"
        obs = self.simple_t.to_newick()
        self.assertEqual(obs, exp)

    def test_adopt(self):
        """Adopt a node!"""
        n1 = TreeNode(name='n1')
        n2 = TreeNode(name='n2')
        n3 = TreeNode(name='n3')

        self.simple_t._adopt(n1)
        self.simple_t.children[-1]._adopt(n2)
        n2._adopt(n3)

        # adopt doesn't update .children
        self.assertEqual(len(self.simple_t.children), 2)

        self.assertIs(n1.parent, self.simple_t)
        self.assertIs(n2.parent, self.simple_t.children[-1])
        self.assertIs(n3.parent, n2)

    def test_remove_node(self):
        """Remove a node by index"""
        n = self.simple_t._remove_node(-1)
        self.assertEqual(n.parent, None)
        self.assertEqual(len(self.simple_t.children), 1)
        self.assertEqual(len(n.children), 2)
        self.assertNotIn(n, self.simple_t.children)

    def test_prune(self):
        """Collapse single descendent nodes"""
        # check the identity case
        cp = self.simple_t.copy()
        self.simple_t.prune()

        gen = zip(cp.traverse(include_self=True),
                  self.simple_t.traverse(include_self=True))

        for a, b in gen:
            self.assertIsNot(a, b)
            self.assertEqual(a.name, b.name)
            self.assertEqual(a.length, b.length)

        # create a single descendent by removing tip 'a'
        n = self.simple_t.children[0]
        n.remove(n.children[0])
        self.simple_t.prune()

        self.assertEqual(len(self.simple_t.children), 2)
        self.assertEqual(self.simple_t.children[0].name, 'i2')
        self.assertEqual(self.simple_t.children[1].name, 'b')

    def test_prune_length(self):
        """Collapse single descendent nodes"""
        # check the identity case
        cp = self.simple_t.copy()
        self.simple_t.prune()

        gen = zip(cp.traverse(include_self=True),
                  self.simple_t.traverse(include_self=True))

        for a, b in gen:
            self.assertIsNot(a, b)
            self.assertEqual(a.name, b.name)
            self.assertEqual(a.length, b.length)

        for n in self.simple_t.traverse():
            n.length = 1.0

        # create a single descendent by removing tip 'a'
        n = self.simple_t.children[0]
        n.remove(n.children[0])
        self.simple_t.prune()

        self.assertEqual(len(self.simple_t.children), 2)
        self.assertEqual(self.simple_t.children[0].name, 'i2')
        self.assertEqual(self.simple_t.children[1].name, 'b')
        self.assertEqual(self.simple_t.children[1].length, 2.0)

    def test_subset(self):
        """subset should return set of leaves that descends from node"""
        t = self.simple_t
        self.assertEqual(t.subset(), frozenset('abcd'))
        c = t.children[0]
        self.assertEqual(c.subset(), frozenset('ab'))
        leaf = c.children[1]
        self.assertEqual(leaf.subset(), frozenset(''))

    def test_subsets(self):
        """subsets should return all subsets descending from a set"""
        t = self.simple_t
        self.assertEqual(t.subsets(), frozenset(
            [frozenset('ab'), frozenset('cd')]))

    def test_is_tip(self):
        """see if we're a tip or not"""
        self.assertFalse(self.simple_t.is_tip())
        self.assertFalse(self.simple_t.children[0].is_tip())
        self.assertTrue(self.simple_t.children[0].children[0].is_tip())

    def test_is_root(self):
        """see if we're at the root or not"""
        self.assertTrue(self.simple_t.is_root())
        self.assertFalse(self.simple_t.children[0].is_root())
        self.assertFalse(self.simple_t.children[0].children[0].is_root())

    def test_root(self):
        """Get the root!"""
        root = self.simple_t
        self.assertIs(root, self.simple_t.root())
        self.assertIs(root, self.simple_t.children[0].root())
        self.assertIs(root, self.simple_t.children[1].children[1].root())

    def test_create_caches_duplicate_tip_names(self):
        with self.assertRaises(DuplicateNodeError):
            TreeNode.from_newick('(a, a)').create_caches()

    def test_find(self):
        """Find a node in a tree"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        exp = t.children[0]
        obs = t.find('c')
        self.assertEqual(obs, exp)

        exp = t.children[0].children[1]
        obs = t.find('b')
        self.assertEqual(obs, exp)

        with self.assertRaises(MissingNodeError):
            t.find('does not exist')

    def test_find_cache_bug(self):
        """First implementation did not force the cache to be at the root"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f,(g,h)f);")
        exp_tip_cache_keys = set(['a', 'b', 'd', 'e', 'g', 'h'])
        exp_non_tip_cache_keys = set(['c', 'f'])
        tip_a = t.children[0].children[0]
        tip_a.create_caches()
        self.assertEqual(tip_a._tip_cache, {})
        self.assertEqual(set(t._tip_cache), exp_tip_cache_keys)
        self.assertEqual(set(t._non_tip_cache), exp_non_tip_cache_keys)
        self.assertEqual(t._non_tip_cache['f'], [t.children[1], t.children[2]])

    def test_find_by_id(self):
        """Find a node by id"""
        t1 = TreeNode.from_newick("((,),(,,));")
        t2 = TreeNode.from_newick("((,),(,,));")

        exp = t1.children[1]
        obs = t1.find_by_id(6)  # right inner node with 3 children
        self.assertEqual(obs, exp)

        exp = t2.children[1]
        obs = t2.find_by_id(6)  # right inner node with 3 children
        self.assertEqual(obs, exp)

        with self.assertRaises(MissingNodeError):
            t1.find_by_id(100)

    def test_find_by_func(self):
        """Find nodes by a function"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        func = lambda x: x.parent == t.find('c')
        exp = ['a', 'b']
        obs = [n.name for n in t.find_by_func(func)]
        self.assertEqual(obs, exp)

    def test_ancestors(self):
        """Get all the ancestors"""
        exp = ['i1', 'root']
        obs = self.simple_t.children[0].children[0].ancestors()
        self.assertEqual([o.name for o in obs], exp)

        exp = ['root']
        obs = self.simple_t.children[0].ancestors()
        self.assertEqual([o.name for o in obs], exp)

        exp = []
        obs = self.simple_t.ancestors()
        self.assertEqual([o.name for o in obs], exp)

    def test_siblings(self):
        """Get the siblings"""
        exp = []
        obs = self.simple_t.siblings()
        self.assertEqual(obs, exp)

        exp = ['i2']
        obs = self.simple_t.children[0].siblings()
        self.assertEqual([o.name for o in obs], exp)

        exp = ['c']
        obs = self.simple_t.children[1].children[1].siblings()
        self.assertEqual([o.name for o in obs], exp)

        self.simple_t.append(TreeNode(name="foo"))
        self.simple_t.append(TreeNode(name="bar"))
        exp = ['i1', 'foo', 'bar']
        obs = self.simple_t.children[1].siblings()
        self.assertEqual([o.name for o in obs], exp)

    def test_ascii_art(self):
        """Make some ascii trees"""
        # unlabeled internal node
        tr = TreeNode.from_newick("(B:0.2,(C:0.3,D:0.4):0.6)F;")
        obs = tr.ascii_art(show_internal=True, compact=False)
        exp = "          /-B\n-F-------|\n         |          /-C\n         "\
              " \\--------|\n                    \\-D"
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=True, compact=True)
        exp = "-F------- /-B\n          \-------- /-C\n                    \-D"
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=False, compact=False)
        exp = "          /-B\n---------|\n         |          /-C\n         "\
              " \\--------|\n                    \\-D"
        self.assertEqual(obs, exp)

    def test_accumulate_to_ancestor(self):
        """Get the distance from a node to its ancestor"""
        t = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;")
        a = t.find('a')
        b = t.find('b')
        exp_to_root = 0.1 + 0.3
        obs_to_root = a.accumulate_to_ancestor(t)
        self.assertEqual(obs_to_root, exp_to_root)

        with self.assertRaises(NoParentError):
            a.accumulate_to_ancestor(b)

    def test_distance(self):
        """Get the distance between two nodes"""
        t = TreeNode.from_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;")
        tips = sorted([n for n in t.tips()], key=lambda x: x.name)

        nptest.assert_almost_equal(tips[0].distance(tips[0]), 0.0)
        nptest.assert_almost_equal(tips[0].distance(tips[1]), 0.3)
        nptest.assert_almost_equal(tips[0].distance(tips[2]), 1.3)
        with self.assertRaises(NoLengthError):
            tips[0].distance(tips[3])

        nptest.assert_almost_equal(tips[1].distance(tips[0]), 0.3)
        nptest.assert_almost_equal(tips[1].distance(tips[1]), 0.0)
        nptest.assert_almost_equal(tips[1].distance(tips[2]), 1.4)
        with self.assertRaises(NoLengthError):
            tips[1].distance(tips[3])

        self.assertEqual(tips[2].distance(tips[0]), 1.3)
        self.assertEqual(tips[2].distance(tips[1]), 1.4)
        self.assertEqual(tips[2].distance(tips[2]), 0.0)
        with self.assertRaises(NoLengthError):
            tips[2].distance(tips[3])

    def test_lowest_common_ancestor(self):
        """TreeNode lowestCommonAncestor should return LCA for set of tips"""
        t1 = TreeNode.from_newick("((a,(b,c)d)e,f,(g,h)i)j;")
        t2 = t1.copy()
        t3 = t1.copy()
        t4 = t1.copy()
        input1 = ['a']  # return self
        input2 = ['a', 'b']  # return e
        input3 = ['b', 'c']  # return d
        input4 = ['a', 'h', 'g']  # return j
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
        obs_1 = t_mul.lowest_common_ancestor(['b', 'c'])
        obs_2 = t_mul.lowest_common_ancestor(['g', 'h'])
        self.assertEqual(obs_1, exp_1)
        self.assertEqual(obs_2, exp_2)

        # empty case
        with self.assertRaises(ValueError):
            t1.lowest_common_ancestor([])

    def test_get_max_distance(self):
        """get_max_distance should get max tip distance across tree"""
        tree = TreeNode.from_newick(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;")
        dist, nodes = tree.get_max_distance()
        nptest.assert_almost_equal(dist, 1.6)
        self.assertEqual(sorted([n.name for n in nodes]), ['b', 'e'])

    def test_set_max_distance(self):
        """set_max_distance sets MaxDistTips across tree"""
        tree = TreeNode.from_newick(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;")
        tree._set_max_distance()
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 1.6)
        self.assertEqual(sorted([tip_a[1].name, tip_b[1].name]), ['b', 'e'])

    def test_shear(self):
        """Shear the nodes"""
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        obs = t.shear(['G', 'M']).to_newick(with_distances=True)
        exp = '(G:3.0,M:3.7);'
        self.assertEqual(obs, exp)

    def test_compare_tip_distances(self):
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        t2 = TreeNode.from_newick('(((H:1,G:1,O:1):2,R:3):1,X:4);')
        obs = t.compare_tip_distances(t2)
        # note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

    def test_compare_tip_distances_sample(self):
        t = TreeNode.from_newick('((H:1,G:1):2,(R:0.5,M:0.7):3);')
        t2 = TreeNode.from_newick('(((H:1,G:1,O:1):2,R:3):1,X:4);')
        obs = t.compare_tip_distances(t2, sample=3, shuffle_f=sorted)
        # note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

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
        exp = DistanceMatrix(np.array([[0, 2.0, 6.7],
                                       [2.0, 0, 6.7],
                                       [6.7, 6.7, 0.0]]), ['H', 'G', 'M'])

        obs = t.tip_tip_distances(endpoints=names)
        self.assertEqual(obs, exp)

        obs = t.tip_tip_distances(endpoints=nodes)
        self.assertEqual(obs, exp)

    def test_neighbors(self):
        """Get neighbors of a node"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        exp = t.children
        obs = t.neighbors()
        self.assertEqual(obs, exp)

        exp = t.children[0].children + [t]
        obs = t.children[0].neighbors()
        self.assertEqual(obs, exp)

        exp = [t.children[0].children[0]] + [t]
        obs = t.children[0].neighbors(ignore=t.children[0].children[1])
        self.assertEqual(obs, exp)

        exp = [t.children[0]]
        obs = t.children[0].children[0].neighbors()
        self.assertEqual(obs, exp)

    def test_has_children(self):
        """Test if has children"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f);")
        self.assertTrue(t.has_children())
        self.assertTrue(t.children[0].has_children())
        self.assertTrue(t.children[1].has_children())
        self.assertFalse(t.children[0].children[0].has_children())
        self.assertFalse(t.children[0].children[1].has_children())
        self.assertFalse(t.children[1].children[0].has_children())
        self.assertFalse(t.children[1].children[1].has_children())

    def test_tips(self):
        """Tip traversal of tree"""
        exp = ['a', 'b', 'c', 'd']
        obs = [n.name for n in self.simple_t.tips()]
        self.assertEqual(obs, exp)
        obs2 = [n.name for n in self.simple_t.traverse(False, False)]
        self.assertEqual(obs2, exp)

    def test_pre_and_postorder(self):
        """Pre and post order traversal of the tree"""
        exp = ['root', 'i1', 'a', 'b', 'i1', 'i2', 'c', 'd', 'i2', 'root']
        obs = [n.name for n in self.simple_t.pre_and_postorder()]
        self.assertEqual(obs, exp)
        obs2 = [n.name for n in self.simple_t.traverse(True, True)]
        self.assertEqual(obs2, exp)

    def test_pre_and_postorder_no_children(self):
        t = TreeNode('brofist')

        # include self
        exp = ['brofist']
        obs = [n.name for n in t.pre_and_postorder()]
        self.assertEqual(obs, exp)

        # do not include self
        obs = list(t.pre_and_postorder(include_self=False))
        self.assertEqual(obs, [])

    def test_levelorder(self):
        """Test level order traversal of the tree"""
        exp = ['root', 'i1', 'i2', 'a', 'b', 'c', 'd']
        obs = [n.name for n in self.simple_t.levelorder()]
        self.assertEqual(obs, exp)

    def test_index_tree(self):
        """index_tree should produce correct index and node map"""
        # test for first tree: contains singleton outgroup
        t1 = TreeNode.from_newick('(((a,b),c),(d,e))')
        t2 = TreeNode.from_newick('(((a,b),(c,d)),(e,f))')
        t3 = TreeNode.from_newick('(((a,b,c),(d)),(e,f))')

        id_1, child_1 = t1.index_tree()
        nodes_1 = [n.id for n in t1.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_1, [0, 1, 2, 3, 6, 4, 5, 7, 8])
        self.assertEqual(child_1, [(2, 0, 1), (6, 2, 3), (7, 4, 5), (8, 6, 7)])

        # test for second tree: strictly bifurcating
        id_2, child_2 = t2.index_tree()
        nodes_2 = [n.id for n in t2.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_2, [0, 1, 4, 2, 3, 5, 8, 6, 7, 9, 10])
        self.assertEqual(child_2, [(4, 0, 1), (5, 2, 3), (8, 4, 5), (9, 6, 7),
                                   (10, 8, 9)])

        # test for third tree: contains trifurcation and single-child parent
        id_3, child_3 = t3.index_tree()
        nodes_3 = [n.id for n in t3.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_3, [0, 1, 2, 4, 3, 5, 8, 6, 7, 9, 10])
        self.assertEqual(child_3, [(4, 0, 2), (5, 3, 3), (8, 4, 5), (9, 6, 7),
                                   (10, 8, 9)])

    def test_root_at(self):
        """Form a new root"""
        t = TreeNode.from_newick("(((a,b)c,(d,e)f)g,h)i;")
        with self.assertRaises(TreeError):
            t.root_at(t.find('h'))

        exp = "(a,b,((d,e)f,(h)g)c)root;"
        rooted = t.root_at('c')
        obs = str(rooted)
        self.assertEqual(obs, exp)

    def test_root_at_midpoint(self):
        """Root at the midpoint"""
        tree1 = self.TreeRoot
        for n in tree1.traverse():
            n.length = 1

        result = tree1.root_at_midpoint()
        self.assertEqual(result.distance(result.find('e')), 1.5)
        self.assertEqual(result.distance(result.find('g')), 2.5)
        exp_dist = tree1.tip_tip_distances()
        obs_dist = result.tip_tip_distances()
        self.assertEqual(obs_dist, exp_dist)

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
        self.assertEqual(result, 1 - 2. / 5)

        result = t.compare_subsets(t4, exclude_absent_taxa=True)
        self.assertEqual(result, 1 - 2. / 3)

        result = t.compare_subsets(self.TreeRoot, exclude_absent_taxa=True)
        self.assertEqual(result, 1)

        result = t.compare_subsets(self.TreeRoot)
        self.assertEqual(result, 1)

    def test_compare_rfd(self):
        """compare_rfd should return the Robinson Foulds distance"""
        t = TreeNode.from_newick('((H,G),(R,M));')
        t2 = TreeNode.from_newick('(((H,G),R),M);')
        t4 = TreeNode.from_newick('(((H,G),(O,R)),X);')

        obs = t.compare_rfd(t2)
        exp = 2.0
        self.assertEqual(obs, exp)

        self.assertEqual(t.compare_rfd(t2), t2.compare_rfd(t))

        obs = t.compare_rfd(t2, proportion=True)
        exp = 0.5
        self.assertEqual(obs, exp)

        with self.assertRaises(ValueError):
            t.compare_rfd(t4)

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

        self.assertEqual([(n.name, n.id) for n in t1.traverse()],
                         [(n.name, n.id) for n in t2.traverse()])
        self.assertEqual([(n.name, n.id) for n in t1.traverse()],
                         [(n.name, n.id) for n in t1_copy.traverse()])
        self.assertNotEqual([(n.name, n.id) for n in t1.traverse()],
                            [(n.name, n.id) for n in t3.traverse()])

    def test_assign_ids_index_tree(self):
        """assign_ids and index_tree should assign the same IDs"""
        t1 = TreeNode.from_newick('(((a,b),c),(d,e))')
        t2 = TreeNode.from_newick('(((a,b),(c,d)),(e,f))')
        t3 = TreeNode.from_newick('(((a,b,c),(d)),(e,f))')
        t1_copy = t1.copy()
        t2_copy = t2.copy()
        t3_copy = t3.copy()

        t1.assign_ids()
        t1_copy.index_tree()
        t2.assign_ids()
        t2_copy.index_tree()
        t3.assign_ids()
        t3_copy.index_tree()

        self.assertEqual([n.id for n in t1.traverse()],
                         [n.id for n in t1_copy.traverse()])
        self.assertEqual([n.id for n in t2.traverse()],
                         [n.id for n in t2_copy.traverse()])
        self.assertEqual([n.id for n in t3.traverse()],
                         [n.id for n in t3_copy.traverse()])

    def test_unrooted_deepcopy(self):
        """Do an unrooted_copy"""
        t = TreeNode.from_newick("((a,(b,c)d)e,(f,g)h)i;")
        exp = "(b,c,(a,((f,g)h)e)d)root;"
        obs = t.find('d').unrooted_deepcopy()
        self.assertEqual(str(obs), exp)

        t_ids = {id(n) for n in t.traverse()}
        obs_ids = {id(n) for n in obs.traverse()}

        self.assertEqual(t_ids.intersection(obs_ids), set())

    def test_descending_branch_length(self):
        """Calculate descending branch_length"""
        tr = TreeNode.from_newick("(((A:.1,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H"
                                  ":.4,I:.5)J:1.3)K;")
        tdbl = tr.descending_branch_length()
        sdbl = tr.descending_branch_length(['A', 'E'])
        nptest.assert_almost_equal(tdbl, 8.9)
        nptest.assert_almost_equal(sdbl, 2.2)
        self.assertRaises(ValueError, tr.descending_branch_length,
                          ['A', 'DNE'])
        self.assertRaises(ValueError, tr.descending_branch_length, ['A', 'C'])

        tr = TreeNode.from_newick("(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4"
                                  ",I:.5)J:1.3)K;")
        tdbl = tr.descending_branch_length()
        nptest.assert_almost_equal(tdbl, 8.8)

        tr = TreeNode.from_newick("(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:"
                                  ".5)J:1.3)K;")
        tdbl = tr.descending_branch_length()
        nptest.assert_almost_equal(tdbl, 7.9)

        tr = TreeNode.from_newick("(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:"
                                  ".5)J:1.3)K;")
        tdbl = tr.descending_branch_length(['A', 'D', 'E'])
        nptest.assert_almost_equal(tdbl, 2.1)

        tr = TreeNode.from_newick("(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:."
                                  "4,I:.5)J:1.3)K;")
        tdbl = tr.descending_branch_length(['I', 'D', 'E'])
        nptest.assert_almost_equal(tdbl, 6.6)

        # test with a situation where we have unnamed internal nodes
        tr = TreeNode.from_newick("(((A,B:1.2):.6,(D:.9,E:.6)F):2.4,(H:.4,I:"
                                  ".5)J:1.3);")
        tdbl = tr.descending_branch_length()
        nptest.assert_almost_equal(tdbl, 7.9)

    def test_to_array(self):
        """Convert a tree to arrays"""
        t = TreeNode.from_newick(
            '(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10)')
        id_index, child_index = t.index_tree()
        arrayed = t.to_array()

        self.assertEqual(id_index, arrayed['id_index'])
        self.assertEqual(child_index, arrayed['child_index'])

        exp = np.array([1, 2, 3, 5, 4, 6, 8, 9, 7, 10, np.nan])
        obs = arrayed['length']
        nptest.assert_equal(obs, exp)

        exp = np.array(['a', 'b', 'c', 'd', 'x',
                        'y', 'e', 'f', 'z', 'z', None])
        obs = arrayed['name']
        nptest.assert_equal(obs, exp)

        exp = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        obs = arrayed['id']
        nptest.assert_equal(obs, exp)

    def test_from_file(self):
        """Parse a tree from a file"""
        t_io = StringIO("((a,b)c,(d,e)f)g;")
        t = TreeNode.from_file(t_io)
        self.assertEqual(list('abcdefg'), [n.name for n in t.postorder()])

    def test_linkage_matrix(self):
        # Ensure matches: http://www.southampton.ac.uk/~re1u06/teaching/upgma/
        id_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        linkage = np.asarray([[1.0,  5.0,  1.0,  2.0],
                              [0.0,  3.0,  8.0,  2.0],
                              [6.0,  7.0, 12.5,  3.0],
                              [8.0,  9.0, 16.5,  5.0],
                              [2.0, 10.0, 29.0,  6.0],
                              [4.0, 11.0, 34.0,  7.0]])

        tree = TreeNode.from_linkage_matrix(linkage, id_list)
        self.assertEquals("(E:17.0,(C:14.5,((A:4.0,D:4.0):4.25,(G:6.25,(B:0.5,"
                          "F:0.5):5.75):2.0):6.25):2.5);",
                          tree.to_newick(with_distances=True))

    def test_from_newick_empty(self):
        obs = TreeNode.from_newick('')
        self.assertTrue(obs.name is None)
        self.assertTrue(obs.length is None)
        self.assertTrue(obs.parent is None)
        self.assertEqual(obs.children, [])
        self.assertTrue(obs.id is None)

    def test_to_newick_single_node(self):
        # single node, no name, with semicolon
        obs = TreeNode().to_newick()
        self.assertEqual(obs, ';')

        # single node, no name, without semicolon
        obs = TreeNode().to_newick(semicolon=False)
        self.assertEqual(obs, '')

        # single node, with name, with semicolon
        obs = TreeNode(name='brofist').to_newick()
        self.assertEqual(obs, 'brofist;')

        # single node, with name, without semicolon
        obs = TreeNode(name='brofist').to_newick(semicolon=False)
        self.assertEqual(obs, 'brofist')


class DndTokenizerTests(TestCase):

    """Tests of the DndTokenizer factory function."""

    def test_gdata(self):
        """DndTokenizer should work as expected on real data"""
        exp = \
            ['(', '(', 'xyz', ':', '0.28124', ',', '(', 'def', ':', '0.24498',
             ',', 'mno', ':', '0.03627', ')', ':', '0.17710', ')', ':',
             '0.04870', ',', 'abc', ':', '0.05925', ',', '(', 'ghi', ':',
             '0.06914', ',', 'jkl', ':', '0.13776', ')', ':', '0.09853', ')',
             ';']
        # split it up for debugging on an item-by-item basis
        obs = list(_dnd_tokenizer(sample))
        self.assertEqual(len(obs), len(exp))
        for i, j in zip(obs, exp):
            self.assertEqual(i, j)
        # try it all in one go
        self.assertEqual(list(_dnd_tokenizer(sample)), exp)

    def test_nonames(self):
        """DndTokenizer should work as expected on trees with no names"""
        exp = ['(', '(', ',', ')', ',', '(', ',', ')', ')', ';']
        obs = list(_dnd_tokenizer(no_names))
        self.assertEqual(obs, exp)

    def test_missing_tip_name(self):
        """DndTokenizer should work as expected on trees with a missing name"""
        exp = ['(', '(', 'a', ',', 'b', ')', ',', '(', 'c', ',', ')', ')', ';']
        obs = list(_dnd_tokenizer(missing_tip_name))
        self.assertEqual(obs, exp)

    def test_minimal(self):
        """DndTokenizer should work as expected a minimal tree without names"""
        exp = ['(', ')', ';']
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
        exp.children[0].append(TreeNode())
        exp.children[0].append(TreeNode())
        exp.children[1].append(TreeNode())
        exp.children[1].append(TreeNode())
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
        exp.children[0].append(TreeNode(name='a'))
        exp.children[0].append(TreeNode(name='b'))
        exp.children[1].append(TreeNode(name='c'))
        exp.children[1].append(TreeNode())
        self.assertEqual(str(obs), str(exp))

    def test_gsingle(self):
        """DndParser should produce a single-child TreeNode on minimal data"""
        t = TreeNode.from_newick(single)
        self.assertEqual(len(t), 1)
        child = t[0]
        self.assertEqual(child.name, 'abc')
        self.assertEqual(child.length, 3)
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
        self.assertEqual(len(t[0]), 0)  # first child is terminal
        self.assertEqual(len(t[1]), 2)  # second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0):6.0);')

    def test_gnodedata(self):
        """DndParser should assign name to internal nodes correctly"""
        t = TreeNode.from_newick(nodedata)
        self.assertEqual(len(t), 2)
        self.assertEqual(len(t[0]), 0)  # first child is terminal
        self.assertEqual(len(t[1]), 2)  # second child has two children
        self.assertEqual(str(t), '(abc:3.0,(def:4.0,ghi:5.0)jkl:6.0);')
        info_dict = {}
        for node in t.traverse():
            info_dict[node.name] = node.length
        self.assertEqual(info_dict['abc'], 3.0)
        self.assertEqual(info_dict['def'], 4.0)
        self.assertEqual(info_dict['ghi'], 5.0)
        self.assertEqual(info_dict['jkl'], 6.0)

    def test_data(self):
        """DndParser should work as expected on real data"""
        t = TreeNode.from_newick(sample)
        self.assertEqual(
            str(t), '((xyz:0.28124,(def:0.24498,mno:0.03627):0.1771):0.0487,'
                    'abc:0.05925,(ghi:0.06914,jkl:0.13776):0.09853);')
        tdata = TreeNode.from_newick(node_data_sample, unescape_name=True)
        self.assertEqual(
            str(tdata), "((xyz:0.28124,(def:0.24498,mno:0.03627)A:0.1771)"
                        "B:0.0487,abc:0.05925,(ghi:0.06914,jkl:0.13776)"
                        "C:0.09853);")

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

        self.assertEqual(tree_unesc.name, 'E')
        self.assertEqual(tree_unesc.children[0].name, 'A a')
        self.assertEqual(tree_unesc.children[1].children[0].name, 'B')
        self.assertEqual(tree_unesc.children[1].children[0].length, 1.0)
        self.assertEqual(tree_unesc.children[1].children[1].name, 'C')
        self.assertEqual(tree_unesc.children[2].name, 'D_e')
        self.assertEqual(tree_unesc.children[2].length, 0.5)

        self.assertEqual(tree_esc.name, 'E')
        self.assertEqual(tree_esc.children[0].name, 'A_a')
        self.assertEqual(tree_esc.children[1].children[0].name, 'B')
        self.assertEqual(tree_esc.children[1].children[0].length, 1.0)
        self.assertEqual(tree_esc.children[1].children[1].name, 'C')
        self.assertEqual(tree_esc.children[2].name, "'D_e'")
        self.assertEqual(tree_esc.children[2].length, 0.5)

        reload_test = tree_esc.to_newick(with_distances=True,
                                         escape_name=False)
        obs = TreeNode.from_newick(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True),
                         tree_esc.to_newick(with_distances=True))
        reload_test = tree_unesc.to_newick(with_distances=True,
                                           escape_name=False)
        obs = TreeNode.from_newick(reload_test, unescape_name=False)
        self.assertEqual(obs.to_newick(with_distances=True),
                         tree_unesc.to_newick(with_distances=True))

    def test_DndParser_list(self):
        """Make sure TreeNode.from_newick can handle list of strings"""
        t_str = ["(A_a,(B:1.0,C)", ",'D_e':0.5)E;"]
        tree_unesc = TreeNode.from_newick(t_str, unescape_name=True)

        self.assertEqual(tree_unesc.name, 'E')
        self.assertEqual(tree_unesc.children[0].name, 'A a')
        self.assertEqual(tree_unesc.children[1].children[0].name, 'B')
        self.assertEqual(tree_unesc.children[1].children[0].length, 1.0)
        self.assertEqual(tree_unesc.children[1].children[1].name, 'C')
        self.assertEqual(tree_unesc.children[2].name, 'D_e')
        self.assertEqual(tree_unesc.children[2].length, 0.5)

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
