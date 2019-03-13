# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main
from collections import defaultdict

import numpy as np
import numpy.testing as npt
from scipy.stats import pearsonr

from skbio import DistanceMatrix, TreeNode
from skbio.tree import (DuplicateNodeError, NoLengthError,
                        TreeError, MissingNodeError, NoParentError)
from skbio.util import RepresentationWarning


class TreeNodeSubclass(TreeNode):
    pass


class TreeTests(TestCase):

    def setUp(self):
        """Prep the self"""
        self.simple_t = TreeNode.read(io.StringIO("((a,b)i1,(c,d)i2)root;"))
        nodes = dict([(x, TreeNode(x)) for x in 'abcdefgh'])
        nodes['a'].append(nodes['b'])
        nodes['b'].append(nodes['c'])
        nodes['c'].append(nodes['d'])
        nodes['c'].append(nodes['e'])
        nodes['c'].append(nodes['f'])
        nodes['f'].append(nodes['g'])
        nodes['a'].append(nodes['h'])
        self.TreeRoot = nodes['a']

        def rev_f(items):
            items.reverse()

        def rotate_f(items):
            tmp = items[-1]
            items[1:] = items[:-1]
            items[0] = tmp

        self.rev_f = rev_f
        self.rotate_f = rotate_f
        self.complex_tree = TreeNode.read(io.StringIO(
            "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);"))

    def test_bug_issue_1416(self):
        tree = TreeNode.read(['(((a,b,f,g),c),d);'])
        new_tree = tree.shear(['a', 'b', 'c', 'f'])

        exp = {'a', 'b', 'c', 'f'}
        obs = {n.name for n in new_tree.tips()}

        self.assertEqual(obs, exp)
        self.assertEqual(id(new_tree), id(new_tree.children[0].parent))
        self.assertEqual(id(new_tree), id(new_tree.children[1].parent))

    def test_observed_node_counts(self):
        """returns observed nodes counts given vector of otu observation counts
        """
        # no OTUs observed
        otu_counts = {}
        expected = defaultdict(int)
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)
        # error on zero count(s)
        otu_counts = {'a': 0}
        self.assertRaises(ValueError, self.simple_t.observed_node_counts,
                          otu_counts)
        otu_counts = {'a': 0, 'b': 0, 'c': 0, 'd': 0}
        self.assertRaises(ValueError, self.simple_t.observed_node_counts,
                          otu_counts)

        # all OTUs observed once
        otu_counts = {'a': 1, 'b': 1, 'c': 1, 'd': 1}
        expected = defaultdict(int)
        expected[self.simple_t.find('root')] = 4
        expected[self.simple_t.find('i1')] = 2
        expected[self.simple_t.find('i2')] = 2
        expected[self.simple_t.find('a')] = 1
        expected[self.simple_t.find('b')] = 1
        expected[self.simple_t.find('c')] = 1
        expected[self.simple_t.find('d')] = 1
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)

        # some OTUs observed twice
        otu_counts = {'a': 2, 'b': 1, 'c': 1, 'd': 1}
        expected = defaultdict(int)
        expected[self.simple_t.find('root')] = 5
        expected[self.simple_t.find('i1')] = 3
        expected[self.simple_t.find('i2')] = 2
        expected[self.simple_t.find('a')] = 2
        expected[self.simple_t.find('b')] = 1
        expected[self.simple_t.find('c')] = 1
        expected[self.simple_t.find('d')] = 1
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)

        otu_counts = {'a': 2, 'b': 1, 'c': 1, 'd': 2}
        expected = defaultdict(int)
        expected[self.simple_t.find('root')] = 6
        expected[self.simple_t.find('i1')] = 3
        expected[self.simple_t.find('i2')] = 3
        expected[self.simple_t.find('a')] = 2
        expected[self.simple_t.find('b')] = 1
        expected[self.simple_t.find('c')] = 1
        expected[self.simple_t.find('d')] = 2
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)

        # some OTUs observed, others not observed
        otu_counts = {'a': 2, 'b': 1}
        expected = defaultdict(int)
        expected[self.simple_t.find('root')] = 3
        expected[self.simple_t.find('i1')] = 3
        expected[self.simple_t.find('a')] = 2
        expected[self.simple_t.find('b')] = 1
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)

        otu_counts = {'d': 1}
        expected = defaultdict(int)
        expected[self.simple_t.find('root')] = 1
        expected[self.simple_t.find('i2')] = 1
        expected[self.simple_t.find('d')] = 1
        self.assertEqual(self.simple_t.observed_node_counts(otu_counts),
                         expected)

        # error on non-tips
        otu_counts = {'a': 2, 'e': 1}
        self.assertRaises(MissingNodeError, self.simple_t.observed_node_counts,
                          otu_counts)
        otu_counts = {'a': 2, 'i1': 1}
        self.assertRaises(MissingNodeError, self.simple_t.observed_node_counts,
                          otu_counts)

        # test with another tree
        otu_counts = {}
        expected = defaultdict(int)
        self.assertEqual(self.complex_tree.observed_node_counts(otu_counts),
                         expected)

        otu_counts = {'e': 42, 'f': 1}
        expected[self.complex_tree.root()] = 43
        expected[self.complex_tree.find('int5')] = 43
        expected[self.complex_tree.find('e')] = 42
        expected[self.complex_tree.find('f')] = 1
        self.assertEqual(self.complex_tree.observed_node_counts(otu_counts),
                         expected)

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
        second_tree = TreeNode.read(io.StringIO("(x,y)z;"))
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
        second_tree = TreeNode.read(io.StringIO("(x1,y1)z1;"))
        third_tree = TreeNode.read(io.StringIO("(x2,y2)z2;"))
        first_tree = TreeNode.read(io.StringIO("(x1,y1)z1;"))
        fourth_tree = TreeNode.read(io.StringIO("(x2,y2)z2;"))
        self.simple_t.extend([second_tree, third_tree])

        first_tree.extend(fourth_tree.children)
        self.assertEqual(0, len(fourth_tree.children))
        self.assertEqual(first_tree.children[0].name, 'x1')
        self.assertEqual(first_tree.children[1].name, 'y1')
        self.assertEqual(first_tree.children[2].name, 'x2')
        self.assertEqual(first_tree.children[3].name, 'y2')

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
        self.assertEqual(str(p), ';\n')
        p.name = 'abc'
        self.assertEqual(str(p), 'abc;\n')
        p.length = 3
        self.assertEqual(str(p), 'abc:3;\n')  # don't suppress branch from root
        q = TreeNode()
        p.append(q)
        self.assertEqual(str(p), '()abc:3;\n')
        r = TreeNode()
        q.append(r)
        self.assertEqual(str(p), '(())abc:3;\n')
        r.name = 'xyz'
        self.assertEqual(str(p), '((xyz))abc:3;\n')
        q.length = 2
        self.assertEqual(str(p), '((xyz):2)abc:3;\n')

    def test_pop(self):
        """Pop off a node"""
        second_tree = TreeNode.read(io.StringIO("(x1,y1)z1;"))
        third_tree = TreeNode.read(io.StringIO("(x2,y2)z2;"))
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
        def f(node):
            return node.name in ['b', 'd']

        self.simple_t.remove_deleted(f)
        exp = "((a)i1,(c)i2)root;\n"
        obs = str(self.simple_t)
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

    def test_shear_prune_parent_dropped(self):
        bugtree = "((a,b),((c,d),(e,f)));"
        to_keep = ['c', 'd']
        exp = "(c,d);\n"
        obs = str(TreeNode.read(io.StringIO(bugtree)).shear(to_keep))
        self.assertEqual(obs, exp)

    def test_prune_nested_single_descendent(self):
        bugtree = "(((a,b)));"
        exp = "(a,b);\n"
        t = TreeNode.read(io.StringIO(bugtree))
        t.prune()
        obs = str(t)
        self.assertEqual(obs, exp)

    def test_prune_root_single_desc(self):
        t = TreeNode.read(["((a,b)c)extra;"])
        exp = "(a,b)c;\n"
        t.prune()
        self.assertEqual(str(t), exp)

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

    def test_invalidate_lookup_caches(self):
        root = self.simple_t
        root.create_caches()
        self.assertNotEqual(root._tip_cache, {})
        self.assertNotEqual(root._non_tip_cache, {})
        root.invalidate_caches()
        self.assertEqual(root._tip_cache, {})
        self.assertEqual(root._non_tip_cache, {})

    def test_invalidate_attr_caches(self):
        tree = TreeNode.read(io.StringIO("((a,b,(c,d)e)f,(g,h)i)root;"))

        def f(n):
            return [n.name] if n.is_tip() else []

        tree.cache_attr(f, 'tip_names')
        tree.invalidate_caches()
        for n in tree.traverse(include_self=True):
            self.assertFalse(hasattr(n, 'tip_names'))

    def test_create_caches_duplicate_tip_names(self):
        with self.assertRaises(DuplicateNodeError):
            TreeNode.read(io.StringIO('(a, a);')).create_caches()

    def test_find_all(self):
        t = TreeNode.read(io.StringIO("((a,b)c,((d,e)c)c,(f,(g,h)c)a)root;"))
        exp = [t.children[0],
               t.children[1].children[0],
               t.children[1],
               t.children[2].children[1]]
        obs = t.find_all('c')
        self.assertEqual(obs, exp)

        identity = t.find_all(t)
        self.assertEqual(len(identity), 1)
        self.assertEqual(identity[0], t)

        identity_name = t.find_all('root')
        self.assertEqual(len(identity_name), 1)
        self.assertEqual(identity_name[0], t)

        exp = [t.children[2],
               t.children[0].children[0]]
        obs = t.find_all('a')
        self.assertEqual(obs, exp)

        with self.assertRaises(MissingNodeError):
            t.find_all('missing')

    def test_find(self):
        """Find a node in a tree"""
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f);"))
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
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f,(g,h)f);"))
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
        t1 = TreeNode.read(io.StringIO("((,),(,,));"))
        t2 = TreeNode.read(io.StringIO("((,),(,,));"))

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
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f);"))

        def func(x):
            return x.parent == t.find('c')

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
        tr = TreeNode.read(io.StringIO("(B:0.2,(C:0.3,D:0.4):0.6)F;"))
        obs = tr.ascii_art(show_internal=True, compact=False)
        exp = ("          /-B\n"
               "-F-------|\n"
               "         |          /-C\n"
               "          \\--------|\n"
               "                    \\-D")
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=True, compact=True)
        exp = ("-F------- /-B\n"
               "          \\-------- /-C\n"
               "                    \\-D")
        self.assertEqual(obs, exp)
        obs = tr.ascii_art(show_internal=False, compact=False)
        exp = ("          /-B\n"
               "---------|\n"
               "         |          /-C\n"
               "          \\--------|\n"
               "                    \\-D")
        self.assertEqual(obs, exp)

    def test_ascii_art_with_support(self):
        """Make some ascii trees with support values"""
        tr = TreeNode.read(io.StringIO("(B:0.2,(C:0.3,D:0.4)90:0.6)F;"))
        exp = "          /-B\n-F-------|\n         |          /-C\n         "\
              " \\90------|\n                    \\-D"
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr.assign_supports()
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr = TreeNode.read(io.StringIO("((A,B)75,(C,D)'80:spA');"))
        exp = "                    /-A\n          /75------|\n         |    "\
              "      \\-B\n---------|\n         |          /-C\n          \\"\
              "80:spA--|\n                    \\-D"
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr.assign_supports()
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)

    def test_ascii_art_three_children(self):
        obs = TreeNode.read(io.StringIO('(a,(b,c,d));')).ascii_art()
        self.assertEqual(obs, exp_ascii_art_three_children)

    def test_accumulate_to_ancestor(self):
        """Get the distance from a node to its ancestor"""
        t = TreeNode.read(io.StringIO(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;"))
        a = t.find('a')
        b = t.find('b')
        exp_to_root = 0.1 + 0.3
        obs_to_root = a.accumulate_to_ancestor(t)
        self.assertEqual(obs_to_root, exp_to_root)

        with self.assertRaises(NoParentError):
            a.accumulate_to_ancestor(b)

    def test_distance_nontip(self):
        # example derived from issue #807, credit @wwood
        tstr = "((A:1.0,B:2.0)'g__genus1':3.0)root;"
        tree = TreeNode.read(io.StringIO(tstr))
        self.assertEqual(tree.find('A').distance(tree.find('g__genus1')), 1.0)

    def test_distance(self):
        """Get the distance between two nodes"""
        t = TreeNode.read(io.StringIO(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;"))
        tips = sorted([n for n in t.tips()], key=lambda x: x.name)

        npt.assert_almost_equal(tips[0].distance(tips[0]), 0.0)
        npt.assert_almost_equal(tips[0].distance(tips[1]), 0.3)
        npt.assert_almost_equal(tips[0].distance(tips[2]), 1.3)
        with self.assertRaises(NoLengthError):
            tips[0].distance(tips[3])

        npt.assert_almost_equal(tips[1].distance(tips[0]), 0.3)
        npt.assert_almost_equal(tips[1].distance(tips[1]), 0.0)
        npt.assert_almost_equal(tips[1].distance(tips[2]), 1.4)
        with self.assertRaises(NoLengthError):
            tips[1].distance(tips[3])

        self.assertEqual(tips[2].distance(tips[0]), 1.3)
        self.assertEqual(tips[2].distance(tips[1]), 1.4)
        self.assertEqual(tips[2].distance(tips[2]), 0.0)
        with self.assertRaises(NoLengthError):
            tips[2].distance(tips[3])

    def test_lowest_common_ancestor(self):
        """TreeNode lowestCommonAncestor should return LCA for set of tips"""
        t1 = TreeNode.read(io.StringIO("((a,(b,c)d)e,f,(g,h)i)j;"))
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
        tree = TreeNode.read(io.StringIO(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;"))
        dist, nodes = tree.get_max_distance()
        npt.assert_almost_equal(dist, 1.6)
        self.assertEqual(sorted([n.name for n in nodes]), ['b', 'e'])

    def test_set_max_distance(self):
        """set_max_distance sets MaxDistTips across tree"""
        tree = TreeNode.read(io.StringIO(
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;"))
        tree._set_max_distance()
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 1.6)
        self.assertEqual(sorted([tip_a[1].name, tip_b[1].name]), ['b', 'e'])

    def test_set_max_distance_tie_bug(self):
        """Corresponds to #1077"""
        s = io.StringIO("((a:1,b:1)c:2,(d:3,e:4)f:5)root;")
        t = TreeNode.read(s)

        exp = ((3.0, t.find('a')), (9.0, t.find('e')))

        # the above tree would trigger an exception in max. The central issue
        # was that the data being passed to max were a tuple of tuple:
        # ((left_d, left_n), (right_d, right_n))
        # the call to max would break in this scenario as it would fall onto
        # idx 1 of each tuple to assess the "max".
        t._set_max_distance()

        self.assertEqual(t.MaxDistTips, exp)

    def test_set_max_distance_inplace_modification_bug(self):
        """Corresponds to #1223"""
        s = io.StringIO("((a:1,b:1)c:2,(d:3,e:4)f:5)root;")
        t = TreeNode.read(s)

        exp = [((0.0, t.find('a')), (0.0, t.find('a'))),
               ((0.0, t.find('b')), (0.0, t.find('b'))),
               ((1.0, t.find('a')), (1.0, t.find('b'))),
               ((0.0, t.find('d')), (0.0, t.find('d'))),
               ((0.0, t.find('e')), (0.0, t.find('e'))),
               ((3.0, t.find('d')), (4.0, t.find('e'))),
               ((3.0, t.find('a')), (9.0, t.find('e')))]

        t._set_max_distance()

        self.assertEqual([n.MaxDistTips for n in t.postorder()], exp)

    def test_shear(self):
        """Shear the nodes"""
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        obs = str(t.shear(['G', 'M']))
        exp = '(G:3.0,M:3.7);\n'
        self.assertEqual(obs, exp)

    def test_compare_tip_distances(self):
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        t2 = TreeNode.read(io.StringIO('(((H:1,G:1,O:1):2,R:3):1,X:4);'))
        obs = t.compare_tip_distances(t2)
        # note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

    def test_compare_tip_distances_sample(self):
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        t2 = TreeNode.read(io.StringIO('(((H:1,G:1,O:1):2,R:3):1,X:4);'))
        obs = t.compare_tip_distances(t2, sample=3, shuffle_f=sorted)
        # note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

        # 4 common taxa, still picking H, G, R
        s = '((H:1,G:1):2,(R:0.5,M:0.7,Q:5):3);'
        t = TreeNode.read(io.StringIO(s))
        s3 = '(((H:1,G:1,O:1):2,R:3,Q:10):1,X:4);'
        t3 = TreeNode.read(io.StringIO(s3))
        obs = t.compare_tip_distances(t3, sample=3, shuffle_f=sorted)

    def test_compare_tip_distances_no_common_tips(self):
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        t2 = TreeNode.read(io.StringIO('(((Z:1,Y:1,X:1):2,W:3):1,V:4);'))

        with self.assertRaises(ValueError):
            t.compare_tip_distances(t2)

    def test_compare_tip_distances_single_common_tip(self):
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        t2 = TreeNode.read(io.StringIO('(((R:1,Y:1,X:1):2,W:3):1,V:4);'))

        self.assertEqual(t.compare_tip_distances(t2), 1)
        self.assertEqual(t2.compare_tip_distances(t), 1)

    def test_tip_tip_distances_endpoints(self):
        """Test getting specifc tip distances  with tipToTipDistances"""
        t = TreeNode.read(io.StringIO('((H:1,G:1):2,(R:0.5,M:0.7):3);'))
        nodes = [t.find('H'), t.find('G'), t.find('M')]
        names = ['H', 'G', 'M']
        exp = DistanceMatrix(np.array([[0, 2.0, 6.7],
                                       [2.0, 0, 6.7],
                                       [6.7, 6.7, 0.0]]), ['H', 'G', 'M'])

        obs = t.tip_tip_distances(endpoints=names)
        self.assertEqual(obs, exp)

        obs = t.tip_tip_distances(endpoints=nodes)
        self.assertEqual(obs, exp)

    def test_tip_tip_distances_non_tip_endpoints(self):
        t = TreeNode.read(io.StringIO('((H:1,G:1)foo:2,(R:0.5,M:0.7):3);'))
        with self.assertRaises(ValueError):
            t.tip_tip_distances(endpoints=['foo'])

    def test_tip_tip_distances_no_length(self):
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f);"))
        exp_t = TreeNode.read(io.StringIO("((a:0,b:0)c:0,(d:0,e:0)f:0);"))
        exp_t_dm = exp_t.tip_tip_distances()

        t_dm = npt.assert_warns(RepresentationWarning, t.tip_tip_distances)
        self.assertEqual(t_dm, exp_t_dm)

        for node in t.preorder():
            self.assertIs(node.length, None)

    def test_tip_tip_distances_missing_length(self):
        t = TreeNode.read(io.StringIO("((a,b:6)c:4,(d,e:0)f);"))
        exp_t = TreeNode.read(io.StringIO("((a:0,b:6)c:4,(d:0,e:0)f:0);"))
        exp_t_dm = exp_t.tip_tip_distances()

        t_dm = npt.assert_warns(RepresentationWarning, t.tip_tip_distances)
        self.assertEqual(t_dm, exp_t_dm)

    def test_neighbors(self):
        """Get neighbors of a node"""
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f);"))
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
        t = TreeNode.read(io.StringIO("((a,b)c,(d,e)f);"))
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

    def test_tips_self(self):
        """ See issue #1509 """
        tree = TreeNode.read(['(c, (b,a)x)y;'])
        ts = list(tree.find('c').tips(include_self=True))
        self.assertEqual(len(ts), 1)
        t = ts[0]
        self.assertEqual(t.name, 'c')
        self.assertTrue(t.is_tip())

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

    def test_bifurcate(self):
        t1 = TreeNode.read(io.StringIO('(((a,b),c),(d,e));'))
        t2 = TreeNode.read(io.StringIO('((a,b,c));'))
        t3 = t2.copy()

        t1.bifurcate()
        t2.bifurcate()
        t3.bifurcate(insert_length=0)

        self.assertEqual(str(t1), '(((a,b),c),(d,e));\n')
        self.assertEqual(str(t2), '((c,(a,b)));\n')
        self.assertEqual(str(t3), '((c,(a,b):0));\n')

    def test_bifurcate_with_subclass(self):
        tree = TreeNodeSubclass()
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())

        tree.bifurcate()

        for node in tree.traverse():
            self.assertIs(type(node), TreeNodeSubclass)

    def test_index_tree_single_node(self):
        """index_tree handles single node tree"""
        t1 = TreeNode.read(io.StringIO('root;'))
        id_index, child_index = t1.index_tree()
        self.assertEqual(id_index[0], t1)
        npt.assert_equal(child_index, np.array([[]]))

    def test_index_tree(self):
        """index_tree should produce correct index and node map"""
        # test for first tree: contains singleton outgroup
        t1 = TreeNode.read(io.StringIO('(((a,b),c),(d,e));'))
        t2 = TreeNode.read(io.StringIO('(((a,b),(c,d)),(e,f));'))
        t3 = TreeNode.read(io.StringIO('(((a,b,c),(d)),(e,f));'))

        id_1, child_1 = t1.index_tree()
        nodes_1 = [n.id for n in t1.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_1, [0, 1, 2, 3, 6, 4, 5, 7, 8])
        npt.assert_equal(child_1, np.array([[2, 0, 1], [6, 2, 3], [7, 4, 5],
                                            [8, 6, 7]]))

        # test for second tree: strictly bifurcating
        id_2, child_2 = t2.index_tree()
        nodes_2 = [n.id for n in t2.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_2, [0, 1, 4, 2, 3, 5, 8, 6, 7, 9, 10])
        npt.assert_equal(child_2, np.array([[4, 0, 1], [5, 2, 3],
                                            [8, 4, 5], [9, 6, 7],
                                            [10, 8, 9]]))

        # test for third tree: contains trifurcation and single-child parent
        id_3, child_3 = t3.index_tree()
        nodes_3 = [n.id for n in t3.traverse(self_before=False,
                   self_after=True)]
        self.assertEqual(nodes_3, [0, 1, 2, 4, 3, 5, 8, 6, 7, 9, 10])
        npt.assert_equal(child_3, np.array([[4, 0, 2], [5, 3, 3], [8, 4, 5],
                                            [9, 6, 7], [10, 8, 9]]))

    def test_root_at(self):
        """Form a new root"""
        t = TreeNode.read(io.StringIO("(((a,b)c,(d,e)f)g,h)i;"))
        with self.assertRaises(TreeError):
            t.root_at(t.find('h'))

        exp = "(a,b,((d,e)f,(h)g)c)root;\n"
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

    def test_root_at_midpoint_no_lengths(self):
        # should get same tree back (a copy)
        nwk = '(a,b)c;\n'
        t = TreeNode.read(io.StringIO(nwk))
        obs = t.root_at_midpoint()
        self.assertEqual(str(obs), nwk)

    def test_root_at_midpoint_tie(self):
        nwk = "(((a:1,b:1)c:2,(d:3,e:4)f:5),g:1)root;"
        t = TreeNode.read(io.StringIO(nwk))
        exp = "((d:3,e:4)f:2,((a:1,b:1)c:2,(g:1)):3)root;"
        texp = TreeNode.read(io.StringIO(exp))

        obs = t.root_at_midpoint()

        for o, e in zip(obs.traverse(), texp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

    def test_compare_subsets(self):
        """compare_subsets should return the fraction of shared subsets"""
        t = TreeNode.read(io.StringIO('((H,G),(R,M));'))
        t2 = TreeNode.read(io.StringIO('(((H,G),R),M);'))
        t4 = TreeNode.read(io.StringIO('(((H,G),(O,R)),X);'))

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
        t = TreeNode.read(io.StringIO('((H,G),(R,M));'))
        t2 = TreeNode.read(io.StringIO('(((H,G),R),M);'))
        t4 = TreeNode.read(io.StringIO('(((H,G),(O,R)),X);'))

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
        t1 = TreeNode.read(io.StringIO("(((a,b),c),(e,f),(g));"))
        t2 = TreeNode.read(io.StringIO("(((a,b),c),(e,f),(g));"))
        t3 = TreeNode.read(io.StringIO("((g),(e,f),(c,(a,b)));"))
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
        t1 = TreeNode.read(io.StringIO('(((a,b),c),(d,e));'))
        t2 = TreeNode.read(io.StringIO('(((a,b),(c,d)),(e,f));'))
        t3 = TreeNode.read(io.StringIO('(((a,b,c),(d)),(e,f));'))
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
        t = TreeNode.read(io.StringIO("((a,(b,c)d)e,(f,g)h)i;"))
        exp = "(b,c,(a,((f,g)h)e)d)root;\n"
        obs = t.find('d').unrooted_deepcopy()
        self.assertEqual(str(obs), exp)

        t_ids = {id(n) for n in t.traverse()}
        obs_ids = {id(n) for n in obs.traverse()}

        self.assertEqual(t_ids.intersection(obs_ids), set())

    def test_descending_branch_length(self):
        """Calculate descending branch_length"""
        tr = TreeNode.read(io.StringIO(
            "(((A:.1,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"))
        tdbl = tr.descending_branch_length()
        sdbl = tr.descending_branch_length(['A', 'E'])
        npt.assert_almost_equal(tdbl, 8.9)
        npt.assert_almost_equal(sdbl, 2.2)
        self.assertRaises(ValueError, tr.descending_branch_length,
                          ['A', 'DNE'])
        self.assertRaises(ValueError, tr.descending_branch_length, ['A', 'C'])

        tr = TreeNode.read(io.StringIO(
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"))
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 8.8)

        tr = TreeNode.read(io.StringIO(
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:.5)J:1.3)K;"))
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 7.9)

        tr = TreeNode.read(io.StringIO(
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:.5)J:1.3)K;"))
        tdbl = tr.descending_branch_length(['A', 'D', 'E'])
        npt.assert_almost_equal(tdbl, 2.1)

        tr = TreeNode.read(io.StringIO(
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"))
        tdbl = tr.descending_branch_length(['I', 'D', 'E'])
        npt.assert_almost_equal(tdbl, 6.6)

        # test with a situation where we have unnamed internal nodes
        tr = TreeNode.read(io.StringIO(
            "(((A,B:1.2):.6,(D:.9,E:.6)F):2.4,(H:.4,I:.5)J:1.3);"))
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 7.9)

    def test_to_array(self):
        """Convert a tree to arrays"""
        t = TreeNode.read(io.StringIO(
            '(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10);'))
        id_index, child_index = t.index_tree()
        arrayed = t.to_array()

        self.assertEqual(id_index, arrayed['id_index'])
        npt.assert_equal(child_index, arrayed['child_index'])

        exp = np.array([1, 2, 3, 5, 4, 6, 8, 9, 7, 10, np.nan])
        obs = arrayed['length']
        npt.assert_equal(obs, exp)

        exp = np.array(['a', 'b', 'c', 'd', 'x',
                        'y', 'e', 'f', 'z', 'z', None])
        obs = arrayed['name']
        npt.assert_equal(obs, exp)

        exp = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        obs = arrayed['id']
        npt.assert_equal(obs, exp)

    def test_to_array_attrs(self):
        t = TreeNode.read(io.StringIO(
            '(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10);'))
        id_index, child_index = t.index_tree()
        arrayed = t.to_array(attrs=[('name', object)])

        # should only have id_index, child_index, and name since we specified
        # attrs
        self.assertEqual(len(arrayed), 3)

        self.assertEqual(id_index, arrayed['id_index'])
        npt.assert_equal(child_index, arrayed['child_index'])

        exp = np.array(['a', 'b', 'c', 'd', 'x',
                        'y', 'e', 'f', 'z', 'z', None])
        obs = arrayed['name']
        npt.assert_equal(obs, exp)

        # invalid attrs
        with self.assertRaises(AttributeError):
            t.to_array(attrs=[('name', object), ('brofist', int)])

    def test_to_array_nan_length_value(self):
        t = TreeNode.read(io.StringIO("((a:1, b:2)c:3)root;"))
        indexed = t.to_array(nan_length_value=None)
        npt.assert_equal(indexed['length'],
                         np.array([1, 2, 3, np.nan], dtype=float))
        indexed = t.to_array(nan_length_value=0.0)
        npt.assert_equal(indexed['length'],
                         np.array([1, 2, 3, 0.0], dtype=float))
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed['length'],
                         np.array([1, 2, 3, 42.0], dtype=float))

        t = TreeNode.read(io.StringIO("((a:1, b:2)c:3)root:4;"))
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed['length'],
                         np.array([1, 2, 3, 4], dtype=float))

        t = TreeNode.read(io.StringIO("((a:1, b:2)c)root;"))
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed['length'],
                         np.array([1, 2, 42.0, 42.0], dtype=float))

    def test_from_taxonomy(self):
        input_lineages = {'1': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                          '2': ['a', 'b', 'c', None, None, 'x', 'y'],
                          '3': ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                          '4': ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                          '5': ['h', 'i', 'j', 'k', 'l', 'm', 'n']}
        exp = TreeNode.read(io.StringIO(
            "((((((((1)g)f)e)d,((((2)y)x)))c)b)a,"
            "(((((((3,5)n,(4)q)m)l)k)j)i)h);"))

        root = TreeNode.from_taxonomy(input_lineages.items())

        self.assertIs(type(root), TreeNode)

        self.assertEqual(root.compare_subsets(exp), 0.0)

        root = TreeNodeSubclass.from_taxonomy(input_lineages.items())

        self.assertIs(type(root), TreeNodeSubclass)

    def test_to_taxonomy(self):
        input_lineages = {'1': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                          '2': ['a', 'b', 'c', None, None, 'x', 'y'],
                          '3': ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                          '4': ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                          '5': ['h', 'i', 'j', 'k', 'l', 'm', 'n']}
        tree = TreeNode.from_taxonomy(input_lineages.items())
        exp = sorted(input_lineages.items())
        obs = [(n.name, lin) for n, lin in tree.to_taxonomy(allow_empty=True)]
        self.assertEqual(sorted(obs), exp)

    def test_to_taxonomy_filter(self):
        input_lineages = {'1': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                          '2': ['a', 'b', 'c', None, None, 'x', 'y'],
                          '3': ['h', 'i', 'j', 'k', 'l'],  # test jagged
                          '4': ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                          '5': ['h', 'i', 'j', 'k', 'l', 'm', 'n']}
        tree = TreeNode.from_taxonomy(input_lineages.items())

        def f(node, lin):
            return 'k' in lin or 'x' in lin

        exp = [('2', ['a', 'b', 'c', 'x', 'y']),
               ('3', ['h', 'i', 'j', 'k', 'l']),
               ('4', ['h', 'i', 'j', 'k', 'l', 'm', 'q']),
               ('5', ['h', 'i', 'j', 'k', 'l', 'm', 'n'])]
        obs = [(n.name, lin) for n, lin in tree.to_taxonomy(filter_f=f)]
        self.assertEqual(sorted(obs), exp)

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

        self.assertIs(type(tree), TreeNode)

        self.assertEqual("(E:17.0,(C:14.5,((A:4.0,D:4.0):4.25,(G:6.25,(B:0.5,"
                         "F:0.5):5.75):2.0):6.25):2.5);\n",
                         str(tree))

        tree = TreeNodeSubclass.from_linkage_matrix(linkage, id_list)

        self.assertIs(type(tree), TreeNodeSubclass)

    def test_shuffle_invalid_iter(self):
        shuffler = self.simple_t.shuffle(n=-1)
        with self.assertRaises(ValueError):
            next(shuffler)

    def test_shuffle_n_2(self):
        exp = ["((a,b)i1,(d,c)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((a,b)i1,(d,c)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((a,b)i1,(d,c)i2)root;\n"]

        obs_g = self.simple_t.shuffle(k=2, shuffle_f=self.rev_f, n=np.inf)
        obs = [str(next(obs_g)) for i in range(5)]
        self.assertEqual(obs, exp)

    def test_shuffle_n_none(self):
        exp = ["((d,c)i1,(b,a)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((d,c)i1,(b,a)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n"]
        obs_g = self.simple_t.shuffle(shuffle_f=self.rev_f, n=4)
        obs = [str(next(obs_g)) for i in range(4)]
        self.assertEqual(obs, exp)

    def test_shuffle_complex(self):
        exp = ["(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);\n"]

        obs_g = self.complex_tree.shuffle(shuffle_f=self.rev_f,
                                          names=['c', 'd', 'e', 'f'], n=4)
        obs = [str(next(obs_g)) for i in range(4)]
        self.assertEqual(obs, exp)

    def test_shuffle_names(self):
        exp = ["((c,a)i1,(b,d)i2)root;\n",
               "((b,c)i1,(a,d)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((c,a)i1,(b,d)i2)root;\n"]

        obs_g = self.simple_t.shuffle(names=['a', 'b', 'c'],
                                      shuffle_f=self.rotate_f, n=np.inf)
        obs = [str(next(obs_g)) for i in range(4)]
        self.assertEqual(obs, exp)

    def test_shuffle_raises(self):
        with self.assertRaises(ValueError):
            next(self.simple_t.shuffle(k=1))

        with self.assertRaises(ValueError):
            next(self.simple_t.shuffle(k=5, names=['a', 'b']))

        with self.assertRaises(MissingNodeError):
            next(self.simple_t.shuffle(names=['x', 'y']))

    def test_assign_supports(self):
        """Extract support values of internal nodes."""
        # test nodes with support values alone as labels
        tree = TreeNode.read(['((a,b)75,(c,d)90);'])
        tree.assign_supports()
        node1, node2 = tree.children
        # check if internal nodes are assigned correct support values
        self.assertEqual(node1.support, 75)
        self.assertEqual(node2.support, 90)
        # check if original node names are cleared
        self.assertIsNone(node1.name)
        self.assertIsNone(node2.name)
        # check if support values are not assigned to root and tips
        self.assertIsNone(tree.support)
        for taxon in ('a', 'b', 'c', 'd'):
            self.assertIsNone(tree.find(taxon).support)

        # test nodes with support values and branch lengths
        tree = TreeNode.read(['((a,b)0.85:1.23,(c,d)0.95:4.56);'])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, 0.85)
        self.assertEqual(node2.support, 0.95)

        # test whether integer or float support values can be correctly parsed
        tree = TreeNode.read(['((a,b)75,(c,d)80.0,(e,f)97.5,(g,h)0.95);'])
        tree.assign_supports()
        node1, node2, node3, node4 = tree.children
        self.assertTrue(isinstance(node1.support, int))
        self.assertEqual(node1.support, 75)
        self.assertTrue(isinstance(node2.support, float))
        self.assertEqual(node2.support, 80.0)
        self.assertTrue(isinstance(node3.support, float))
        self.assertEqual(node3.support, 97.5)
        self.assertTrue(isinstance(node4.support, float))
        self.assertEqual(node4.support, 0.95)

        # test support values that are negative or scientific notation (not a
        # common scenario but can happen)
        tree = TreeNode.read(['((a,b)-1.23,(c,d)1.23e-4);'])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, -1.23)
        self.assertEqual(node2.support, 0.000123)

        # test nodes with support and extra label
        tree = TreeNode.read(['((a,b)\'80:X\',(c,d)\'60:Y\');'])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, 80)
        self.assertEqual(node1.name, 'X')
        self.assertEqual(node2.support, 60)
        self.assertEqual(node2.name, 'Y')

        # test nodes without label, with non-numeric label, and with branch
        # length only
        tree = TreeNode.read(['((a,b),(c,d)x,(e,f):1.0);'])
        tree.assign_supports()
        for node in tree.children:
            self.assertIsNone(node.support)

    def test_unpack(self):
        """Unpack an internal node."""
        # test unpacking a node without branch length
        tree = TreeNode.read(['((c,d)a,(e,f)b);'])
        tree.find('b').unpack()
        exp = '((c,d)a,e,f);\n'
        self.assertEqual(str(tree), exp)

        # test unpacking a node with branch length
        tree = TreeNode.read(['((c:2.0,d:3.0)a:1.0,(e:2.0,f:1.0)b:2.0);'])
        tree.find('b').unpack()
        exp = '((c:2.0,d:3.0)a:1.0,e:4.0,f:3.0);'
        self.assertEqual(str(tree).rstrip(), exp)

        # test attempting to unpack root
        tree = TreeNode.read(['((d,e)b,(f,g)c)a;'])
        msg = 'Cannot unpack root.'
        with self.assertRaisesRegex(TreeError, msg):
            tree.find('a').unpack()

        # test attempting to unpack tip
        msg = 'Cannot unpack tip.'
        with self.assertRaisesRegex(TreeError, msg):
            tree.find('d').unpack()

    def test_unpack_by_func(self):
        """Unpack internal nodes of a tree by a function."""
        # unpack internal nodes with branch length <= 1.0
        def func(x):
            return x.length <= 1.0

        # will unpack node 'a', but not tip 'e'
        # will add the branch length of 'a' to its child nodes 'c' and 'd'
        tree = TreeNode.read(['((c:2,d:3)a:1,(e:1,f:2)b:2);'])
        tree.unpack_by_func(func)
        exp = '((e:1.0,f:2.0)b:2.0,c:3.0,d:4.0);'
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack internal nodes with branch length < 2.01
        # will unpack both 'a' and 'b'
        tree = TreeNode.read(['((c:2,d:3)a:1,(e:1,f:2)b:2);'])
        tree.unpack_by_func(lambda x: x.length <= 2.0)
        exp = '(c:3.0,d:4.0,e:3.0,f:4.0);'
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack two nested nodes 'a' and 'c' simultaneously
        tree = TreeNode.read(['(((e:3,f:2)c:1,d:3)a:1,b:4);'])
        tree.unpack_by_func(lambda x: x.length <= 2.0)
        exp = '(b:4.0,d:4.0,e:5.0,f:4.0);'
        self.assertEqual(str(tree).rstrip(), exp)

        # test a complicated scenario (unpacking nodes 'g', 'h' and 'm')
        def func(x):
            return x.length < 2.0
        tree = TreeNode.read(['(((a:1.04,b:2.32,c:1.44)d:3.20,'
                              '(e:3.91,f:2.47)g:1.21)h:1.75,'
                              '(i:4.14,(j:2.06,k:1.58)l:3.32)m:0.77);'])
        tree.unpack_by_func(func)
        exp = ('((a:1.04,b:2.32,c:1.44)d:4.95,e:6.87,f:5.43,i:4.91,'
               '(j:2.06,k:1.58)l:4.09);')
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 75
        def func(x):
            return x.support < 75
        tree = TreeNode.read(['(((a,b)85,(c,d)78)75,(e,(f,g)64)80);'])
        tree.assign_supports()
        tree.unpack_by_func(func)
        exp = '(((a,b)85,(c,d)78)75,(e,f,g)80);'
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 85
        tree = TreeNode.read(['(((a,b)85,(c,d)78)75,(e,(f,g)64)80);'])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support < 85)
        exp = '((a,b)85,c,d,e,f,g);'
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 0.95
        tree = TreeNode.read(['(((a,b)0.97,(c,d)0.98)1.0,(e,(f,g)0.88)0.96);'])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support < 0.95)
        exp = '(((a,b)0.97,(c,d)0.98)1.0,(e,f,g)0.96);'
        self.assertEqual(str(tree).rstrip(), exp)

        # test a case where there are branch lengths, none support values and
        # node labels
        tree = TreeNode.read(['(((a:1.02,b:0.33)85:0.12,(c:0.86,d:2.23)'
                              '70:3.02)75:0.95,(e:1.43,(f:1.69,g:1.92)64:0.20)'
                              'node:0.35)root;'])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support is not None and x.support < 75)
        exp = ('(((a:1.02,b:0.33)85:0.12,c:3.88,d:5.25)75:0.95,'
               '(e:1.43,f:1.89,g:2.12)node:0.35)root;')
        self.assertEqual(str(tree).rstrip(), exp)


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

exp_ascii_art_three_children = r"""          /-a
         |
---------|          /-b
         |         |
          \--------|--c
                   |
                    \-d"""


if __name__ == '__main__':
    main()
