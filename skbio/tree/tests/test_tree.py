# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from collections import defaultdict

import numpy as np
import numpy.testing as npt
import pandas as pd
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
        # a simple tree
        self.simple_t = TreeNode.read(["((a,b)i1,(c,d)i2)root;"])
        #                     /-a
        #           /i1------|
        #          |          \-b
        # -root----|
        #          |          /-c
        #           \i2------|
        #                     \-d

        # another test tree
        nodes = dict([(x, TreeNode(x)) for x in "abcdefgh"])
        nodes["a"].append(nodes["b"])
        nodes["b"].append(nodes["c"])
        nodes["c"].append(nodes["d"])
        nodes["c"].append(nodes["e"])
        nodes["c"].append(nodes["f"])
        nodes["f"].append(nodes["g"])
        nodes["a"].append(nodes["h"])
        self.TreeRoot = nodes["a"]
        # (((d,e,(g)f)c)b,h)a;
        #                               /-d
        #                              |
        #           /b------- /c-------|--e
        #          |                   |
        # -a-------|                    \f------- /-g
        #          |
        #           \-h

        self.complex_tree = TreeNode.read([
            "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);"])

    def test_gops(self):
        """Basic TreeNode operations should work as expected."""
        p = TreeNode()
        self.assertEqual(str(p), ";\n")
        p.name = "abc"
        self.assertEqual(str(p), "abc;\n")
        p.length = 3
        self.assertEqual(str(p), "abc:3;\n")  # don"t suppress branch from root
        q = TreeNode()
        p.append(q)
        self.assertEqual(str(p), "()abc:3;\n")
        r = TreeNode()
        q.append(r)
        self.assertEqual(str(p), "(())abc:3;\n")
        r.name = "xyz"
        self.assertEqual(str(p), "((xyz))abc:3;\n")
        q.length = 2
        self.assertEqual(str(p), "((xyz):2)abc:3;\n")

    def test_iter(self):
        """iter wraps children."""
        exp = ["i1", "i2"]
        obs = [n.name for n in self.simple_t]
        self.assertEqual(obs, exp)

    # ------------------------------------------------
    # Tree copying
    # ------------------------------------------------

    def test_copy(self):
        """Copy a tree."""
        t = self.simple_t
        t.children[0].length = 1.2
        t.children[1].children[0].length = 0.5
        cp = t.copy()
        for obs, exp in zip(cp.traverse(), t.traverse()):
            self.assertIsNot(obs, exp)
            self.assertEqual(obs.name, exp.name)
            self.assertEqual(obs.length, exp.length)

        # deep vs shallow copy
        t.dummy = [1, [2, 3], 4]
        cp = t.copy()
        cp.dummy[1].append(0)
        self.assertListEqual(t.dummy[1], [2, 3])
        cp = t.copy(deep=False)
        cp.dummy[1].append(0)
        self.assertListEqual(t.dummy[1], [2, 3, 0])

        # with attribute cache
        t.cache_attr(lambda n: 1, "node_count", sum)
        cp = t.copy()
        self.assertFalse(hasattr(cp, "node_count"))

    def test_deepcopy(self):
        t = self.simple_t
        t.dummy = [1, [2, 3], 4]
        cp = t.deepcopy()
        cp.dummy[1].append(0)
        self.assertListEqual(t.dummy[1], [2, 3])

    def test__copy__(self):
        t = self.simple_t
        t.dummy = [1, [2, 3], 4]
        cp = self.simple_t.__copy__()
        for obs, exp in zip(cp.traverse(), t.traverse()):
            self.assertIsNot(obs, exp)
            self.assertEqual(obs.name, exp.name)
            self.assertEqual(obs.length, exp.length)
        cp.dummy[1].append(0)
        self.assertListEqual(t.dummy[1], [2, 3, 0])

    def test__deepcopy__(self):
        t = self.simple_t
        t.dummy = [1, [2, 3], 4]
        cp = self.simple_t.__deepcopy__({})
        for obs, exp in zip(cp.traverse(), t.traverse()):
            self.assertIsNot(obs, exp)
            self.assertEqual(obs.name, exp.name)
            self.assertEqual(obs.length, exp.length)
        cp.dummy[1].append(0)
        self.assertListEqual(t.dummy[1], [2, 3])

    def test_subtree(self):
        """Make a copy of a subtree."""
        with self.assertRaises(NotImplementedError):
            self.simple_t.children[0].subtree()

    # ------------------------------------------------
    # Tree navigation
    # ------------------------------------------------

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

    def test_has_children(self):
        """Test if has children"""
        t = TreeNode.read(["((a,b)c,(d,e)f);"])
        self.assertTrue(t.has_children())
        self.assertTrue(t.children[0].has_children())
        self.assertTrue(t.children[1].has_children())
        self.assertFalse(t.children[0].children[0].has_children())
        self.assertFalse(t.children[0].children[1].has_children())
        self.assertFalse(t.children[1].children[0].has_children())
        self.assertFalse(t.children[1].children[1].has_children())

    def test_root(self):
        """Get the root!"""
        t = self.simple_t
        self.assertIs(t, self.simple_t.root())
        self.assertIs(t, self.simple_t.children[0].root())
        self.assertIs(t, self.simple_t.children[1].children[1].root())

    def test_ancestors(self):
        """Get all the ancestors"""
        exp = ["i1", "root"]
        obs = self.simple_t.children[0].children[0].ancestors()
        self.assertEqual([o.name for o in obs], exp)

        exp = ["root"]
        obs = self.simple_t.children[0].ancestors()
        self.assertEqual([o.name for o in obs], exp)

        exp = []
        obs = self.simple_t.ancestors()
        self.assertEqual([o.name for o in obs], exp)

    def test_siblings(self):
        """Get siblings of a node."""
        exp = []
        obs = self.simple_t.siblings()
        self.assertEqual(obs, exp)

        exp = ["i2"]
        obs = self.simple_t.children[0].siblings()
        self.assertEqual([o.name for o in obs], exp)

        exp = ["c"]
        obs = self.simple_t.children[1].children[1].siblings()
        self.assertEqual([o.name for o in obs], exp)

        self.simple_t.append(TreeNode(name="foo"))
        self.simple_t.append(TreeNode(name="bar"))
        exp = ["i1", "foo", "bar"]
        obs = self.simple_t.children[1].siblings()
        self.assertEqual([o.name for o in obs], exp)

    def test_neighbors(self):
        """Get neighbors of a node"""
        t = TreeNode.read(["((a,b)c,(d,e)f);"])
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

    def test_lowest_common_ancestor(self):
        """TreeNode lowestCommonAncestor should return LCA for set of tips"""
        t1 = TreeNode.read(["((a,(b,c)d)e,f,(g,h)i)j;"])
        t2 = t1.copy()
        t3 = t1.copy()
        t4 = t1.copy()
        input1 = ["a"]  # return self
        input2 = ["a", "b"]  # return e
        input3 = ["b", "c"]  # return d
        input4 = ["a", "h", "g"]  # return j
        exp1 = t1.find("a")
        exp2 = t2.find("e")
        exp3 = t3.find("d")
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
        exp_1 = t_mul.find("d")
        exp_2 = t_mul.find("i")
        obs_1 = t_mul.lowest_common_ancestor(["b", "c"])
        obs_2 = t_mul.lowest_common_ancestor(["g", "h"])
        self.assertEqual(obs_1, exp_1)
        self.assertEqual(obs_2, exp_2)

        # root included
        t_root = TreeNode.read(["(a,b)c;"])
        obs = t_root.lowest_common_ancestor(["a", "c"])
        self.assertIs(obs, t_root)

        # empty case
        with self.assertRaises(ValueError):
            t1.lowest_common_ancestor([])

    # ------------------------------------------------
    # Tree traversal
    # ------------------------------------------------

    def test_preorder(self):
        """Preorder traversal of the tree"""
        exp = ["root", "i1", "a", "b", "i2", "c", "d"]
        obs = [n.name for n in self.simple_t.preorder()]
        self.assertEqual(obs, exp)

        exp = ["i1", "a", "b", "i2", "c", "d"]
        obs = [n.name for n in self.simple_t.preorder(include_self=False)]
        self.assertEqual(obs, exp)

    def test_postorder(self):
        """Postorder traversal of the tree"""
        exp = ["a", "b", "i1", "c", "d", "i2", "root"]
        obs = [n.name for n in self.simple_t.postorder()]
        self.assertEqual(obs, exp)

        exp = ["a", "b", "i1", "c", "d", "i2"]
        obs = [n.name for n in self.simple_t.postorder(include_self=False)]
        self.assertEqual(obs, exp)

    def test_pre_and_postorder(self):
        """Pre and post order traversal of the tree"""
        exp = ["root", "i1", "a", "b", "i1", "i2", "c", "d", "i2", "root"]
        obs = [n.name for n in self.simple_t.pre_and_postorder()]
        self.assertEqual(obs, exp)
        obs2 = [n.name for n in self.simple_t.traverse(True, True)]
        self.assertEqual(obs2, exp)

    def test_pre_and_postorder_no_children(self):
        t = TreeNode("brofist")

        # include self
        exp = ["brofist"]
        obs = [n.name for n in t.pre_and_postorder()]
        self.assertEqual(obs, exp)

        # do not include self
        obs = list(t.pre_and_postorder(include_self=False))
        self.assertEqual(obs, [])

    def test_levelorder(self):
        """Level order traversal of the tree"""
        exp = ["root", "i1", "i2", "a", "b", "c", "d"]
        obs = [n.name for n in self.simple_t.levelorder()]
        self.assertEqual(obs, exp)

        exp = ["i1", "i2", "a", "b", "c", "d"]
        obs = [n.name for n in self.simple_t.levelorder(include_self=False)]
        self.assertEqual(obs, exp)

    def test_tips(self):
        """Tip traversal of tree"""
        exp = ["a", "b", "c", "d"]
        obs = [n.name for n in self.simple_t.tips()]
        self.assertEqual(obs, exp)
        obs2 = [n.name for n in self.simple_t.traverse(False, False)]
        self.assertEqual(obs2, exp)

    def test_tips_self(self):
        """ See issue #1509 """
        tree = TreeNode.read(["(c,(b,a)x)y;"])
        ts = list(tree.find("c").tips(include_self=True))
        self.assertEqual(len(ts), 1)
        t = ts[0]
        self.assertEqual(t.name, "c")
        self.assertTrue(t.is_tip())

    # ------------------------------------------------
    # Tree manipulation
    # ------------------------------------------------

    def test_append(self):
        """Add a node to children."""
        t = self.simple_t

        # append a single tip
        t.create_caches()
        t.append(TreeNode(name="n1"))
        self.assertEqual(len((c := t.children)), 3)
        self.assertEqual(c[-1].name, "n1")
        self.assertIs(c[-1].parent, t)
        self.assertFalse(hasattr(t, "_tip_cache"))

        # append an entire tree
        t2 = TreeNode.read(["(x,y)z;"])
        t.append(t2)
        self.assertEqual(len((c := t.children)), 4)
        self.assertEqual(c[0].name, "i1")
        self.assertEqual(c[1].name, "i2")
        self.assertEqual(c[2].name, "n1")
        self.assertEqual(c[3].name, "z")
        self.assertEqual(c[3].children[0].name, "x")
        self.assertEqual(c[3].children[1].name, "y")
        self.assertEqual(t2.parent, t)

        # move a clade from another tree, keep cache
        t3 = TreeNode.read(["((x2,x3)o,(x4,x5)p)q;"])
        n3 = t3.find("o")
        t2.create_caches()
        t2.append(n3, uncache=False)
        self.assertEqual(len((c2 := t2.children)), 3)
        self.assertEqual(c2[0].name, "x")
        self.assertEqual(c2[1].name, "y")
        self.assertEqual(c2[2].name, "o")
        self.assertEqual(c2[2].children[0].name, "x2")
        self.assertEqual(c2[2].children[1].name, "x3")
        self.assertIs(c2[2].parent, t2)
        self.assertEqual(len(t3.children), 1)
        self.assertEqual(t3.children[0].name, "p")
        self.assertTrue(hasattr(t, "_tip_cache"))
        self.assertTrue(hasattr(t3, "_tip_cache"))

    def test_extend(self):
        """Add arbitrary number of nodes to children."""

        # add two clades
        t = self.simple_t
        t.create_caches()
        t2 = TreeNode.read(["(x1,y1)z1;"])
        t3 = TreeNode.read(["(x2,y2)z2;"])
        t.extend([t2, t3])
        self.assertIs(t2.parent, t)
        self.assertIs(t3.parent, t)
        self.assertEqual(len(c := t.children), 4)
        self.assertEqual(c[0].name, "i1")
        self.assertEqual(c[1].name, "i2")
        self.assertEqual(c[2].name, "z1")
        self.assertEqual(c[3].name, "z2")
        self.assertEqual(c[2].children[0].name, "x1")
        self.assertEqual(c[2].children[1].name, "y1")
        self.assertEqual(c[3].children[0].name, "x2")
        self.assertEqual(c[3].children[1].name, "y2")
        self.assertFalse(hasattr(t, "_tip_cache"))

        # move all children, keep cache
        t1 = TreeNode.read(["(x1,y1)z1;"])
        t4 = TreeNode.read(["(x2,y2)z2;"])
        t1.create_caches()
        t4.create_caches()
        t1.extend(t4.children, uncache=False)
        self.assertEqual(len(t4.children), 0)
        self.assertEqual(len(c := t1.children), 4)
        self.assertEqual(c[0].name, "x1")
        self.assertEqual(c[1].name, "y1")
        self.assertEqual(c[2].name, "x2")
        self.assertEqual(c[3].name, "y2")
        self.assertIs(c[2].parent, t1)
        self.assertIs(c[3].parent, t1)
        self.assertTrue(hasattr(t1, "_tip_cache"))
        self.assertTrue(hasattr(t4, "_tip_cache"))

        # empty input
        t.extend([])
        self.assertEqual(len(t.children), 4)

    def test_insert(self):
        "Insert a node into the branch connecting self and its parent."
        # insert a new node into a branch with no length
        t = self.simple_t
        node = t.find("i1")
        node.insert(TreeNode("x"))
        obs = t.find("x")
        self.assertTrue(obs.parent is t)
        self.assertTrue(node.parent is obs)
        self.assertIn(obs, t.children)
        self.assertIn(node, obs.children)
        self.assertIsNone(obs.length)
        self.assertIsNone(node.length)

        msg = "Distance is provided but branch has no length."
        with self.assertRaisesRegex(ValueError, msg):
            node.insert(TreeNode("x"), distance=1.0)

        msg = "Self has no parent."
        with self.assertRaisesRegex(NoParentError, msg):
            t.insert(TreeNode("x"))

        # insert an existing clade into a branch with length
        t = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5,g:1)h;"])
        donor_t = TreeNode.read(["((x:1,y:1)m:1.5,(z:1,w:1)n:0.5,l:2.5);"])
        t.find("c").insert(donor_t.find("m"))
        obs = t.find("m")
        self.assertTrue(obs.parent is t)
        self.assertTrue(t.find("c").parent is obs)
        self.assertNotIn(obs, donor_t.children)
        self.assertEqual(obs.length, 1)
        self.assertEqual(t.find("c").length, 1)

        t.find("d").insert(donor_t.find("n"), 2)
        obs = t.find("n")
        self.assertTrue(obs.parent is t.find("f"))
        self.assertTrue(t.find("d").parent is obs)
        self.assertEqual(obs.length, 1)
        self.assertEqual(t.find("d").length, 2)

        msg = "Distance cannot exceed branch length."
        with self.assertRaisesRegex(ValueError, msg):
            t.find("c").insert(TreeNode("x"), 20)

        # with branch support, keep cache
        t = TreeNode.read(["(((a,b)90)d);"])
        t.assign_supports()
        t.create_caches()
        n = t.lca(["a", "b"])
        t.lca(["a", "b"]).insert(TreeNode("x"), uncache=False)
        self.assertTrue(hasattr(t, "_tip_cache"))
        self.assertEqual(n.parent.name, "x")
        self.assertEqual(n.parent.parent.name, "d")
        self.assertIs(n.parent.support, 90)

        # with custom branch attribute
        t = TreeNode.read(["(((a,b)c)d);"])
        n = t.find("c")
        n.battr = 1  # branch attribute
        n.nattr = 2  # node attribute
        n.insert(TreeNode("x"), branch_attrs=["battr"])
        self.assertEqual(t.find("x").battr, 1)
        self.assertFalse(hasattr(t.find("x"), "nattr"))

    def test_pop(self):
        """Pop off a node by index."""
        t = self.simple_t
        t.extend([
            TreeNode.read(["(x1,y1)z1;"]),
            TreeNode.read(["(x2,y2)z2;"])
        ])

        # pop last child (default)
        t.create_caches()
        z2 = t.pop()
        self.assertEqual(z2.name, "z2")
        self.assertIsNone(z2.parent)
        self.assertEqual(len(t.children), 3)
        self.assertNotIn(z2, t.children)
        self.assertEqual(z2.children[0].name, "x2")
        self.assertEqual(z2.children[1].name, "y2")
        self.assertFalse(hasattr(t, "_tip_cache"))

        # pop first child, keep cache
        t.create_caches()
        i1 = t.pop(0, uncache=False)
        self.assertEqual(i1.name, "i1")
        self.assertIsNone(i1.parent)
        self.assertEqual(len(t.children), 2)
        self.assertNotIn(i1, t.children)
        self.assertEqual(i1.children[0].name, "a")
        self.assertEqual(i1.children[1].name, "b")
        self.assertTrue(hasattr(t, "_tip_cache"))

        # check remaining tree structure
        self.assertEqual(t.children[0].name, "i2")
        self.assertEqual(t.children[1].name, "z1")

    def test_remove(self):
        """Remove a node."""
        t = self.simple_t

        # a node that can be removed
        parent = t.find("i1")
        child = t.find("a")
        self.assertTrue(parent.remove(child))
        self.assertIsNone(child.parent)
        self.assertEqual(len(parent.children), 1)
        self.assertNotIn(child, parent.children)
        self.assertEqual(parent.children[0].name, "b")
        self.assertFalse(hasattr(t, "_tip_cache"))

        # nodes that cannot be removed
        self.assertFalse(t.remove(child))
        self.assertFalse(t.remove(t.find("c")))
        self.assertFalse(t.remove(TreeNode()))

        # keep cache
        t.create_caches()
        self.assertTrue(t.remove(parent, uncache=False))
        self.assertTrue(hasattr(t, "_tip_cache"))

    def test_remove_by_func(self):
        """Remove nodes by function."""
        t = self.simple_t

        def f(node):
            return node.name in ("b", "d")

        # def function
        t.remove_by_func(f)
        self.assertEqual(str(t), "((a)i1,(c)i2)root;\n")

        # lambda function, keep cache
        t.create_caches()
        t.remove_by_func(lambda x: x.name == "a", uncache=False)
        self.assertEqual(str(t), "(i1,(c)i2)root;\n")
        self.assertTrue(hasattr(t, "_tip_cache"))

    def test_prune(self):
        """Collapse single-child nodes."""
        t = self.simple_t.copy()

        # nothing to collapse because tree is bifurcating
        cp = t.copy()
        t.prune()
        for obs, exp in zip(t.traverse(), cp.traverse()):
            self.assertEqual(obs.name, exp.name)
            self.assertEqual(obs.length, exp.length)

        # create a single child by removing tip "a"
        n = t.children[0]
        n.remove(n.children[0])
        t.prune()
        self.assertEqual(len(t.children), 2)
        self.assertEqual(t.children[0].name, "i2")
        self.assertEqual(t.children[1].name, "b")
        self.assertIsNone(t.children[1].length)

        # tree with branch lengths
        t = self.simple_t.copy()
        for n in t.traverse():
            n.length = 1
        n = t.children[0]
        n.remove(n.children[0])
        t.prune()
        self.assertEqual(len(t.children), 2)
        self.assertEqual(t.children[0].name, "i2")
        self.assertEqual(t.children[1].name, "b")
        self.assertAlmostEqual(t.children[1].length, 2)

        # root has single child
        t = TreeNode.read(["((a,b)c)extra;"])
        t.prune()
        self.assertEqual(str(t), "(a,b)c;\n")

        # nested single-child nodes
        t = TreeNode.read(["(((a,b)));"])
        t.prune()
        self.assertEqual(str(t), "(a,b);\n")

        # single-child root and branch lengths
        t = TreeNode.read(["((((a:1)b:2)c:3)d:4)e:5;"])
        t.prune()
        self.assertEqual(t.name, "a")
        self.assertAlmostEqual(t.length, 15)

        # keep cache
        t = TreeNode.read(["((((a)b)c)d)e;"])
        t.create_caches()
        t.prune(uncache=False)
        self.assertEqual(t.name, "a")
        self.assertTrue(hasattr(t, "_tip_cache"))

        # prune from internal node
        t = TreeNode.read(["((((a:1)b:2)c:3)d:4)e:5;"])
        n = t.children[0].children[0]
        n.prune()
        self.assertEqual(t.name, "e")
        self.assertEqual(len(t.children), 1)
        c = t.children[0]
        self.assertEqual(c.name, "d")
        self.assertAlmostEqual(c.length, 4)
        self.assertEqual(len(c.children), 1)
        self.assertIs(c.children[0], n)
        self.assertIs(n.parent, c)
        self.assertAlmostEqual(n.length, 6)

    def test_shear(self):
        """Shear tree to keep given tips."""
        # LCA is root, and root is retained
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        obs = str(t.shear(["G", "M"]))
        exp = "(G:3.0,M:3.7);\n"
        self.assertEqual(obs, exp)

        # in place
        obs = t.shear(["G", "M"], inplace=True)
        self.assertEqual(str(t), exp)
        self.assertIsNone(obs)

        # keep cache
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        t.create_caches()
        obs = t.shear(["G", "M"], inplace=True, uncache=False)
        self.assertEqual(str(t), exp)
        self.assertTrue(hasattr(t, "_tip_cache"))

        # LCA is internal, and root is dropped
        t = TreeNode.read(["((a,b),((c,d),(e,f)));"])
        obs = str(t.shear(["c", "d"]))
        exp = "(c,d);\n"
        self.assertEqual(obs, exp)

        # don't prune
        obs = str(t.shear(["c", "d"], prune=False))
        exp = "(((c,d)));\n"
        self.assertEqual(obs, exp)

        # in place
        exp = "(c,d);\n"
        t.shear(["c", "d"], inplace=True)
        self.assertEqual(str(t), exp)

        # polytomy (issue 1416)
        t = TreeNode.read(["(((a,b,f,g),c),d);"])
        obs = str(t.shear(["a", "b", "c", "f"]))
        exp = "((a,b,f),c);\n"
        self.assertEqual(obs, exp)

        # in place
        t.shear(["a", "b", "c", "f"], inplace=True)
        self.assertEqual(str(obs), exp)

        # complex example; shear from internal node
        t = TreeNode.read(["((((((a,b),(c,d)),(e,f)):1,(g,(h,i))):2,((j,k),l)),m);"])
        n = t.children[0].children[0]

        # LCA is root of subtree
        names = {"a", "b", "e", "h"}
        obs = n.shear(names)
        self.assertAlmostEqual(obs.length, 2)
        self.assertSetEqual(obs.subset(), names)

        # LCA is internal node of subtree (check branch length)
        names = {"a", "b", "e"}
        obs = n.shear(names)
        self.assertAlmostEqual(obs.length, 3)
        self.assertSetEqual(obs.subset(), names)

        # in place (check branch length and connection)
        n.shear(names, inplace=True)
        self.assertAlmostEqual(n.length, 3)
        self.assertIs(n.parent, t.children[0])
        self.assertSetEqual(n.subset(), names)

        # name not found
        msg = "Names are not a subset of the tree."
        with self.assertRaisesRegex(ValueError, msg):
            TreeNode.read(["(a,b)c;"]).shear(["a", "c"])

        # non-strict mode
        obs = str(TreeNode.read(["(a,b)c;"]).shear(["a", "c"], strict=False))
        exp = "a;\n"
        self.assertEqual(obs, exp)

    def test_unpack(self):
        """Unpack an internal node."""
        # test unpacking a node without branch length
        tree = TreeNode.read(["((c,d)a,(e,f)b);"])
        tree.find("b").unpack()
        exp = "((c,d)a,e,f);\n"
        self.assertEqual(str(tree), exp)
        self.assertFalse(hasattr(tree, "_tip_cache"))

        # keep cache
        tree = TreeNode.read(["((c,d)a,(e,f)b);"])
        tree.find("b").unpack(uncache=False)
        self.assertEqual(str(tree), exp)
        self.assertTrue(hasattr(tree, "_tip_cache"))

        # test unpacking a node with branch length
        tree = TreeNode.read(["((c:2.0,d:3.0)a:1.0,(e:2.0,f:1.0)b:2.0);"])
        tree.find("b").unpack()
        exp = "((c:2.0,d:3.0)a:1.0,e:4.0,f:3.0);"
        self.assertEqual(str(tree).rstrip(), exp)

        # test attempting to unpack root
        tree = TreeNode.read(["((d,e)b,(f,g)c)a;"])
        msg = "Cannot unpack root."
        with self.assertRaisesRegex(TreeError, msg):
            tree.find("a").unpack()

        # test attempting to unpack tip
        msg = "Cannot unpack tip."
        with self.assertRaisesRegex(TreeError, msg):
            tree.find("d").unpack()

    def test_unpack_by_func(self):
        """Unpack internal nodes of a tree by a function."""
        # unpack internal nodes with branch length <= 1.0
        def func(x):
            return x.length <= 1.0

        # will unpack node "a", but not tip "e"
        # will add the branch length of "a" to its child nodes "c" and "d"
        tree = TreeNode.read(["((c:2,d:3)a:1,(e:1,f:2)b:2);"])
        tree.unpack_by_func(func)
        exp = "((e:1.0,f:2.0)b:2.0,c:3.0,d:4.0);"
        self.assertEqual(str(tree).rstrip(), exp)

        # keep cache
        tree = TreeNode.read(["((c:2,d:3)a:1,(e:1,f:2)b:2);"])
        tree.create_caches()
        tree.unpack_by_func(func, uncache=False)
        self.assertEqual(str(tree).rstrip(), exp)
        self.assertTrue(hasattr(tree, "_tip_cache"))

        # unpack internal nodes with branch length < 2.01
        # will unpack both "a" and "b"
        tree = TreeNode.read(["((c:2,d:3)a:1,(e:1,f:2)b:2);"])
        tree.unpack_by_func(lambda x: x.length <= 2.0)
        exp = "(c:3.0,d:4.0,e:3.0,f:4.0);"
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack two nested nodes "a" and "c" simultaneously
        tree = TreeNode.read(["(((e:3,f:2)c:1,d:3)a:1,b:4);"])
        tree.unpack_by_func(lambda x: x.length <= 2.0)
        exp = "(b:4.0,d:4.0,e:5.0,f:4.0);"
        self.assertEqual(str(tree).rstrip(), exp)

        # test a complicated scenario (unpacking nodes "g", "h" and "m")
        def func(x):
            return x.length < 2.0
        tree = TreeNode.read(["(((a:1.04,b:2.32,c:1.44)d:3.20,"
                              "(e:3.91,f:2.47)g:1.21)h:1.75,"
                              "(i:4.14,(j:2.06,k:1.58)l:3.32)m:0.77);"])
        tree.unpack_by_func(func)
        exp = ("((a:1.04,b:2.32,c:1.44)d:4.95,e:6.87,f:5.43,i:4.91,"
               "(j:2.06,k:1.58)l:4.09);")
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 75
        def func(x):
            return x.support < 75
        tree = TreeNode.read(["(((a,b)85,(c,d)78)75,(e,(f,g)64)80);"])
        tree.assign_supports()
        tree.unpack_by_func(func)
        exp = "(((a,b)85,(c,d)78)75,(e,f,g)80);"
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 85
        tree = TreeNode.read(["(((a,b)85,(c,d)78)75,(e,(f,g)64)80);"])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support < 85)
        exp = "((a,b)85,c,d,e,f,g);"
        self.assertEqual(str(tree).rstrip(), exp)

        # unpack nodes with support < 0.95
        tree = TreeNode.read(["(((a,b)0.97,(c,d)0.98)1.0,(e,(f,g)0.88)0.96);"])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support < 0.95)
        exp = "(((a,b)0.97,(c,d)0.98)1.0,(e,f,g)0.96);"
        self.assertEqual(str(tree).rstrip(), exp)

        # test a case where there are branch lengths, none support values and
        # node labels
        tree = TreeNode.read(["(((a:1.02,b:0.33)85:0.12,(c:0.86,d:2.23)"
                              "70:3.02)75:0.95,(e:1.43,(f:1.69,g:1.92)64:0.20)"
                              "node:0.35)root;"])
        tree.assign_supports()
        tree.unpack_by_func(lambda x: x.support is not None and x.support < 75)
        exp = ("(((a:1.02,b:0.33)85:0.12,c:3.88,d:5.25)75:0.95,"
               "(e:1.43,f:1.89,g:2.12)node:0.35)root;")
        self.assertEqual(str(tree).rstrip(), exp)

    def test_bifurcate(self):
        t1 = TreeNode.read(["(((a,b),c),(d,e));"])
        t1.create_caches()
        t1.bifurcate()
        self.assertEqual(str(t1), "(((a,b),c),(d,e));\n")
        self.assertFalse(hasattr(t1, "_tip_cache"))

        # keep cache
        t1 = TreeNode.read(["(((a,b),c),(d,e));"])
        t1.create_caches()
        t1.bifurcate(uncache=False)
        self.assertEqual(str(t1), "(((a,b),c),(d,e));\n")
        self.assertTrue(hasattr(t1, "_tip_cache"))

        # with and without insert length
        t2 = TreeNode.read(["((a,b,c));"])
        t3 = t2.copy()
        t2.bifurcate()
        self.assertEqual(str(t2), "((c,(a,b)));\n")
        t3.bifurcate(insert_length=0)
        self.assertEqual(str(t3), "((c,(a,b):0));\n")
        
        # bifurcate with subclass
        tree = TreeNodeSubclass()
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())
        tree.append(TreeNodeSubclass())

        tree.bifurcate()

        for node in tree.traverse():
            self.assertIs(type(node), TreeNodeSubclass)

    def test_shuffle(self):
        # default behavior: all tips are shuffled, only one tree is yielded
        # shuffling method is stochastic and so is the result
        t = self.simple_t.copy()
        obs = list(t.shuffle())
        self.assertEqual(len(obs), 1)
        self.assertSetEqual(obs[0].subset(), set('abcd'))

        # specify a random generator to make results deterministic
        rng = np.random.default_rng(42)
        t = self.simple_t.copy()
        obs = str(next(t.shuffle(shuffle_f=rng)))
        exp = "((d,c)i1,(b,a)i2)root;\n"
        self.assertEqual(obs, exp)

        # can also specify a random seed; result is the same
        t = self.simple_t.copy()
        obs = str(next(t.shuffle(shuffle_f=42)))
        self.assertEqual(obs, exp)

        # can also supply a function; result is the same
        t = self.simple_t.copy()
        f = np.random.default_rng(42).shuffle
        obs = str(next(t.shuffle(shuffle_f=f)))
        self.assertEqual(obs, exp)

        # yield a row of 5 trees
        rng = np.random.default_rng(42)
        t = self.simple_t.copy()
        obs = list(map(str, t.shuffle(shuffle_f=rng, n=5)))
        self.assertEqual(len(obs), 5)
        exp = ["((d,c)i1,(b,a)i2)root;\n",
               "((a,b)i1,(d,c)i2)root;\n",
               "((a,c)i1,(d,b)i2)root;\n",
               "((d,b)i1,(a,c)i2)root;\n",
               "((a,c)i1,(d,b)i2)root;\n"]
        self.assertListEqual(obs, exp)

        # yield infinitely
        rng = np.random.default_rng(42)
        t = self.simple_t.copy()
        gen = t.shuffle(shuffle_f=rng, n=None)
        obs = [str(next(gen)) for i in range(100)]
        self.assertListEqual(obs[:5], exp)

        # define two simple, non-stochastic shuffling functions
        def rev_f(items):
            items.reverse()

        def rotate_f(items):
            tmp = items[-1]
            items[1:] = items[:-1]
            items[0] = tmp

        # apply a simple function
        t = self.simple_t.copy()
        obs = str(next(t.shuffle(shuffle_f=rev_f)))
        exp = "((d,c)i1,(b,a)i2)root;\n"
        self.assertEqual(obs, exp)

        # specify names to shuffle
        t = self.simple_t.copy()
        obs = list(map(str, t.shuffle(names=list("abc"), shuffle_f=rotate_f, n=4)))
        exp = ["((c,a)i1,(b,d)i2)root;\n",
               "((b,c)i1,(a,d)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((c,a)i1,(b,d)i2)root;\n"]
        self.assertListEqual(obs, exp)

        # specify number of names to shuffle
        t = self.simple_t.copy()
        obs = list(map(str, t.shuffle(k=2, shuffle_f=rev_f, n=5)))
        exp = ["((a,b)i1,(d,c)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((a,b)i1,(d,c)i2)root;\n",
               "((a,b)i1,(c,d)i2)root;\n",
               "((a,b)i1,(d,c)i2)root;\n"]
        self.assertListEqual(obs, exp)

        # a complex example
        obs = list(map(str, self.complex_tree.shuffle(
            shuffle_f=rev_f, names=["c", "d", "e", "f"], n=4)))
        exp = ["(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);\n",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);\n"]
        self.assertListEqual(obs, exp)

        # invalid number of iterations
        t = self.simple_t.copy()
        with self.assertRaises(ValueError):
            next(t.shuffle(n=-1))

        # invalid number of names to shuffle
        with self.assertRaises(ValueError):
            next(t.shuffle(k=1))

        # k and names conflict
        with self.assertRaises(ValueError):
            next(t.shuffle(k=5, names=["a", "b"]))

        # tip names not found
        with self.assertRaises(MissingNodeError):
            next(t.shuffle(names=["x", "y"]))

    # ------------------------------------------------
    # Tree rerooting
    # ------------------------------------------------

    def test_unroot(self):
        """Convert a rooted tree into unrooted."""
        # default behavior
        t = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        t.unroot()
        exp = "(a,b,(d,e)f)c;\n"
        self.assertEqual(str(t), exp)

        # choose the other side
        t = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        t.unroot(side=1)
        exp = "((a,b)c,d,e)f;\n"
        self.assertEqual(str(t), exp)

        # start from internal node
        t = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        t.find("c").unroot()
        exp = "(a,b,(d,e)f)c;\n"
        self.assertEqual(str(t), exp)
        self.assertFalse(hasattr(t, "_tip_cache"))

        # keep caches
        t = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        t.find("c").unroot(uncache=False)
        exp = "(a,b,(d,e)f)c;\n"
        self.assertEqual(str(t), exp)
        self.assertTrue(hasattr(t, "_tip_cache"))

        # with branch lengths
        t = TreeNode.read(["((a:2.0,b:1.5)c:0.5,(d:1.0,e:1.2)f:0.3)g;"])
        t.unroot()
        exp = "(a:2.0,b:1.5,(d:1.0,e:1.2)f:0.8)c;\n"
        self.assertEqual(str(t), exp)

        # other child has no branch length
        t = TreeNode.read(["((a,b)c:1.0,(d,e)f)g;"])
        t.unroot()
        exp = "(a,b,(d,e)f:1.0)c;\n"
        self.assertEqual(str(t), exp)

        # first child is a tip
        t = TreeNode.read(["(a,(b,c)d)e;"])
        t.unroot()
        exp = "(a,b,c)d;\n"
        self.assertEqual(str(t), exp)

        # both children are tips
        t = TreeNode.read(["(a,b)c;"])
        t.unroot()
        exp = "(b)a;\n"
        self.assertEqual(str(t), exp)

        # tree is already unrooted
        t = TreeNode.read(["(a,b,(d,e)f)c;"])
        t.unroot()
        exp = "(a,b,(d,e)f)c;\n"
        self.assertEqual(str(t), exp)

    def test_unrooted_copy(self):
        tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        node = tree.find("d")

        # name as branch label (default behavior, but will change in the
        # future)
        obs = node.unrooted_copy()
        exp = "(b,c,(a,((f,g)h)e)d)root;\n"
        self.assertEqual(str(obs), exp)

        # name as node label
        obs = node.unrooted_copy(branch_attrs={"length"})
        exp = "(b,c,(a,((f,g)h)i)e)d;\n"
        self.assertEqual(str(obs), exp)

        # name the new root node (only when it doesn't have one)
        obs = node.unrooted_copy(root_name="hello")
        exp = "(b,c,(a,((f,g)h)e)d)hello;\n"
        self.assertEqual(str(obs), exp)

        obs = node.unrooted_copy(branch_attrs={"length"}, root_name="hello")
        exp = "(b,c,(a,((f,g)h)i)e)d;\n"
        self.assertEqual(str(obs), exp)

        # tree caches are cleaned
        tree.cache_attr(lambda n: 1, "node_count", sum)
        self.assertTrue(hasattr(node, "node_count"))
        obs = node.unrooted_copy()
        self.assertFalse(hasattr(obs, "node_count"))

        # transfer branch support to opposite node
        tree = TreeNode.read(["((a,b)90,(c,d)90);"])
        node = tree.find("a")
        obs = node.unrooted_copy(branch_attrs={"support", "length"})
        exp = "((b,((c,d)90))90)a;\n"
        self.assertEqual(str(obs), exp)

        tree.assign_supports()
        obs = node.unrooted_copy(branch_attrs={"support", "length"})
        exp = "((b,((c,d)90)90))a;\n"
        self.assertEqual(str(obs), exp)

        # retain custom attributes
        tree = TreeNode.read(["(((a,b)c,d)e,f)g;"])
        tree.find("c").dummy = "this"
        tree.find("e").dummy = "that"
        obs = tree.find("c").unrooted_copy(branch_attrs={"length"})
        exp = "(a,b,(d,(f)g)e)c;\n"
        self.assertEqual(str(obs), exp)
        self.assertEqual(obs.dummy, "this")
        self.assertEqual(obs.find("e").dummy, "that")
        self.assertIsNone(getattr(obs.find("d"), "dummy", None))

        # deep vs shallow copy
        tree = TreeNode.read(["(((a,b)c,d)e,f)g;"])
        tree.find("c").dummy = [1, [2, 3], 4]
        tcopy = tree.unrooted_copy(deep=True)
        tcopy.find("c").dummy[1].append(0)
        self.assertListEqual(tree.find("c").dummy[1], [2, 3])

        tcopy = tree.unrooted_copy()
        tcopy.find("c").dummy[1].append(0)
        self.assertListEqual(tree.find("c").dummy[1], [2, 3, 0])

    def test_unrooted_deepcopy(self):
        t = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        exp = "(b,c,(a,((f,g)h)e)d)root;\n"
        obs = t.find("d").unrooted_deepcopy()
        self.assertEqual(str(obs), exp)

        t_ids = {id(n) for n in t.traverse()}
        obs_ids = {id(n) for n in obs.traverse()}

        self.assertEqual(t_ids.intersection(obs_ids), set())

    def test_unrooted_move(self):
        t = TreeNode.read(["(((a:1,b:1)c:1,(d:1,e:1)f:2)g:0.5,(h:1,i:1)j:0.5)k;"])
        tcopy = t.copy()
        obs = tcopy.find("c")
        obs.unrooted_move()
        self.assertFalse(hasattr(tcopy, "_tip_cache"))
        exp = TreeNode.read(["(a:1,b:1,((d:1,e:1)f:2,((h:1,i:1)j:0.5)k:0.5)g:1)c;"])
        self.assertTrue(obs.is_root())
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

        # keep caches
        tcopy = t.copy()
        obs = tcopy.find("c")
        obs.unrooted_move(uncache=False)
        self.assertTrue(hasattr(tcopy, "_tip_cache"))

    def test_root_at(self):
        """Root tree at a given node."""
        t = TreeNode.read(["(((a,b)c,(d,e)f)g,h)i;"])

        # original behavior (name as branch label); deprecated
        obs = str(t.root_at("c"))
        exp = "(a,b,((d,e)f,(h)g)c)root;\n"
        self.assertEqual(obs, exp)

        # root at internal node
        obs = str(t.root_at("c", branch_attrs=[]))
        exp = "(a,b,((d,e)f,(h)i)g)c;\n"
        self.assertEqual(obs, exp)

        # root at self
        obs = str(t.find("c").root_at(branch_attrs=[]))
        self.assertEqual(obs, exp)

        # root at tip (and input node instead of name)
        obs = str(t.root_at(t.find("h"), branch_attrs=[]))
        exp = "((((a,b)c,(d,e)f)g)i)h;\n"
        self.assertEqual(obs, exp)

        # root at root (no change)
        obs = str(t.root_at("i", branch_attrs=[]))
        self.assertEqual(obs, str(t))

        # in-place rooting
        n = t.copy().find("c")
        obs = n.root_at(inplace=True)
        exp = "(a,b,((d,e)f,(h)g)c)root;\n"
        self.assertIs(n, obs)
        self.assertEqual(str(obs), exp)

    def test_root_at_above(self):
        """Root tree at the branch above a given node."""
        # no branch length
        t = TreeNode.read(["(((a,b)c,(d,e)f)g,h)i;"])
        obs = str(t.root_at("c", above=True, branch_attrs=[]))
        exp = "((a,b)c,((d,e)f,(h)i)g)root;\n"
        self.assertEqual(obs, exp)

        # in-place rooting
        n = t.find("c")
        obs = t.root_at(n, above=True, branch_attrs=[], inplace=True)
        self.assertIs(n.parent, obs)
        self.assertEqual(str(obs), exp)

        # root at midpoint of branch
        t = TreeNode.read(["(((a,b)c:1.0,(d,e)f)g,h)i;"])
        obs = str(t.root_at("c", above=True, branch_attrs=[]))
        exp = "((a,b)c:0.5,((d,e)f,(h)i)g:0.5)root;\n"
        self.assertEqual(obs, exp)

        # root at specific position
        t = TreeNode.read(["(((a,b)c:1.0,(d,e)f)g,h)i;"])
        obs = str(t.root_at("c", above=0.4, branch_attrs=[]))
        exp = "((a,b)c:0.4,((d,e)f,(h)i)g:0.6)root;\n"
        self.assertEqual(obs, exp)

        # with branch support
        t = TreeNode.read(["(((a,b)'90:c',(d,e)'80:f')g,h)i;"])
        t.assign_supports()
        obs = str(t.root_at("c", above=True, branch_attrs=[]))
        exp = "((a,b)'90:c',((d,e)'80:f',(h)i)'90:g')root;\n"
        self.assertEqual(obs, exp)

    def test_root_at_reset(self):
        """Root tree while resetting original root."""
        t = TreeNode.read(["(((a,b)c,(d,e)f)g,h)i;"])

        # unroot tree prior to rerooting
        obs = str(t.root_at("c", reset=True, branch_attrs=[]))
        exp = "(a,b,((d,e)f,h)g)c;\n"
        self.assertEqual(obs, exp)

        # root at a basal node (which will be avoided during unrooting)
        obs = str(t.root_at("g", reset=True, branch_attrs=[]))
        exp = "((a,b)c,(d,e)f,h)g;\n"
        self.assertEqual(obs, exp)

        # in-place rooting
        n = t.find("g")
        obs = t.root_at(n, reset=True, branch_attrs=[], inplace=True)
        self.assertIs(n, obs)
        self.assertEqual(str(obs), exp)

        # tree is already unrooted
        t = TreeNode.read(["((a,b)c,d,e)f;"])
        obs = str(t.root_at("c", branch_attrs=[], reset=True))
        exp = str(t.root_at("c", branch_attrs=[]))
        self.assertEqual(obs, exp)

    def test_root_at_midpoint(self):
        """Root tree at the midpoint"""
        t = self.TreeRoot
        for n in t.traverse():
            n.length = 1

        # g and h are farthest apart, by 5, therefore root should be
        # 2.5 away from h, i.e., midpoint between b and c
        result = t.root_at_midpoint()
        self.assertEqual(result.distance(result.find("e")), 1.5)
        self.assertEqual(result.distance(result.find("g")), 2.5)
        exp_dist = t.tip_tip_distances()
        obs_dist = result.tip_tip_distances()
        self.assertEqual(obs_dist, exp_dist)

        # in-place rerooting
        b = t.find("b")
        result = t.root_at_midpoint(inplace=True)
        self.assertIs(b.parent, result)
        self.assertEqual(result.tip_tip_distances(), exp_dist)

    def test_root_at_midpoint_no_lengths(self):
        # should get same tree back (a copy)
        nwk = "(a,b)c;\n"
        t = TreeNode.read([nwk])
        obs = t.root_at_midpoint()
        self.assertEqual(str(obs), nwk)

    def test_root_at_midpoint_tie(self):
        # original behavior (name as branch label); deprecated
        t = TreeNode.read(["(((a:1,b:1)c:2,(d:3,e:4)f:5),g:1)root;"])
        obs = t.root_at_midpoint()
        exp = TreeNode.read(["((d:3,e:4)f:2,((a:1,b:1)c:2,(g:1)):3)root;"])
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

        t = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5,g:1)h;"])
        # farthest tip-to-tip distance is 12 (a or b to e)
        # therefore new root should be 2 above f
        obs = t.root_at_midpoint(branch_attrs=[])
        exp = TreeNode.read(["((d:3,e:4)f:2,((a:1,b:1)c:2,g:1)h:3)root;"])
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

        # no root name
        obs = t.root_at_midpoint(branch_attrs=[], root_name=None)
        self.assertIsNone(obs.name)

        # with branch support
        t = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)'80:f':5,g:1)h;"])
        t.assign_supports()
        obs = t.root_at_midpoint(branch_attrs=[])
        exp = TreeNode.read(["((d:3,e:4)'80:f':2,((a:1,b:1)c:2,g:1)'80:h':3)root;"])
        exp.assign_supports()
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)
            self.assertEqual(o.support, e.support)

    def test_root_at_midpoint_node(self):
        t = TreeNode.read(["(((a:2,b:3)c:1,d:1)e:1,f:3)g;"])
        # farthest tip-to-tip distance is 8 (b - c - e - f)
        # therefore new root should be at e
        obs = t.root_at_midpoint(branch_attrs=[])
        exp = TreeNode.read(["((a:2.0,b:3.0)c:1.0,d:1.0,(f:3.0)g:1.0)e;"])
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

        # remove original root
        obs = t.root_at_midpoint(branch_attrs=[], reset=True)
        exp = TreeNode.read(["((a:2.0,b:3.0)c:1.0,d:1.0,f:4.0)e;"])
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.name, e.name)
            self.assertEqual(o.length, e.length)

    def test_root_by_outgroup(self):
        tree = TreeNode.read(["((((a,b),(c,d)),(e,f)),g);"])

        # outgroup is monophyletic
        obs = str(tree.root_by_outgroup(["a", "b"]))
        exp = "((a,b),((c,d),((e,f),g)));\n"
        self.assertEqual(obs, exp)

        # outgroup is monophyletic after rotating
        obs = str(tree.root_by_outgroup(["e", "f", "g"]))
        exp = "(((e,f),g),((c,d),(b,a)));\n"
        self.assertEqual(obs, exp)

        # outgroup is a single taxon
        obs = str(tree.root_by_outgroup(["a"]))
        exp = "(a,(b,((c,d),((e,f),g))));\n"
        self.assertEqual(obs, exp)

        # outgroup is not monophyletic
        msg = "Outgroup is not monophyletic in the tree."
        with self.assertRaisesRegex(TreeError, msg):
            tree.root_by_outgroup(["a", "c"])

        # outgroup has extra taxa
        msg = "Outgroup is not a proper subset of taxa in the tree."
        with self.assertRaisesRegex(TreeError, msg):
            tree.root_by_outgroup(["a", "b", "x"])

        # outgroup is not in tree
        with self.assertRaisesRegex(TreeError, msg):
            tree.root_by_outgroup(["x", "y"])

        # outgroup is the whole tree
        with self.assertRaisesRegex(TreeError, msg):
            tree.root_by_outgroup("abcdefg")

        # generate unrooted tree
        obs = str(tree.root_by_outgroup(["a", "b"], above=False))
        exp = "(a,b,((c,d),((e,f),g)));\n"
        self.assertEqual(obs, exp)

        # keep old root node
        obs = str(tree.root_by_outgroup(["a", "b"], reset=False))
        exp = "((a,b),((c,d),((e,f),(g))));\n"
        self.assertEqual(obs, exp)

        # specify root name
        obs = str(tree.root_by_outgroup(["a", "b"], root_name="root"))
        exp = "((a,b),((c,d),((e,f),g)))root;\n"
        self.assertEqual(obs, exp)

        # in-place rooting
        lca = tree.lca(["a", "b"])
        self.assertIsNot(obs, lca.parent)
        obs = tree.root_by_outgroup(["a", "b"], root_name="root", inplace=True)
        self.assertEqual(str(obs), exp)
        self.assertIs(obs, lca.parent)

        # transfer branch support
        tree = TreeNode.read(["((((a,b)80,(c,d)),(e,f)),g);"])
        tree.assign_supports()
        obs = str(tree.root_by_outgroup(["a", "b"]))
        exp = "((a,b)80,((c,d),((e,f),g))80);\n"
        self.assertEqual(obs, exp)

        # transfer custom branch attribute
        tree = TreeNode.read(["((((a,b),(c,d))x,(e,f)),g);"])
        obs = str(tree.root_by_outgroup(["a", "b"], branch_attrs=["name"]))
        exp = "((a,b),((c,d),((e,f),g)x));\n"
        self.assertEqual(obs, exp)

    # ------------------------------------------------
    # Tree metrics
    # ------------------------------------------------

    def test_count(self):
        """Get node counts"""
        exp = 7
        obs = self.simple_t.count()
        self.assertEqual(obs, exp)

        exp = 4
        obs = self.simple_t.count(tips=True)
        self.assertEqual(obs, exp)

    def test_subset(self):
        """subset should return set of leaves that descends from node"""
        t = self.simple_t
        self.assertEqual(t.subset(), frozenset("abcd"))
        c = t.children[0]
        self.assertEqual(c.subset(), frozenset("ab"))
        leaf = c.children[1]
        self.assertEqual(leaf.subset(), frozenset(""))

    def test_subsets(self):
        """subsets should return all subsets descending from a set"""
        t = self.simple_t
        self.assertEqual(t.subsets(), frozenset(
            [frozenset("ab"), frozenset("cd")]))

    def test_assign_supports(self):
        """Extract support values of internal nodes."""
        # test nodes with support values alone as labels
        tree = TreeNode.read(["((a,b)75,(c,d)90);"])
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
        for taxon in ("a", "b", "c", "d"):
            self.assertIsNone(tree.find(taxon).support)

        # test nodes with support values and branch lengths
        tree = TreeNode.read(["((a,b)0.85:1.23,(c,d)0.95:4.56);"])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, 0.85)
        self.assertEqual(node2.support, 0.95)

        # test whether integer or float support values can be correctly parsed
        tree = TreeNode.read(["((a,b)75,(c,d)80.0,(e,f)97.5,(g,h)0.95);"])
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
        tree = TreeNode.read(["((a,b)-1.23,(c,d)1.23e-4);"])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, -1.23)
        self.assertEqual(node2.support, 0.000123)

        # test nodes with support and extra label
        tree = TreeNode.read(["((a,b)'80:X',(c,d)'60:Y');"])
        tree.assign_supports()
        node1, node2 = tree.children
        self.assertEqual(node1.support, 80)
        self.assertEqual(node1.name, "X")
        self.assertEqual(node2.support, 60)
        self.assertEqual(node2.name, "Y")

        # test nodes without label, with non-numeric label, and with branch
        # length only
        tree = TreeNode.read(["((a,b),(c,d)x,(e,f):1.0);"])
        tree.assign_supports()
        for node in tree.children:
            self.assertIsNone(node.support)

    def test_is_bifurcating(self):
        """Check if tree is bifurcating."""
        t = self.simple_t
        self.assertTrue(t.is_bifurcating())
        t = TreeNode.read(["((a,b,c),(d,e))root;"])
        self.assertFalse(t.is_bifurcating())
        t = TreeNode.read(["((((a,b)c)d)e,f)root;"])
        self.assertTrue(t.is_bifurcating())
        self.assertFalse(t.is_bifurcating(strict=True))

    def test_observed_node_counts(self):
        """returns observed nodes counts given vector of observed taxon counts
        """
        t = self.simple_t

        # no taxon observed
        taxon_counts = {}
        exp = defaultdict(int)
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        # error on zero count(s)
        taxon_counts = {"a": 0}
        self.assertRaises(ValueError, t.observed_node_counts, taxon_counts)
        taxon_counts = {"a": 0, "b": 0, "c": 0, "d": 0}
        self.assertRaises(ValueError, t.observed_node_counts, taxon_counts)

        # all taxa observed once
        taxon_counts = {"a": 1, "b": 1, "c": 1, "d": 1}
        exp = defaultdict(int)
        exp[t.find("root")] = 4
        exp[t.find("i1")] = 2
        exp[t.find("i2")] = 2
        exp[t.find("a")] = 1
        exp[t.find("b")] = 1
        exp[t.find("c")] = 1
        exp[t.find("d")] = 1
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        # some taxa observed twice
        taxon_counts = {"a": 2, "b": 1, "c": 1, "d": 1}
        exp = defaultdict(int)
        exp[t.find("root")] = 5
        exp[t.find("i1")] = 3
        exp[t.find("i2")] = 2
        exp[t.find("a")] = 2
        exp[t.find("b")] = 1
        exp[t.find("c")] = 1
        exp[t.find("d")] = 1
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        taxon_counts = {"a": 2, "b": 1, "c": 1, "d": 2}
        exp = defaultdict(int)
        exp[t.find("root")] = 6
        exp[t.find("i1")] = 3
        exp[t.find("i2")] = 3
        exp[t.find("a")] = 2
        exp[t.find("b")] = 1
        exp[t.find("c")] = 1
        exp[t.find("d")] = 2
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        # some taxa observed, others not observed
        taxon_counts = {"a": 2, "b": 1}
        exp = defaultdict(int)
        exp[t.find("root")] = 3
        exp[t.find("i1")] = 3
        exp[t.find("a")] = 2
        exp[t.find("b")] = 1
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        taxon_counts = {"d": 1}
        exp = defaultdict(int)
        exp[t.find("root")] = 1
        exp[t.find("i2")] = 1
        exp[t.find("d")] = 1
        self.assertEqual(t.observed_node_counts(taxon_counts), exp)

        # error on non-tips
        taxon_counts = {"a": 2, "e": 1}
        self.assertRaises(MissingNodeError, t.observed_node_counts, taxon_counts)
        taxon_counts = {"a": 2, "i1": 1}
        self.assertRaises(MissingNodeError, t.observed_node_counts, taxon_counts)

        # test with another tree
        t2 = self.complex_tree
        taxon_counts = {}
        exp = defaultdict(int)
        self.assertEqual(t2.observed_node_counts(taxon_counts), exp)

        taxon_counts = {"e": 42, "f": 1}
        exp[t2.root()] = 43
        exp[t2.find("int5")] = 43
        exp[t2.find("e")] = 42
        exp[t2.find("f")] = 1
        self.assertEqual(t2.observed_node_counts(taxon_counts), exp)

    def test_accumulate_to_ancestor(self):
        """Get the distance from a node to its ancestor"""
        t = TreeNode.read([
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;"])
        a = t.find("a")
        b = t.find("b")
        exp_to_root = 0.1 + 0.3
        obs_to_root = a.accumulate_to_ancestor(t)
        self.assertEqual(obs_to_root, exp_to_root)

        with self.assertRaises(NoParentError):
            a.accumulate_to_ancestor(b)

    def test_descending_branch_length(self):
        """Calculate descending branch_length"""
        tr = TreeNode.read([
            "(((A:.1,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tdbl = tr.descending_branch_length()
        sdbl = tr.descending_branch_length(["A", "E"])
        npt.assert_almost_equal(tdbl, 8.9)
        npt.assert_almost_equal(sdbl, 2.2)
        self.assertRaises(ValueError, tr.descending_branch_length,
                          ["A", "DNE"])
        self.assertRaises(ValueError, tr.descending_branch_length, ["A", "C"])

        tr = TreeNode.read([
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 8.8)

        tr = TreeNode.read([
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 7.9)

        tr = TreeNode.read([
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tdbl = tr.descending_branch_length(["A", "D", "E"])
        npt.assert_almost_equal(tdbl, 2.1)

        tr = TreeNode.read([
            "(((A,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tdbl = tr.descending_branch_length(["I", "D", "E"])
        npt.assert_almost_equal(tdbl, 6.6)

        # test with a situation where we have unnamed internal nodes
        tr = TreeNode.read([
            "(((A,B:1.2):.6,(D:.9,E:.6)F):2.4,(H:.4,I:.5)J:1.3);"])
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 7.9)

        # issue 1847
        tr = TreeNode.read([
            "(((A:.1,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        tr.length = 1
        tdbl = tr.descending_branch_length()
        npt.assert_almost_equal(tdbl, 8.9)

    def test_distance_nontip(self):
        # example derived from issue #807, credit @wwood
        tstr = "((A:1.0,B:2.0)'g__genus1':3.0)root;"
        tree = TreeNode.read([tstr])
        self.assertEqual(tree.find("A").distance(tree.find("g__genus1")), 1.0)

    def test_distance(self):
        """Get the distance between two nodes"""
        t = TreeNode.read(["((a:0.1,b:0.2)c:0.3,(d:0.4,e)f:0.5)root;"])
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

    def test_get_max_distance(self):
        """get_max_distance should get max tip distance across tree"""
        tree = TreeNode.read([
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;"])
        dist, nodes = tree.get_max_distance()
        npt.assert_almost_equal(dist, 1.6)
        self.assertEqual(sorted([n.name for n in nodes]), ["b", "e"])

    def test_set_max_distance(self):
        """set_max_distance sets MaxDistTips across tree"""
        tree = TreeNode.read([
            "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.5)f:0.6)root;"])
        tree._set_max_distance()
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 1.6)
        self.assertEqual(sorted([tip_a[1].name, tip_b[1].name]), ["b", "e"])

    def test_set_max_distance_tie_bug(self):
        """Corresponds to #1077"""
        t = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5)root;"])
        exp = ((3.0, t.find("a")), (9.0, t.find("e")))

        # the above tree would trigger an exception in max. The central issue
        # was that the data being passed to max were a tuple of tuple:
        # ((left_d, left_n), (right_d, right_n))
        # the call to max would break in this scenario as it would fall onto
        # idx 1 of each tuple to assess the "max".
        t._set_max_distance()

        self.assertEqual(t.MaxDistTips, exp)

    def test_set_max_distance_inplace_modification_bug(self):
        """Corresponds to #1223"""
        t = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5)root;"])

        exp = [((0.0, t.find("a")), (0.0, t.find("a"))),
               ((0.0, t.find("b")), (0.0, t.find("b"))),
               ((1.0, t.find("a")), (1.0, t.find("b"))),
               ((0.0, t.find("d")), (0.0, t.find("d"))),
               ((0.0, t.find("e")), (0.0, t.find("e"))),
               ((3.0, t.find("d")), (4.0, t.find("e"))),
               ((3.0, t.find("a")), (9.0, t.find("e")))]

        t._set_max_distance()

        self.assertEqual([n.MaxDistTips for n in t.postorder()], exp)

    def test_tip_tip_distances_endpoints(self):
        """Test getting specifc tip distances  with tipToTipDistances"""
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        nodes = [t.find("H"), t.find("G"), t.find("M")]
        names = ["H", "G", "M"]
        exp = DistanceMatrix(np.array([[0, 2.0, 6.7],
                                       [2.0, 0, 6.7],
                                       [6.7, 6.7, 0.0]]), ["H", "G", "M"])

        obs = t.tip_tip_distances(endpoints=names)
        self.assertEqual(obs, exp)

        obs = t.tip_tip_distances(endpoints=nodes)
        self.assertEqual(obs, exp)

    def test_tip_tip_distances_non_tip_endpoints(self):
        t = TreeNode.read(["((H:1,G:1)foo:2,(R:0.5,M:0.7):3);"])
        with self.assertRaises(ValueError):
            t.tip_tip_distances(endpoints=["foo"])

    def test_tip_tip_distances_no_length(self):
        t = TreeNode.read(["((a,b)c,(d,e)f);"])
        exp_t = TreeNode.read(["((a:0,b:0)c:0,(d:0,e:0)f:0);"])
        exp_t_dm = exp_t.tip_tip_distances()

        t_dm = npt.assert_warns(RepresentationWarning, t.tip_tip_distances)
        self.assertEqual(t_dm, exp_t_dm)

        for node in t.preorder():
            self.assertIs(node.length, None)

    def test_tip_tip_distances_missing_length(self):
        t = TreeNode.read(["((a,b:6)c:4,(d,e:0)f);"])
        exp_t = TreeNode.read(["((a:0,b:6)c:4,(d:0,e:0)f:0);"])
        exp_t_dm = exp_t.tip_tip_distances()

        t_dm = npt.assert_warns(RepresentationWarning, t.tip_tip_distances)
        self.assertEqual(t_dm, exp_t_dm)

    def test_compare_rfd(self):
        """compare_rfd should return the Robinson Foulds distance"""
        t = TreeNode.read(["((H,G),(R,M));"])
        t2 = TreeNode.read(["(((H,G),R),M);"])
        t4 = TreeNode.read(["(((H,G),(O,R)),X);"])

        obs = t.compare_rfd(t2)
        exp = 2.0
        self.assertEqual(obs, exp)

        self.assertEqual(t.compare_rfd(t2), t2.compare_rfd(t))

        obs = t.compare_rfd(t2, proportion=True)
        exp = 0.5
        self.assertEqual(obs, exp)

        with self.assertRaises(ValueError):
            t.compare_rfd(t4)

    def test_compare_subsets(self):
        """compare_subsets should return the fraction of shared subsets"""
        t = TreeNode.read(["((H,G),(R,M));"])
        t2 = TreeNode.read(["(((H,G),R),M);"])
        t4 = TreeNode.read(["(((H,G),(O,R)),X);"])

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

    def test_compare_tip_distances(self):
        # default behavior
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        t2 = TreeNode.read(["(((H:1,G:1,O:1):2,R:3):1,X:4);"])
        obs = t.compare_tip_distances(t2)
        # note: common taxa are H, G, R (only)
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

        # sample a subset of taxa
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        t2 = TreeNode.read(["(((H:1,G:1,O:1):2,R:3):1,X:4);"])
        obs = t.compare_tip_distances(t2, sample=3)
        # Note: common taxa are H, G, R (only), all of which are selected, despite that
        # the default shuffling function is stochastic.
        m1 = np.array([[0, 2, 6.5], [2, 0, 6.5], [6.5, 6.5, 0]])
        m2 = np.array([[0, 2, 6], [2, 0, 6], [6, 6, 0]])
        r = pearsonr(m1.flat, m2.flat)[0]
        self.assertAlmostEqual(obs, (1 - r) / 2)

        # 4 common taxa, custom shuffling function, still picking H, G, R
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7,Q:5):3);"])
        t3 = TreeNode.read(["(((H:1,G:1,O:1):2,R:3,Q:10):1,X:4);"])
        obs = t.compare_tip_distances(t3, sample=3, shuffle_f=sorted)

        # no common taxa
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        t2 = TreeNode.read(["(((Z:1,Y:1,X:1):2,W:3):1,V:4);"])
        with self.assertRaises(ValueError):
            t.compare_tip_distances(t2)

        # single common taxon
        t = TreeNode.read(["((H:1,G:1):2,(R:0.5,M:0.7):3);"])
        t2 = TreeNode.read(["(((R:1,Y:1,X:1):2,W:3):1,V:4);"])
        self.assertEqual(t.compare_tip_distances(t2), 1)
        self.assertEqual(t2.compare_tip_distances(t), 1)

    # ------------------------------------------------
    # Tree indexing and searching
    # ------------------------------------------------

    def test_has_caches(self):
        """Check if tree has caches."""
        t = self.simple_t
        self.assertTupleEqual(t.has_caches(), (None, False))
        t.find("a")
        self.assertTupleEqual(t.has_caches(), (None, True))
        t.cache_attr(lambda n: 1, "node_count", sum)
        self.assertTupleEqual(t.has_caches(), ({"node_count"}, True))
        t.clear_caches()
        self.assertTupleEqual(t.has_caches(), (None, False))

    def test_clear_caches(self):
        # delete lookup caches
        t = TreeNode.read(["((a:1.2,b:1.6)c:0.3,(d:0.8,e:1.0)f:0.6)g;"])
        t.create_caches()
        keys = ("_tip_cache", "_non_tip_cache")
        for key in keys:
            self.assertTrue(hasattr(t, key))
        t.clear_caches(attr=False)
        for key in keys:
            self.assertFalse(hasattr(t, key))

        # delete all attribute caches
        t.cache_attr(lambda n: [n.name] if n.is_tip() else [], "tip_names")
        delattr(t.children[0], "tip_names")
        t.clear_caches(lookup=False)
        self.assertFalse(hasattr(t, "_registered_caches"))
        for node in t.traverse(include_self=True):
            self.assertFalse(hasattr(node, "tip_names"))

        # delete individual attribute caches
        t.cache_attr(lambda n: 1, "node_count", sum)
        t.cache_attr(lambda n: n.length or 0.0, "total_length", sum)
        t.clear_caches(attr="node_count")
        self.assertTrue(hasattr(t, "_registered_caches"))
        self.assertNotIn("node_count", t._registered_caches)
        self.assertIn("total_length", t._registered_caches)
        for node in t.traverse(include_self=True):
            self.assertFalse(hasattr(node, "node_count"))
            self.assertTrue(hasattr(node, "total_length"))
        delattr(t.children[1], "total_length")
        t.clear_caches(attr="total_length")
        self.assertFalse(hasattr(t, "_registered_caches"))
        for node in t.traverse(include_self=True):
            self.assertFalse(hasattr(node, "total_length"))

    def test_cache_attr(self):
        # cache names of all descending tips
        t = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        f = lambda n: [n.name] if n.is_tip() else []
        t.cache_attr(f, "tip_names")
        self.assertIn("tip_names", t._registered_caches)
        self.assertListEqual(t.tip_names, list("abde"))
        self.assertListEqual(t.children[0].tip_names, list("ab"))
        self.assertListEqual(t.children[1].tip_names, list("de"))

        # don't register as cache
        t.clear_caches()
        t.cache_attr(f, "tip_names", register=False)
        self.assertFalse(hasattr(t, "_registered_caches"))
        self.assertListEqual(t.tip_names, list("abde"))

        # tuple instead of list
        t.cache_attr(f, "tip_names", tuple)
        self.assertTupleEqual(t.tip_names, tuple("abde"))

        # set and frozenset
        t.cache_attr(f, "tip_names", set)
        self.assertIs(type(t.tip_names), set)
        self.assertSetEqual(t.tip_names, set("abde"))
        t.cache_attr(f, "tip_names", frozenset)
        self.assertIs(type(t.tip_names), frozenset)
        self.assertSetEqual(t.tip_names, set("abde"))

        # cache number of nodes per clade
        t = TreeNode.read(["((a:1.2,b:1.6)c:0.3,(d:0.8,e:1.0)f:0.6)g;"])
        f = lambda n: 1
        t.cache_attr(f, "node_count", sum)
        self.assertEqual(t.node_count, 7)
        self.assertEqual(t.children[0].node_count, 3)
        self.assertEqual(t.children[1].node_count, 3)

        # cache total branch length per clade
        t.clear_caches()
        f = lambda n: n.length or 0.0
        t.cache_attr(f, "total_length", sum)
        self.assertAlmostEqual(t.total_length, 5.5)
        self.assertAlmostEqual(t.children[0].total_length, 3.1)
        self.assertAlmostEqual(t.children[1].total_length, 2.4)

        # cache accumulative distance from tips using a custom function
        t.clear_caches()
        dist_f = lambda x: np.array(x.length or 0.0, ndmin=1)
        comb_f = lambda prev, curr: np.concatenate(prev) + curr if prev else curr
        t.cache_attr(dist_f, "accu_dist", comb_f)
        npt.assert_almost_equal(t.accu_dist, np.array([1.5, 1.9, 1.4, 1.6]))
        npt.assert_almost_equal(t.children[0].accu_dist, np.array([1.5, 1.9]))

        # cache and combine accumulative distance using a custom function
        t.clear_caches()

        def depths_f(node):
            if node.is_tip():
                return [0.0]
            else:
                return [y + (x.length or 0.0) for x in node.children for y in x.depths]

        t.cache_attr(depths_f, "depths", None)
        for obs, exp in zip(t.depths, [1.5, 1.9, 1.4, 1.6]):
            self.assertAlmostEqual(obs, exp)
        for obs, exp in zip(t.children[0].depths, [1.2, 1.6]):
            self.assertAlmostEqual(obs, exp)

        # invalid cache type
        msg = "Cache type is invalid."
        with self.assertRaisesRegex(TypeError, msg):
            t.cache_attr(sum, "missing", "invalid")

    def test_assign_ids(self):
        """Assign IDs to the tree"""
        t1 = TreeNode.read(["(((a,b),c),(e,f),(g));"])
        t2 = TreeNode.read(["(((a,b),c),(e,f),(g));"])
        t3 = TreeNode.read(["((g),(e,f),(c,(a,b)));"])
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
        t1 = TreeNode.read(["(((a,b),c),(d,e));"])
        t2 = TreeNode.read(["(((a,b),(c,d)),(e,f));"])
        t3 = TreeNode.read(["(((a,b,c),(d)),(e,f));"])
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

    def test_index_tree(self):
        """index_tree should produce correct index and node map"""
        # test for first tree: contains singleton outgroup
        t1 = TreeNode.read(["(((a,b),c),(d,e));"])
        t2 = TreeNode.read(["(((a,b),(c,d)),(e,f));"])
        t3 = TreeNode.read(["(((a,b,c),(d)),(e,f));"])

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

        # test for single-node tree
        t1 = TreeNode.read(["root;"])
        id_index, child_index = t1.index_tree()
        self.assertEqual(id_index[0], t1)
        npt.assert_equal(child_index, np.array([[]]))

    def test_create_caches(self):
        t = TreeNode.read(["(((a,b)x,(c,d)x,e),(f,g)y)root;"])

        # create a lookup table for a fresh tree
        t.create_caches()
        self.assertEqual(t._tip_cache["a"].name, "a")
        self.assertEqual(len(t._non_tip_cache["x"]), 2)
        self.assertEqual(len(t._non_tip_cache["y"][0].children), 2)

        # skip re-creating a lookup table if the tree already has it
        t._tip_cache["a"] = None
        t.create_caches()
        self.assertIsNone(t._tip_cache["a"])

        # can create a lookup table for the entire tree from any node
        node = t.find("b")
        t.clear_caches(attr=False)
        node.create_caches()
        self.assertEqual(t._tip_cache["c"].name, "c")
        self.assertListEqual([
            x.name for x in t._non_tip_cache["y"][0].children], ["f", "g"])

        # raise if duplicate tip names found
        msg = "Duplicate tip name 'a' found."
        with self.assertRaisesRegex(DuplicateNodeError, msg):
            TreeNode.read(["(a,a);"]).create_caches()

    def test_find(self):
        t = TreeNode.read(["((a,b)c,(d,e)f);"])

        # find an internal node
        exp = t.children[0]
        obs = t.find("c")
        self.assertEqual(obs, exp)

        # find a tip
        exp = t.children[0].children[1]
        obs = t.find("b")
        self.assertEqual(obs, exp)

        # input is node
        obs = t.find(exp)
        self.assertIs(obs, exp)

        # name not found
        msg = "Node 'missing' is not found."
        with self.assertRaisesRegex(MissingNodeError, msg):
            t.find("missing")

    def test_find_all(self):
        t = TreeNode.read(["((a,b)c,((d,e)c)c,(f,(g,h)c)a)root;"])

        # find all nodes with a given name
        exp = [t.children[2],
               t.children[0].children[0]]
        obs = t.find_all("a")
        self.assertEqual(obs, exp)

        exp = [t.children[0],
               t.children[1].children[0],
               t.children[1],
               t.children[2].children[1]]
        obs = t.find_all("c")
        self.assertEqual(obs, exp)

        # input is TreeNode
        obs = t.find_all(t.children[0])
        self.assertEqual(obs, exp)

        # find root itself
        obs = t.find_all("root")
        self.assertEqual(len(obs), 1)
        self.assertIs(obs[0], t)

        # node not found
        msg = "Node 'missing' is not found."
        with self.assertRaisesRegex(MissingNodeError, msg):
            t.find_all("missing")

    def test_find_by_id(self):
        """Find a node by id"""
        t1 = TreeNode.read(["((,),(,,));"])
        t2 = TreeNode.read(["((,),(,,));"])

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
        t = TreeNode.read(["((a,b)c,(d,e)f);"])

        def func(x):
            return x.parent == t.find("c")

        exp = ["a", "b"]
        obs = [n.name for n in t.find_by_func(func)]
        self.assertEqual(obs, exp)

    # ------------------------------------------------
    # Tree visualization
    # ------------------------------------------------

    def test_ascii_art(self):
        """Make some ascii trees"""
        # unlabeled internal node
        tr = TreeNode.read(["(B:0.2,(C:0.3,D:0.4):0.6)F;"])
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
        tr = TreeNode.read(["(B:0.2,(C:0.3,D:0.4)90:0.6)F;"])
        exp = "          /-B\n-F-------|\n         |          /-C\n         "\
              " \\90------|\n                    \\-D"
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr.assign_supports()
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr = TreeNode.read(["((A,B)75,(C,D)'80:spA');"])
        exp = "                    /-A\n          /75------|\n         |    "\
              "      \\-B\n---------|\n         |          /-C\n          \\"\
              "80:spA--|\n                    \\-D"
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)
        tr.assign_supports()
        obs = tr.ascii_art(show_internal=True, compact=False)
        self.assertEqual(obs, exp)

    def test_ascii_art_three_children(self):
        obs = TreeNode.read(["(a,(b,c,d));"]).ascii_art()
        self.assertEqual(obs, exp_ascii_art_three_children)

    # ------------------------------------------------
    # Format conversion
    # ------------------------------------------------

    def test_from_linkage_matrix(self):
        # Ensure matches: http://www.southampton.ac.uk/~re1u06/teaching/upgma/
        id_list = ["A", "B", "C", "D", "E", "F", "G"]
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

    def test_from_taxonomy(self):
        lineages = [("1", ["a", "b", "c", "d", "e", "f", "g"]),
                    ("2", ["a", "b", "c", None, None, "x", "y"]),
                    ("3", ["h", "i", "j", "k", "l", "m", "n"]),
                    ("4", ["h", "i", "j", "k", "l", "m", "q"]),
                    ("5", ["h", "i", "j", "k", "l", "m", "n"])]
        exp = TreeNode.read([
            "((((((((1)g)f)e)d,((((2)y)x)))c)b)a,"
            "(((((((3,5)n,(4)q)m)l)k)j)i)h);"])

        # input as 2-element tuples
        obs = TreeNode.from_taxonomy(lineages)
        self.assertIs(type(obs), TreeNode)
        self.assertEqual(obs.compare_subsets(exp), 0.0)

        obs = TreeNodeSubclass.from_taxonomy(lineages)
        self.assertIs(type(obs), TreeNodeSubclass)

        # input as dictionary
        dict_ = dict(lineages)
        obs = TreeNode.from_taxonomy(dict_)
        self.assertEqual(obs.compare_subsets(exp), 0.0)

        # input as data frame
        df_ = pd.DataFrame([x[1] for x in lineages], [x[0] for x in lineages])
        obs = TreeNode.from_taxonomy(df_)
        self.assertEqual(obs.compare_subsets(exp), 0.0)

    def test_to_taxonomy(self):
        input_lineages = {"1": ["a", "b", "c", "d", "e", "f", "g"],
                          "2": ["a", "b", "c", None, None, "x", "y"],
                          "3": ["h", "i", "j", "k", "l", "m", "n"],
                          "4": ["h", "i", "j", "k", "l", "m", "q"],
                          "5": ["h", "i", "j", "k", "l", "m", "n"]}
        tree = TreeNode.from_taxonomy(input_lineages.items())
        exp = sorted(input_lineages.items())
        obs = [(n.name, lin) for n, lin in tree.to_taxonomy(allow_empty=True)]
        self.assertEqual(sorted(obs), exp)

    def test_to_taxonomy_filter(self):
        input_lineages = {"1": ["a", "b", "c", "d", "e", "f", "g"],
                          "2": ["a", "b", "c", None, None, "x", "y"],
                          "3": ["h", "i", "j", "k", "l"],  # test jagged
                          "4": ["h", "i", "j", "k", "l", "m", "q"],
                          "5": ["h", "i", "j", "k", "l", "m", "n"]}
        tree = TreeNode.from_taxonomy(input_lineages.items())

        def f(node, lin):
            return "k" in lin or "x" in lin

        exp = [("2", ["a", "b", "c", "x", "y"]),
               ("3", ["h", "i", "j", "k", "l"]),
               ("4", ["h", "i", "j", "k", "l", "m", "q"]),
               ("5", ["h", "i", "j", "k", "l", "m", "n"])]
        obs = [(n.name, lin) for n, lin in tree.to_taxonomy(filter_f=f)]
        self.assertEqual(sorted(obs), exp)

    def test_from_taxdump(self):
        # same example as in skbio.io.format.taxdump
        nodes = pd.DataFrame([
            [1,       1,      "no rank"],
            [2,       131567, "superkingdom"],
            [543,     91347,  "family"],
            [548,     570,    "species"],
            [561,     543,    "genus"],
            [562,     561,    "species"],
            [570,     543,    "genus"],
            [620,     543,    "genus"],
            [622,     620,    "species"],
            [766,     28211,  "order"],
            [1224,    2,      "phylum"],
            [1236,    1224,   "class"],
            [28211,   1224,   "class"],
            [91347,   1236,   "order"],
            [118884,  1236,   "no rank"],
            [126792,  36549,  "species"],
            [131567,  1,      "no rank"],
            [585056,  562,    "no rank"],
            [1038927, 562,    "no rank"],
            [2580236, 488338, "species"]],
            columns=["tax_id", "parent_tax_id", "rank"]).set_index("tax_id")

        names = pd.DataFrame([
            [1, "root", np.nan, "scientific name"],
            [2, "Bacteria", "Bacteria <bacteria>", "scientific name"],
            [2, "eubacteria", np.nan, "genbank common name"],
            [543, "Enterobacteriaceae", np.nan, "scientific name"],
            [548, "Klebsiella aerogenes", np.nan, "scientific name"],
            [561, "Escherichia", np.nan, "scientific name"],
            [562, "\"Bacillus coli\" Migula 1895", np.nan, "authority"],
            [562, "Escherichia coli", np.nan, "scientific name"],
            [562, "Escherichia/Shigella coli", np.nan, "equivalent name"],
            [570, "Donovania", np.nan, "synonym"],
            [570, "Klebsiella", np.nan, "scientific name"],
            [620, "Shigella", np.nan, "scientific name"],
            [622, "Shigella dysenteriae", np.nan, "scientific name"],
            [766, "Rickettsiales", np.nan, "scientific name"],
            [1224, "Proteobacteria", np.nan, "scientific name"],
            [1236, "Gammaproteobacteria", np.nan, "scientific name"],
            [28211, "Alphaproteobacteria", np.nan, "scientific name"],
            [91347, "Enterobacterales", np.nan, "scientific name"],
            [118884, "unclassified Gammaproteobacteria", np.nan,
             "scientific name"],
            [126792, "Plasmid pPY113", np.nan, "scientific name"],
            [131567, "cellular organisms", np.nan, "scientific name"],
            [585056, "Escherichia coli UMN026", np.nan, "scientific name"],
            [1038927, "Escherichia coli O104:H4", np.nan, "scientific name"],
            [2580236, "synthetic Escherichia coli Syn61", np.nan,
             "scientific name"]],
             columns=["tax_id", "name_txt", "unique_name",
                      "name_class"]).set_index("tax_id")

        # nodes without names (use tax_id as name)
        obs = TreeNode.from_taxdump(nodes)
        exp = ("(((((((((585056,1038927)562)561,(548)570,(622)620)543)91347,"
               "118884)1236,(766)28211)1224)2)131567)1;")
        self.assertEqual(str(obs).rstrip(), exp)
        self.assertEqual(obs.count(), 18)
        self.assertEqual(obs.count(tips=True), 6)

        # default scenario (nodes and names)
        obs = TreeNode.from_taxdump(nodes, names)

        # check tree is in same size
        self.assertEqual(obs.count(), 18)
        self.assertEqual(obs.count(tips=True), 6)

        # check id, name and rank are correctly set at root
        self.assertEqual(obs.id, 1)
        self.assertEqual(obs.name, "root")
        self.assertEqual(obs.rank, "no rank")

        # check an internal node
        node = obs.find("Enterobacteriaceae")
        self.assertEqual(node.id, 543)
        self.assertEqual(node.rank, "family")

        # check its children (which should preserve input order)
        self.assertEqual(len(node.children), 3)
        self.assertListEqual([x.name for x in node.children], [
            "Escherichia", "Klebsiella", "Shigella"])

        # check that non-scientific name isn"t used
        with self.assertRaises(MissingNodeError):
            obs.find("Donovania")

        # name as a dictionary
        names = names[names["name_class"] == "scientific name"][
            "name_txt"].to_dict()
        obs = TreeNode.from_taxdump(nodes, names)
        self.assertEqual(obs.count(), 18)
        self.assertEqual(obs.name, "root")
        self.assertEqual(obs.find("Enterobacteriaceae").id, 543)

        # nodes has no top level
        nodes = pd.DataFrame([
            [2, 1, "A"],
            [3, 2, "B"],
            [1, 3, "C"]],
            columns=["tax_id", "parent_tax_id", "rank"]).set_index("tax_id")
        with self.assertRaises(ValueError) as ctx:
            TreeNode.from_taxdump(nodes)
        self.assertEqual(str(ctx.exception), "There is no top-level node.")

        # nodes has more than one top level
        nodes = pd.DataFrame([
            [1, 1, "A"],
            [2, 2, "B"],
            [3, 3, "C"]],
            columns=["tax_id", "parent_tax_id", "rank"]).set_index("tax_id")
        with self.assertRaises(ValueError) as ctx:
            TreeNode.from_taxdump(nodes)
        self.assertEqual(str(
            ctx.exception), "There are more than one top-level node.")

    def test_to_array(self):
        """Convert a tree to arrays"""
        t = TreeNode.read([
            "(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10);"])
        id_index, child_index = t.index_tree()
        arrayed = t.to_array()

        self.assertEqual(id_index, arrayed["id_index"])
        npt.assert_equal(child_index, arrayed["child_index"])

        exp = np.array([1, 2, 3, 5, 4, 6, 8, 9, 7, 10, np.nan])
        obs = arrayed["length"]
        npt.assert_equal(obs, exp)

        exp = np.array(["a", "b", "c", "d", "x",
                        "y", "e", "f", "z", "z", None])
        obs = arrayed["name"]
        npt.assert_equal(obs, exp)

        exp = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        obs = arrayed["id"]
        npt.assert_equal(obs, exp)

    def test_to_array_attrs(self):
        t = TreeNode.read([
            "(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10);"])
        id_index, child_index = t.index_tree()
        arrayed = t.to_array(attrs=[("name", object)])

        # should only have id_index, child_index, and name since we specified
        # attrs
        self.assertEqual(len(arrayed), 3)

        self.assertEqual(id_index, arrayed["id_index"])
        npt.assert_equal(child_index, arrayed["child_index"])

        exp = np.array(["a", "b", "c", "d", "x",
                        "y", "e", "f", "z", "z", None])
        obs = arrayed["name"]
        npt.assert_equal(obs, exp)

        # invalid attrs
        with self.assertRaises(AttributeError):
            t.to_array(attrs=[("name", object), ("brofist", int)])

    def test_to_array_nan_length_value(self):
        t = TreeNode.read(["((a:1, b:2)c:3)root;"])
        indexed = t.to_array(nan_length_value=None)
        npt.assert_equal(indexed["length"],
                         np.array([1, 2, 3, np.nan], dtype=float))
        indexed = t.to_array(nan_length_value=0.0)
        npt.assert_equal(indexed["length"],
                         np.array([1, 2, 3, 0.0], dtype=float))
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed["length"],
                         np.array([1, 2, 3, 42.0], dtype=float))

        t = TreeNode.read(["((a:1, b:2)c:3)root:4;"])
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed["length"],
                         np.array([1, 2, 3, 4], dtype=float))

        t = TreeNode.read(["((a:1, b:2)c)root;"])
        indexed = t.to_array(nan_length_value=42.0)
        npt.assert_equal(indexed["length"],
                         np.array([1, 2, 42.0, 42.0], dtype=float))


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

empty = "();"
single = "(abc:3);"
double = "(abc:3, def:4);"
onenest = "(abc:3, (def:4, ghi:5):6 );"
nodedata = "(abc:3, (def:4, ghi:5)jkl:6 );"

exp_ascii_art_three_children = r"""          /-a
         |
---------|          /-b
         |         |
          \--------|--c
                   |
                    \-d"""


if __name__ == "__main__":
    main()
