# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import unittest

import numpy as np
import numpy.testing as npt

from skbio import TreeNode
from skbio.tree import BPTree
from skbio.io import NewickFormatError
from skbio.io.format.newick import (
    _newick_to_tree_node, _tree_node_to_newick, _newick_sniffer,
    _newick_to_bp, _bp_to_newick)


class TestNewick(unittest.TestCase):
    def _assert_node_equal(self, n1, n2):
        self.assertEqual(n1.name, n2.name)
        self.assertEqual(n1.length, n2.length)
        self.assertEqual(len(n1.children), len(n2.children))

    def _assert_equal(self, n1, n2):
        def name(x):
            return (str(x.name),
                    float(x.length) if x.length is not None else 0,
                    len(x.children))
        self._assert_node_equal(n1, n2)
        for c1, c2 in zip(sorted(n1.children, key=name),
                          sorted(n2.children, key=name)):
            self.assertTrue(c1.parent is n1)
            self.assertTrue(c2.parent is n2)
            self._assert_equal(c1, c2)

    def _setup_tree(self, kwargs_list):
        trees = []
        for kwargs in kwargs_list:
            trees.append(TreeNode(**kwargs))

        trees[4].extend([trees[2], trees[3]])
        trees[5].extend([trees[0], trees[1], trees[4]])

        return trees[5]

    def _setup_linked_list(self, kwargs_list):
        last_node = None
        for idx, kwargs in enumerate(kwargs_list):
            new_node = TreeNode(**kwargs)

            if last_node is not None:
                new_node.append(last_node)
            last_node = new_node
        return last_node

    def _setup_balanced_binary(self, kwargs_list):
        trees = []
        for kwargs in kwargs_list:
            trees.append(TreeNode(**kwargs))

        trees[0].extend([trees[2], trees[3]])
        trees[1].extend([trees[4], trees[5]])
        trees[6].extend([trees[0], trees[1]])
        return trees[6]

    def setUp(self):
        # Using the factory functions above, we will construct different tree
        # instances. Each tree is expected to serialize to the first newick
        # string in the list. Each string in the list is expected to
        # deserialize into an equivilent rotation of the constructed instance.
        tree_blank = (self._setup_tree([
            {}, {}, {}, {}, {}, {}
        ]), [
            "(,,(,));\n",
            "(,(,),);",
            "((,),,);",
            "   ((,[ this is a comment ])      ,    ,   )    ;  ",
            "((,[ i_can_do_this[0] or escape unmatched '[ ]),[more words],);",
        ])

        tree_leaves_named = (self._setup_tree([
            {'name': 'a_'},
            {'name': 'b'},
            {'name': 'c'},
            {'name': 'd'},
            {},
            {}
        ]), [
            "('a_',b,(c,d));\n",
            "(b,(c,d),'a_');",
            "(b\n,'a_'\n  ,(d \t,c) )  ;",
        ])

        tree_all_named = (self._setup_tree([
            {'name': 'a'},
            {'name': 'b'},
            {'name': 'c'},
            {'name': '[whaaat!\']'},
            {'name': 'e'},
            {'name': 'f'}
        ]), [
            "(a,b,(c,'[whaaat!'']')e)f;\n",
            "(b,(c,'[whaaat!'']')e,a)f;",
            "(b,[comment] \na,('[whaaat!'']',c)e)f;",
        ])

        tree_all_but_root_distances = (self._setup_tree([
            {'length': 0.1},
            {'length': 0.2},
            {'length': 0.3},
            {'length': 0.4},
            {'length': 0.5},
            {}
        ]), [
            "(:0.1,:0.2,(:0.3,:0.4):0.5);\n",
            "(:0.2,(:0.3,:0.4):0.5,:0.1);",
            "(:0.2,:0.1,(:0.4,:0.3):0.5);",
        ])

        tree_all_distances = (self._setup_tree([
            {'length': 0.1},
            {'length': 0.2},
            {'length': 0.3},
            {'length': 0.4},
            {'length': 0.5},
            {'length': 0.0}
        ]), [
            "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;\n",
            "(:0.2,(:0.3,:0.4):0.5,:0.1):0.0;",
            "(:0.2,\n:0.1,(:0.4,\n:0.3):0.5)\n:0.0;",
        ])

        tree_all_leaves_named_with_distances = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de d', 'length': 0.4},
            {'length': 0.5},
            {'length': 0.0}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de_d:0.4):0.5):0.0;\n",
            "('b_a''':0.2,(c:0.3,'de d':0.4):0.5,a:0.1):0.0;",
            "('b_a''':0.2,a:0.1,('de d'[why not]:0.4,c:0.3):0.5):0.0;",
        ])

        tree_all_leaves_named_with_distances_no_root = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de  d', 'length': 0.4},
            {'length': 0.5},
            {}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de__d:0.4):0.5);\n",
            "('b_a''':0.2\n[comment ahoy]\n,(c:0.3,'de  d':0.4):0.5,a:0.1);",
            "('b_a''':0.2,a:0.1,(de__d:0.4,c:0.3):0.5);"
        ])

        tree_all = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de\' d', 'length': 0.4},
            {'name': 'e', 'length': 0.5},
            {'name': 'f', 'length': 0.0}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de''_d:0.4)e:0.5)f:0.0;\n",
            "('b_a''':0.2,(c:0.3,de''_d:0.4)e:0.5,a:0.1)f:0.0;",
            "((de''_d:0.4, c:0.3)e:0.5, 'b_a''':0.2, a:0.1)f:0.0;"
        ])

        balanced_blank = (self._setup_balanced_binary([
            {}, {}, {}, {}, {}, {}, {}
        ]), [
            "((,),(,));\n",
        ])

        balanced_named = (self._setup_balanced_binary([
            {'name': 'a'},
            {'name': 'b'},
            {'name': 'c'},
            {'name': 'd'},
            {'name': 'e'},
            {'name': 'f'},
            {'name': 'g'}
        ]), [
            "((c,d)a,(e,f)b)g;\n",
        ])

        balanced_distances = (self._setup_balanced_binary([
            {'length': 1.0},
            {'length': 2.0},
            {'length': 3.0},
            {'length': 4.0},
            {'length': 5.0},
            {'length': 6.0},
            {'length': 0.0}
        ]), [
            "((:3.0,:4.0):1.0,(:5.0,:6.0):2.0):0.0;\n",
        ])

        blanaced_all = (self._setup_balanced_binary([
            {'name': 'a', 'length': 1.0},
            {'name': 'b', 'length': 2.0},
            {'name': 'c', 'length': 3.0},
            {'name': 'd', 'length': 4.0},
            {'name': 'e', 'length': 5.0},
            {'name': 'f:f\'f', 'length': 6.0},
            {'name': 'g', 'length': 0.0}
        ]), [
            "((c:3.0,d:4.0)a:1.0,(e:5.0,'f:f''f':6.0)b:2.0)g:0.0;\n",
        ])

        linked_list_blank = (self._setup_linked_list([
            {}, {}, {}, {}, {}
        ]), [
            "(((())));\n",
            "[(((())));](((())));",
            "[[(((())));](((())));](((())));\t\t\n"
        ])

        linked_list_named = (self._setup_linked_list([
            {'name': 'aaa'},
            {'name': 'b_a\''},
            {'name': 'c'},
            {'name': 'de d'},
            {'name': 'e'},
        ]), [
            "((((aaa)'b_a''')c)de_d)e;\n"
        ])

        inked_list_distances = (self._setup_linked_list([
            {'length': 0.4},
            {'length': 0.3},
            {'length': 0.2},
            {'length': 0.1},
            {'length': 0.0},
        ]), [
            "((((:0.4):0.3):0.2):0.1):0.0;\n",
            "((((:0.4)[not a label]:0.3):0.2):0.1):0.0;\t\t\n"
        ])

        linked_list_all = (self._setup_linked_list([
            {'name': 'a', 'length': 0.4},
            {'name': 'b_a\'', 'length': 0.3},
            {'name': 'c', 'length': 0.2},
            {'name': 'de d', 'length': 0.1},
            {'name': 'eee', 'length': 0.0},
        ]), [
            "((((a:0.4)'b_a''':0.3)c:0.2)de_d:0.1)eee:0.0;\n"
        ])

        single_empty = (TreeNode(), [";\n", "[comment about the root"
                                     " and its properties];"])
        single_named = (TreeNode(name='athing'), ["athing;\n"])
        single_distance = (TreeNode(length=200.0), [":200.0;\n"])
        single_all = (TreeNode(name='[a]', length=200.0), ["'[a]':200.0;\n"])

        self.trees_newick_lists = [
            tree_blank,
            tree_leaves_named,
            tree_all_named,
            tree_all_but_root_distances,
            tree_all_distances,
            tree_all_leaves_named_with_distances,
            tree_all_leaves_named_with_distances_no_root,
            tree_all,
            balanced_blank,
            balanced_named,
            balanced_distances,
            blanaced_all,
            linked_list_blank,
            linked_list_named,
            inked_list_distances,
            linked_list_all,
            single_empty,
            single_named,
            single_distance,
            single_all
        ]

        # Invalid newick strings and list of error fragments that should be
        # a part of the error message when read.
        self.invalid_newicks = [
            ("", ['root']),
            ("This is not a newick file.", ['whitespace', 'label']),
            ("((();", ['Parenthesis', 'unbalanced']),
            ("(,,,)(,);\n", ['unnested', 'children']),
            ("(()());", ['unnested', 'children']),
            ("(():,,)", ['length']),
            ("[][[]('comment is the gotcha':0.2,,);", ['unbalanced', 'root']),
            ("#SampleID\tHeaderA\tHeaderB\n0\t'yellow'\t0.45;", ['whitespace',
                                                                 'label']),
            ("))();", ['Parenthesis', 'unbalanced']),
            ("((,,),((,,));", ['Parenthesis', 'unbalanced']),
            ("\n".join([",".join(str(i) for i in range(100))
                       for _ in range(100)]), ['whitespace', 'label'])
        ]

    def test_newick_to_tree_node_valid_files(self):
        for tree, newicks in self.trees_newick_lists:
            for newick in newicks:
                fh = io.StringIO(newick)
                read_tree = _newick_to_tree_node(fh)

                self._assert_equal(tree, read_tree)

                fh.close()

    def test_newick_to_tree_node_invalid_files(self):
        for invalid, error_fragments in self.invalid_newicks:
            fh = io.StringIO(invalid)
            with self.assertRaises(NewickFormatError) as cm:
                _newick_to_tree_node(fh)
            for frag in error_fragments:
                self.assertIn(frag, str(cm.exception))
            fh.close()

    def test_tree_node_to_newick(self):
        for tree, newicks in self.trees_newick_lists:
            newick = newicks[0]
            fh = io.StringIO()
            _tree_node_to_newick(tree, fh)

            self.assertEqual(newick, fh.getvalue())

            fh.close()

    def test_roundtrip(self):
        for tree, newicks in self.trees_newick_lists:
            newick = newicks[0]
            fh = io.StringIO(newick)
            tree = _newick_to_tree_node(fh)
            fh2 = io.StringIO()
            _tree_node_to_newick(tree, fh2)
            fh2.seek(0)
            tree2 = _newick_to_tree_node(fh2)

            self.assertEqual(newick, fh2.getvalue())
            self._assert_equal(tree, tree2)

            fh.close()
            fh2.close()

    def test_newick_to_tree_node_convert_underscores(self):
        fh = io.StringIO('(_:0.1, _a, _b)__;')
        tree = _newick_to_tree_node(fh, convert_underscores=False)
        fh2 = io.StringIO()
        _tree_node_to_newick(tree, fh2)
        self.assertEqual(fh2.getvalue(), "('_':0.1,'_a','_b')'__';\n")
        fh2.close()
        fh.close()

    def test_newick_sniffer_valid_files(self):
        for _, newicks in self.trees_newick_lists:
            for newick in newicks:
                fh = io.StringIO(newick)
                self.assertEqual(_newick_sniffer(fh), (True, {}))
                fh.close()

    def test_newick_sniffer_invalid_files(self):
        for invalid, _ in self.invalid_newicks:
            fh = io.StringIO(invalid)
            self.assertEqual(_newick_sniffer(fh), (False, {}))
            fh.close()


class TestNewickBPTree(unittest.TestCase):
    # A diverse battery of valid newick strings covering multifurcation,
    # single-descendant chains, unnamed/unlengthed nodes, quoted labels
    # containing structural characters, and exponent lengths.
    battery = [
        "((a:2,b):1,(c:4,d)y:20,e)r;",
        "(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;",
        "(((a)b)c,((d)e)f)r;",
        "((a,b)c,d,(e))r;",
        "((a,(b,c):5)'d','e; foo':10,((f))g)r;",
        "(('foo,bar':1,baz:2)x:3)r;",
        "(('foo(b)ar':1,baz:2)x:3)r;",
        '((foo"bar":1,baz:2)x:3)r;',
        "((b:3)a:2)root:1;",
        "((A:.01, B:.01)D:.01, C:.01)r;",
        "(a:1e-3,b:2.5E2)r;",
        "((:0.25,:0.5):0.0,x)r;",
    ]

    def _read_bp(self, s):
        return _newick_to_bp(io.StringIO(s))

    def test_exact_arrays(self):
        # A subset of the balanced-parentheses regression oracle, driving the
        # registry reader that delegates to the compiled parser.
        in_ = "((a:2,b):1,(c:4,d)y:20,e)r;"
        exp_bp = [1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0]
        exp_n = ['r', None, 'a', None, 'b', None, None, 'y', 'c', None, 'd',
                 None, None, 'e', None, None]
        exp_l = [0, 1, 2, 0, 0, 0, 0, 20, 4, 0, 0, 0, 0, 0, 0, 0]

        obs = self._read_bp(in_)
        npt.assert_equal(obs.data, np.asarray(exp_bp, dtype=np.uint8))
        for i, (e_n, e_l) in enumerate(zip(exp_n, exp_l)):
            self.assertEqual(obs.name(i), e_n)
            self.assertEqual(obs.length(i), e_l)

    def test_exact_arrays_singledesc(self):
        in_ = "(((a)b)c,((d)e)f)r;"
        exp_bp = [1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]
        exp_n = ['r', 'c', 'b', 'a', None, None, None, 'f', 'e', 'd', None,
                 None, None, None]
        obs = self._read_bp(in_)
        npt.assert_equal(obs.data, np.asarray(exp_bp, dtype=np.uint8))
        for i, e_n in enumerate(exp_n):
            self.assertEqual(obs.name(i), e_n)

    def test_write_roundtrip(self):
        for s in self.battery:
            bp = self._read_bp(s)
            buf = io.StringIO()
            _bp_to_newick(bp, buf)
            buf.seek(0)
            obs = self._read_bp(buf.getvalue())
            npt.assert_equal(obs.data, bp.data, err_msg="round-trip topo %r" % s)
            for i in range(bp.data.size):
                self.assertEqual(obs.name(i), bp.name(i))
                self.assertAlmostEqual(obs.length(i), bp.length(i))

    def test_read_via_registry_and_sniff(self):
        s = "((a:2,b):1,(c:4,d)y:20,e)r;"
        # explicit format
        exp = self._read_bp(s)
        obs = BPTree.read([s], format="newick")
        npt.assert_equal(obs.data, exp.data)
        # format inference (the newick sniffer dispatches to the BPTree reader
        # because into=BPTree)
        obs2 = BPTree.read([s])
        npt.assert_equal(obs2.data, exp.data)

    def test_write_via_registry_default_format(self):
        s = "((a:2,b):1,(c:4,d)y:20,e)r;"
        bp = BPTree.read([s], format="newick")
        buf = io.StringIO()
        bp.write(buf)  # default_write_format == "newick"
        buf.seek(0)
        obs = BPTree.read(buf)
        npt.assert_equal(obs.data, bp.data)

    def test_invalid_missing_semicolon(self):
        with self.assertRaises(NewickFormatError):
            _newick_to_bp(io.StringIO("((a,b)c,d)r"))

    def test_invalid_single_node(self):
        with self.assertRaises(NewickFormatError):
            _newick_to_bp(io.StringIO("i:1;"))

    def test_invalid_unbalanced(self):
        with self.assertRaises(NewickFormatError):
            _newick_to_bp(io.StringIO("((a,b)c,d));"))


if __name__ == '__main__':
    unittest.main()
