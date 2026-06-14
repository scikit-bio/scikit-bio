# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import json
import io
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from skbio.tree import parse_newick, to_skbio_treenode, write_newick, parse_jplace


def get_data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)


class NewickTests(TestCase):
    def test_parse_newick_simple_edge_numbers(self):
        in_ = '((A:.01{0}, B:.01{1})D:.01{3}, C:.01{4}) {5};'
        exp_sk = '((A:.01, B:.01)D:.01, C:.01);'
        obs = parse_newick(in_)
        obs_sk = to_skbio_treenode(obs)
        exp_sk = skbio.TreeNode.read([exp_sk])
        self.assertEqual(obs_sk.compare_rfd(exp_sk), 0)

        self.assertEqual(obs.edge(2), 0)
        self.assertEqual(obs.edge(4), 1)
        self.assertEqual(obs.edge(1), 3)
        self.assertEqual(obs.edge(7), 4)
        self.assertEqual(obs.edge(0), 5)
        self.assertEqual(obs.edge_from_number(0), 2)
        self.assertEqual(obs.edge_from_number(1), 4)
        self.assertEqual(obs.edge_from_number(3), 1)
        self.assertEqual(obs.edge_from_number(4), 7)
        self.assertEqual(obs.edge_from_number(5), 0)

    def _compare_newick(self, obs, exp):
        a = skbio.TreeNode.read([obs])
        b = skbio.TreeNode.read([exp])
        self.assertEqual(a.compare_rfd(b), 0)
        npt.assert_equal(a.cophenet().data, b.cophenet().data)

    def test_write_newick_cases(self):
        tests = ['((foo"bar":1,baz:2)x:3)r;',
                 "(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;",
                 "(((a)b)c,((d)e)f)r;",
                 "((a,(b,c):5)'d','e; foo':10,((f))g)r;"]

        for test in tests:
            buf = io.StringIO()
            obs = write_newick(parse_newick(test), buf, False)
            buf.seek(0)
            obs = buf.read()
            self._compare_newick(obs, test)

    def test_write_newick_edges(self):
        test_a = '((foo"bar":1{0},baz:2{1})x:3{2})r;'
        test_b = "(((a)b)c,((d)e)f)r;"

        buf = io.StringIO()
        obs = write_newick(parse_newick(test_a), buf, True)
        buf.seek(0)
        obs = to_skbio_treenode(parse_newick(buf.read()))
        self.assertEqual(obs.find('foo"bar"').edge_num, 0)
        self.assertEqual(obs.find('baz').edge_num, 1)
        self.assertEqual(obs.find('x').edge_num, 2)

        buf = io.StringIO()
        obs = write_newick(parse_newick(test_b), buf, True)
        buf.seek(0)
        obs = to_skbio_treenode(parse_newick(buf.read()))
        for o in obs.traverse():
            self.assertEqual(o.edge_num, 0)

    def test_parse_newick_singlenode_bug(self):
        test = 'i:1;'
        with self.assertRaises(ValueError):
            parse_newick(test)

    def test_parse_newick_no_semicolon_bug(self):
        test = "((h:1, i:1, j:1, k:1, l: 1),(e:1,f:1),(n:1,o:1,p:1))a:1"
        with self.assertRaises(ValueError):
            parse_newick(test)

        test = "((h:1, i:1, j:1, k:1, l: 1),(e:1,f:1),(n:1,o:1,p:1))a:1;\n"
        parse_newick(test)

    def test_write_newick_underscore_bug(self):
        test = "(((a)b)'c_foo',((d)e)f)r;"
        buf = io.StringIO()
        obs = write_newick(parse_newick(test), buf, False)
        buf.seek(0)
        self.assertIn("'c_foo'", test)

    def test_parse_newick_nested_quotes(self):
        in_ = '((foo"bar":1,baz:2)x:3)r;'
        exp = skbio.TreeNode.read([in_])
        obs = to_skbio_treenode(parse_newick(in_))
        self.assertEqual(obs.compare_rfd(exp), 0.0)

    def test_parse_newick_with_commas(self):
        in_ = "(('foo,bar':1,baz:2)x:3)r;"
        exp = skbio.TreeNode.read([in_])
        obs = to_skbio_treenode(parse_newick(in_))
        self.assertEqual(obs.compare_rfd(exp), 0.0)

    def test_parse_newick_with_parens(self):
        in_ = "(('foo(b)ar':1,baz:2)x:3)r;"
        exp = skbio.TreeNode.read([in_])
        obs = to_skbio_treenode(parse_newick(in_))
        self.assertEqual(obs.compare_rfd(exp), 0.0)

    def test_parse_newick(self):
        in_ = "((a:2,b):1,(c:4,d)y:20,e)r;"

        exp_bp = [1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0]
        exp_n = ['r', None, 'a', None, 'b', None, None, 'y', 'c', None, 'd',
                 None, None, 'e', None, None]
        exp_l = [0, 1, 2, 0, 0, 0, 0, 20, 4, 0, 0, 0, 0, 0, 0, 0]

        obs_bp = parse_newick(in_)

        npt.assert_equal(obs_bp.B, np.asarray(exp_bp, dtype=bool))

        for i, (e_n, e_l) in enumerate(zip(exp_n, exp_l)):
            self.assertEqual(obs_bp.name(i), e_n)
            self.assertEqual(obs_bp.length(i), e_l)

    def test_parse_newick_complex(self):
        in_ = "(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;"
        exp_bp = [1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1,
                  1, 0, 1, 0, 0, 0, 0, 1, 0, 0]
        exp_n = ['r', 'l', 'c', 'a', None, 'b', None, None, 'd', None, None,
                 'e', None, None, 'k', 'f', None, 'g', None, 'j', 'h', None,
                 'i', None, None, None, None, 'm', None, None]
        exp_l = [0, 0, 6, 1, 0, 2.5, 0, 0, 8, 0, 0, 0, 0, 0, 1.2, 0, 0, 0, 0,
                 1, 1, 0, 2, 0, 0, 0, 0, 2, 0, 0]

        obs_bp = parse_newick(in_)

        npt.assert_equal(obs_bp.B, np.asarray(exp_bp, dtype=bool))

        for i, (e_n, e_l) in enumerate(zip(exp_n, exp_l)):
            self.assertEqual(obs_bp.name(i), e_n)
            self.assertEqual(obs_bp.length(i), e_l)

    def test_parse_newick_singledesc(self):
        in_ = "(((a)b)c,((d)e)f)r;"
        exp_bp = [1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]
        exp_n = ['r', 'c', 'b', 'a', None, None, None, 'f', 'e', 'd', None,
                 None, None, None]

        obs_bp = parse_newick(in_)

        npt.assert_equal(obs_bp.B, np.asarray(exp_bp, dtype=bool))

        for i, e_n in enumerate(exp_n):
            self.assertEqual(obs_bp.name(i), e_n)

    def test_parse_newick_unnamed_singledesc(self):
        in_ = "((a,b)c,d,(e))r;"

        exp_bp = [1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0]
        exp_n = ['r', 'c', 'a', None, 'b', None, None, 'd', None, None, 'e',
                 None, None, None]

        obs_bp = parse_newick(in_)

        npt.assert_equal(obs_bp.B, np.asarray(exp_bp, dtype=bool))

        for i, e_n in enumerate(exp_n):
            self.assertEqual(obs_bp.name(i), e_n)

    def test_parse_newick_name_with_semicolon(self):
        in_ = "((a,(b,c):5)'d','e; foo':10,((f))g)r;"

        exp_bp = [1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0]
        exp_n = ['r', 'd', 'a', None, None, 'b', None, 'c', None, None, None,
                 'e; foo', None, 'g', None, 'f', None, None, None, None]
        exp_l = [0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0]

        obs_bp = parse_newick(in_)

        npt.assert_equal(obs_bp.B, np.asarray(exp_bp, dtype=bool))

        for i, (e_n, e_l) in enumerate(zip(exp_n, exp_l)):
            self.assertEqual(obs_bp.name(i), e_n)
            self.assertEqual(obs_bp.length(i), e_l)


class JPlaceParseTests(TestCase):
    def setUp(self):
        self.jplacedata = open(get_data_path('200/placement.jplace')).read()
        no_edge_numbers = re.sub(r"{\d+}", '',
                                 json.loads(self.jplacedata)['tree'])
        self.tree = skbio.TreeNode.read([no_edge_numbers])

    def test_place_jplace_square_braces(self):
        jplacedata = json.loads(self.jplacedata)
        treestr = jplacedata['tree']
        treestr = re.sub(r"{(\d+)}", r"[\1]", treestr)
        jplacedata['tree'] = treestr
        jplacedata = json.dumps(jplacedata)

        exp_tree = self.tree
        obs_df, obs_tree = parse_jplace(jplacedata)
        obs_tree = to_skbio_treenode(obs_tree)
        self.assertEqual(obs_tree.compare_rfd(exp_tree), 0)
        for n in obs_tree.traverse(include_self=False):
            self.assertTrue(n.edge_num >= 0)

    def test_parse_jplace_simple(self):
        columns = ['fragment', 'edge_num', 'likelihood', 'like_weight_ratio',
                   'distal_length', 'pendant_length']
        exp_df = [["82", 361, 0.01013206496780672, 1, 0.02652932626620403,
                   0.039354548684623215],
                  ["99", 308, 0.04520741687623886, 1, 0.11020044356641526,
                   0.06550337922097477],
                  ["43", 309, 0.04054866161921744, 1, 0.010712923050783987,
                   0.020946988900520196],
                  ["195", 277, 0.01918907908397749, 1, 0.03065741838803451,
                   0.04513513498399864],
                  ["162", 55, 0.01758935282545493, 1, 0.0033199487685078776,
                   0.05388735804976052],
                  ["56", 81, 0.2366882303770561, 1, 0.04172580852519453,
                   0.0007060238727097983],
                  ["91", 105, 0.0001863393767883581, 1, 0.04578898721138839,
                   0.08655004339151215],
                  ["174", 89, 0.01216463967379211, 1, 0.04707020642820376,
                   0.045206727542450205],
                  ["5", 143, 0.012162345471765756, 1, 0.023797389484252734,
                   0.10447375403452556],
                  ["55", 139, 0.09563944060686769, 1, 0.014593217782258146,
                   0.04537214236560885]]
        exp_df = pd.DataFrame(exp_df, columns=columns)
        exp_tree = self.tree
        obs_df, obs_tree = parse_jplace(self.jplacedata)
        obs_tree = to_skbio_treenode(obs_tree)
        pdt.assert_frame_equal(obs_df, exp_df)
        self.assertEqual(obs_tree.compare_rfd(exp_tree), 0)

    def test_parse_jplace_multiple_per_fragment(self):
        columns = ['fragment', 'edge_num', 'likelihood', 'like_weight_ratio',
                   'distal_length', 'pendant_length']
        exp_df = [["82", 361, 0.01013206496780672, 1, 0.02652932626620403,
                   0.039354548684623215],
                  ["99", 308, 0.04520741687623886, 1, 0.11020044356641526,
                   0.06550337922097477],
                  ["99", 309, 0.04520741687623886, 1, 0.11020044356641526,
                   0.00550337922097477],
                  ["55", 139, 0.09563944060686769, 1, 0.014593217782258146,
                   0.04537214236560885],
                  ["55", 138, 0.09563944060686769, 10, 0.014593217782258146,
                   0.04537214236560885]]
        exp_df = pd.DataFrame(exp_df, columns=columns)

        data = json.loads(self.jplacedata)
        keep = []
        for placement in data['placements']:
            if placement['n'][0] == '82':
                keep.append(placement)
            elif placement['n'][0] == '99':
                placement['p'].append([309, 0.04520741687623886, 1,
                                       0.11020044356641526,
                                       0.00550337922097477])
                keep.append(placement)
            elif placement['n'][0] == '55':
                placement['p'].append([138, 0.09563944060686769, 10,
                                       0.014593217782258146,
                                       0.04537214236560885])
                keep.append(placement)
        data['placements'] = keep
        data = json.dumps(data)
        exp_tree = self.tree
        obs_df, obs_tree = parse_jplace(data)
        obs_tree = to_skbio_treenode(obs_tree)
        pdt.assert_frame_equal(obs_df, exp_df)
        self.assertEqual(obs_tree.compare_rfd(exp_tree), 0)

    def test_parse_newick_linear_tree(self):
        test = '((b:3)a:2)root:1;'
        topology = parse_newick(test)
        skbio_tree = to_skbio_treenode(topology)
        self.assertEqual(skbio_tree.name, "root")
        self.assertEqual([n.name for n in skbio_tree.children], ["a"])
        self.assertEqual([n.name for n in skbio_tree.non_tips()], ["a"])
        self.assertEqual([n.name for n in skbio_tree.tips()], ["b"])


if __name__ == '__main__':
    main()
