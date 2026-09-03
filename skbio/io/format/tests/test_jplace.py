# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import re
import json
import unittest

import pandas as pd
import pandas.testing as pdt

import skbio
import skbio.io
from skbio.io import JplaceFormatError
from skbio.tree import TreeNode, BPTree
from skbio.tree.bp import parse_jplace
from skbio.io.format.jplace import _jplace_sniffer


def _fixture_path():
    # Reuse the jplace fixture that ships with the BP tests rather than
    # duplicating a ~10 KB file under this test directory.
    return os.path.join(
        os.path.dirname(skbio.__file__),
        "tree",
        "tests",
        "bp",
        "data",
        "200",
        "placement.jplace",
    )


COLUMNS = [
    "fragment",
    "edge_num",
    "likelihood",
    "like_weight_ratio",
    "distal_length",
    "pendant_length",
]


def _edges(bp):
    return [bp.edge(i) for i in range(bp.data.size)]


def _synthetic_jplace(tree):
    # A minimal, valid jplace document with a single placement on edge 0, used
    # to exercise reference-tree Newick grammar (comments, quoting, underscores).
    return json.dumps(
        {
            "tree": tree,
            "placements": [{"p": [[0, -1.0, 1.0, 0.0, 0.0]], "n": ["frag"]}],
            "fields": ["edge_num", "likelihood", "like_weight_ratio",
                       "distal_length", "pendant_length"],
            "version": 3,
            "metadata": {},
        }
    )


class TestJplaceReader(unittest.TestCase):
    """The jplace registry reader returns the reference tree as a BPTree and
    drops the placement records."""

    def setUp(self):
        self.path = _fixture_path()
        with open(self.path) as fh:
            self.jplacedata = fh.read()
        no_edge_numbers = re.sub(r"{\d+}", "", json.loads(self.jplacedata)["tree"])
        self.tree = TreeNode.read([no_edge_numbers])  # topology oracle

    def _read_bp(self, data, **kwargs):
        return skbio.io.read([data], into=BPTree, format="jplace", **kwargs)

    def test_read_into_bptree(self):
        bp = self._read_bp(self.jplacedata)
        self.assertIsInstance(bp, BPTree)
        self.assertEqual(TreeNode.from_bptree(bp).compare_rfd(self.tree), 0)

    def test_read_preserves_edge_numbers(self):
        # The {} edge numbers are parsed onto the tree (placements are dropped).
        bp = self._read_bp(self.jplacedata)
        self.assertTrue(any(e != 0 for e in _edges(bp)))

    def test_read_from_path(self):
        bp = skbio.io.read(self.path, into=BPTree, format="jplace")
        self.assertIsInstance(bp, BPTree)
        self.assertEqual(TreeNode.from_bptree(bp).compare_rfd(self.tree), 0)

    def test_read_via_descriptor(self):
        bp = BPTree.read(self.path, format="jplace")
        self.assertIsInstance(bp, BPTree)

    def test_read_via_descriptor_inferred_format(self):
        # The sniffer routes a jplace file to the jplace reader for into=BPTree.
        bp = BPTree.read(self.path)
        self.assertIsInstance(bp, BPTree)
        self.assertEqual(TreeNode.from_bptree(bp).compare_rfd(self.tree), 0)

    def test_square_braces_treated_as_comment(self):
        # On the grammar branch [] is a Newick comment: converting the canonical
        # {N} to [N] drops the edge numbers (edges become 0) while leaving the
        # topology intact.
        data = json.loads(self.jplacedata)
        data["tree"] = re.sub(r"{(\d+)}", r"[\1]", data["tree"])
        bp = self._read_bp(json.dumps(data))
        self.assertEqual(TreeNode.from_bptree(bp).compare_rfd(self.tree), 0)
        self.assertTrue(all(e == 0 for e in _edges(bp)))

    def test_missing_required_key_raises(self):
        data = json.loads(self.jplacedata)
        del data["fields"]
        with self.assertRaises(JplaceFormatError):
            self._read_bp(json.dumps(data))

    def test_missing_n_raises(self):
        data = json.loads(self.jplacedata)
        for placement in data["placements"]:
            placement.pop("n", None)
        with self.assertRaises(JplaceFormatError):
            self._read_bp(json.dumps(data))

    def test_fields_without_edge_num_raises(self):
        data = json.loads(self.jplacedata)
        data["fields"] = [f for f in data["fields"] if f != "edge_num"]
        with self.assertRaises(JplaceFormatError):
            self._read_bp(json.dumps(data))


class TestJplaceWriter(unittest.TestCase):
    """Writing a BPTree emits a jplace document carrying the tree (with its
    edge numbers) and an empty placements list."""

    def setUp(self):
        self.path = _fixture_path()
        self.bp = skbio.io.read(self.path, into=BPTree, format="jplace")

    def _round_trip(self, bp):
        buf = io.StringIO()
        skbio.io.write(bp, into=buf, format="jplace")
        raw = buf.getvalue()
        return raw, skbio.io.read([raw], into=BPTree, format="jplace")

    def test_write_round_trip_tree_and_edges(self):
        _, back = self._round_trip(self.bp)
        self.assertEqual(
            TreeNode.from_bptree(self.bp).compare_rfd(TreeNode.from_bptree(back)),
            0,
        )
        self.assertEqual(_edges(back), _edges(self.bp))

    def test_write_output_is_valid_jplace(self):
        raw, _ = self._round_trip(self.bp)
        doc = json.loads(raw)
        self.assertEqual(
            sorted(doc.keys()),
            ["fields", "metadata", "placements", "tree", "version"],
        )
        self.assertEqual(doc["placements"], [])
        self.assertIn("edge_num", doc["fields"])
        self.assertEqual(doc["version"], 3)
        # the embedded tree carries {} edge numbers
        self.assertRegex(doc["tree"], r"\{\d+\}")
        # and the document is recognized as jplace
        self.assertEqual(skbio.io.sniff([raw])[0], "jplace")

    def test_write_via_descriptor(self):
        buf = io.StringIO()
        self.bp.write(buf, format="jplace")
        doc = json.loads(buf.getvalue())
        self.assertTrue(doc["tree"].rstrip().endswith(";"))
        self.assertEqual(doc["placements"], [])

    def test_write_edgeless_tree_is_degenerate_but_readable(self):
        # A tree with no edge numbers writes {0} on every edge (no unique
        # numbers); topology still round-trips.
        bp = skbio.io.read(["((a,b)c,d)r;"], into=BPTree, format="newick")
        buf = io.StringIO()
        bp.write(buf, format="jplace")
        doc = json.loads(buf.getvalue())
        self.assertNotIn("{1}", doc["tree"])
        back = skbio.io.read([buf.getvalue()], into=BPTree, format="jplace")
        self.assertEqual(
            TreeNode.from_bptree(bp).compare_rfd(TreeNode.from_bptree(back)),
            0,
        )


class TestJplaceReferenceTreeGrammar(unittest.TestCase):
    """Reference-tree Newick grammar seen through the jplace BPTree path."""

    def _read(self, tree):
        return skbio.io.read([_synthetic_jplace(tree)], into=BPTree,
                             format="jplace")

    def _names(self, bp):
        return {bp.name(i) for i in range(bp.data.size)}

    def test_underscores_preserved(self):
        # jplace taxa are accessions; '_' is significant and must not become a
        # space (the jplace path reads with convert_underscores=False).
        bp = self._read("((A_1:.1{0}, B_2:.1{1})C:.1{2}, D:.1{3}){4};")
        names = self._names(bp)
        self.assertIn("A_1", names)
        self.assertIn("B_2", names)
        self.assertNotIn("A 1", names)

    def test_comment_in_reference_tree_ignored(self):
        # A [comment] on an edge is stripped; the {} edge numbers still resolve.
        bp = self._read("((A:.1{0}[note], B:.1{1})C:.1{2}, D:.1{3}){4};")
        self.assertIn(0, {bp.edge(i) for i in range(bp.data.size)})
        self.assertIn("A", self._names(bp))

    def test_special_char_label_round_trip(self):
        # A whitespace label is quoted in the written tree and re-reads intact.
        bp = self._read("(('a b':.1{0}, B:.1{1})C:.1{2}, D:.1{3}){4};")
        self.assertIn("a b", self._names(bp))
        buf = io.StringIO()
        bp.write(buf, format="jplace")
        self.assertIn("'a b'", json.loads(buf.getvalue())["tree"])
        back = skbio.io.read([buf.getvalue()], into=BPTree, format="jplace")
        self.assertIn("a b", self._names(back))


class TestJplaceSniffer(unittest.TestCase):
    def setUp(self):
        self.path = _fixture_path()

    def test_sniffer_positive(self):
        self.assertEqual(skbio.io.sniff(self.path), ("jplace", {}))

    def test_sniffer_rejects_newick(self):
        fh = io.StringIO("((a:1,b:2)c:3)r;")
        self.assertEqual(_jplace_sniffer(fh), (False, {}))

    def test_sniffer_rejects_json_missing_keys(self):
        fh = io.StringIO(json.dumps({"tree": "(a,b)r;", "placements": []}))
        self.assertEqual(_jplace_sniffer(fh), (False, {}))

    def test_sniffer_rejects_non_json(self):
        fh = io.StringIO("this is not json at all")
        self.assertEqual(_jplace_sniffer(fh), (False, {}))

    def test_newick_does_not_sniff_as_jplace(self):
        # A newick file must be identified as newick, never jplace.
        fh = io.StringIO("((D, E)B, (F, G)C)A;")
        self.assertEqual(skbio.io.sniff(fh)[0], "newick")

    def test_jplace_does_not_sniff_as_newick(self):
        # The embedded newick tree must not cause the whole jplace document to
        # be mistaken for newick.
        self.assertEqual(skbio.io.sniff(self.path)[0], "jplace")


class TestParseJplace(unittest.TestCase):
    """The iow ``parse_jplace(data) -> (DataFrame, BPTree)`` helper remains the
    way to access placements alongside the reference tree."""

    def setUp(self):
        self.path = _fixture_path()
        with open(self.path) as fh:
            self.jplacedata = fh.read()
        no_edge_numbers = re.sub(r"{\d+}", "", json.loads(self.jplacedata)["tree"])
        self.tree = TreeNode.read([no_edge_numbers])

    def test_returns_df_and_tree(self):
        df, tree = parse_jplace(self.jplacedata)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIsInstance(tree, BPTree)

    def test_simple(self):
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
        exp_df = pd.DataFrame(exp_df, columns=COLUMNS)
        df, tree = parse_jplace(self.jplacedata)
        pdt.assert_frame_equal(df, exp_df)
        self.assertEqual(TreeNode.from_bptree(tree).compare_rfd(self.tree), 0)

    def test_multiple_per_fragment(self):
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
        exp_df = pd.DataFrame(exp_df, columns=COLUMNS)

        data = json.loads(self.jplacedata)
        keep = []
        for placement in data["placements"]:
            if placement["n"][0] == "82":
                keep.append(placement)
            elif placement["n"][0] == "99":
                placement["p"].append([309, 0.04520741687623886, 1,
                                       0.11020044356641526,
                                       0.00550337922097477])
                keep.append(placement)
            elif placement["n"][0] == "55":
                placement["p"].append([138, 0.09563944060686769, 10,
                                       0.014593217782258146,
                                       0.04537214236560885])
                keep.append(placement)
        data["placements"] = keep

        df, tree = parse_jplace(json.dumps(data))
        pdt.assert_frame_equal(df, exp_df)
        self.assertEqual(TreeNode.from_bptree(tree).compare_rfd(self.tree), 0)

    def test_matches_registry_tree(self):
        _, tree = parse_jplace(self.jplacedata)
        bp = skbio.io.read([self.jplacedata], into=BPTree, format="jplace")
        self.assertEqual(
            TreeNode.from_bptree(tree).compare_rfd(TreeNode.from_bptree(bp)),
            0,
        )


if __name__ == "__main__":
    unittest.main()
