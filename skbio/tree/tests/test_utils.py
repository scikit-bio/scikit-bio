# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio import TreeNode
from skbio.tree._exception import DuplicateNodeError, MissingNodeError
from skbio.tree._utils import _tree_to_lnkmat, _validate_taxa_and_tree


class UtilsTests(TestCase):

    def test_validate_taxa_and_tree(self):
        # basic valid input
        tree = TreeNode.read([
            "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75)"
            ":1.25):0.0)root;"])
        taxa = ["OTU1", "OTU2", "OTU3"]
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree))

        # all tips observed
        taxa = ["OTU1", "OTU2", "OTU3", "OTU4", "OTU5"]
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree))

        # no tips observed
        taxa = []
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree))

    def test_validate_taxa_and_tree_invalid_input(self):
        # tree has duplicated tip ids
        tree = TreeNode.read([
            "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU2:0.75)"
            ":1.25):0.0)root;"])
        taxa = ["OTU1", "OTU2", "OTU3"]
        self.assertRaises(DuplicateNodeError, _validate_taxa_and_tree, taxa, tree)
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree, unique=False))

        # unrooted tree as input
        tree = TreeNode.read(["((OTU1:0.1,OTU2:0.2):0.3,OTU3:0.5,OTU4:0.7);"])
        taxa = ["OTU1", "OTU2", "OTU3"]
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree))
        self.assertRaises(ValueError, _validate_taxa_and_tree, taxa, tree, rooted=True)

        # taxa has duplicated ids
        tree = TreeNode.read([
            "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75)"
            ":1.25):0.0)root;"])
        taxa = ["OTU1", "OTU2", "OTU2"]
        self.assertRaises(ValueError, _validate_taxa_and_tree, taxa, tree)
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree, unique=False))

        # tree with no branch lengths
        tree = TreeNode.read(["((((OTU1,OTU2),OTU3)),(OTU4,OTU5));"])
        taxa = ["OTU1", "OTU2", "OTU3"]
        self.assertIsNone(_validate_taxa_and_tree(taxa, tree))
        self.assertRaises(ValueError, _validate_taxa_and_tree, taxa, tree,
                          lengths=True)

        # tree missing some branch lengths
        tree = TreeNode.read([
            "(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75)"
            ":1.25):0.0)root;"])
        taxa = ["OTU1", "OTU2", "OTU3"]
        self.assertRaises(ValueError, _validate_taxa_and_tree, taxa, tree,
                          lengths=True)

        # taxa not present in tree
        tree = TreeNode.read([
            "(((((OTU1:0.25,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75)"
            ":1.25):0.0)root;"])
        taxa = ["OTU1", "OTU2", "OTU32"]
        self.assertRaises(MissingNodeError, _validate_taxa_and_tree, taxa, tree)

    def test_tree_to_lnkmat(self):
        tree = TreeNode.read(["(((a,b),(c,d)),(e,(f,g)));"])

        # default taxa
        obs = _tree_to_lnkmat(tree)
        exp = np.array([[ 0,  1],
                        [ 2,  3],
                        [ 7,  8],
                        [ 5,  6],
                        [ 4, 10],
                        [ 9, 11]], dtype=np.intp)
        self.assertEqual(obs.dtype, np.intp)
        npt.assert_array_equal(obs, exp)

        # taxa in postorder
        taxa = "abcdefg"
        self.assertEqual(taxa, "".join(x.name for x in tree.tips()))
        obs = _tree_to_lnkmat(tree, taxa)
        npt.assert_array_equal(obs, exp)

        # supply output buffer
        out = np.zeros((6, 2), dtype=np.float32)
        obs = _tree_to_lnkmat(tree, out=out)
        self.assertIs(obs, out)
        self.assertEqual(obs.dtype, np.float32)
        npt.assert_array_equal(obs, exp)

        # taxa in alternative order
        taxa = "cebfagd"
        obs = _tree_to_lnkmat(tree, taxa)
        exp = np.array([[ 4,  2],
                        [ 0,  6],
                        [ 7,  8],
                        [ 3,  5],
                        [ 1, 10],
                        [ 9, 11]])
        npt.assert_array_equal(obs, exp)

        msg = "Taxa contain duplicates."
        with self.assertRaises(ValueError) as cm:
            _tree_to_lnkmat(tree, "abcdcba")
        self.assertEqual(str(cm.exception), msg)

        msg = "`out` must have a shape of (6, 2)."
        with self.assertRaises(ValueError) as cm:
            _tree_to_lnkmat(tree, out=np.empty((5, 3)))
        self.assertEqual(str(cm.exception), msg)

        msg = "Tip name 'f' is absent from taxa."
        with self.assertRaises(ValueError) as cm:
            _tree_to_lnkmat(tree, "abcde")
        self.assertEqual(str(cm.exception), msg)

        msg = "One or more taxa are absent from the tree."
        with self.assertRaises(ValueError) as cm:
            _tree_to_lnkmat(tree, "abcdefgx")
        self.assertEqual(str(cm.exception), msg)

        msg = "Tree must be strictly bifurcating."
        for tree in (TreeNode.read(["((a,b,c),d);"]),
                     TreeNode.read(["((a),(b,c));"])):
            with self.assertRaises(ValueError) as cm:
                _tree_to_lnkmat(tree)
            self.assertEqual(str(cm.exception), msg)

        msg = "Tree contains duplicate tip names."
        for tree, taxa in zip(
            (TreeNode.read(["((a,a),b);"]),
             TreeNode.read(["((a,b),(a,c));"]),
             TreeNode.read(["((a,b),(a,c));"])),
            ("abc", "abcd", "abc"),
        ):
            with self.assertRaises(ValueError) as cm:
                _tree_to_lnkmat(tree, taxa)
            self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    main()
