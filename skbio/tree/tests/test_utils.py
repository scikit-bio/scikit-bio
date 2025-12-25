# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy.testing as npt

from skbio.tree import TreeNode
from skbio.tree._exception import DuplicateNodeError, MissingNodeError
from skbio.tree._utils import _validate_taxa_and_tree


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


if __name__ == "__main__":
    main()
