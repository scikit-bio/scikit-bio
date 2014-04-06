#!/usr/bin/env python
"""Unit tests for the skbio.util.trie module"""
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from itertools import izip

from skbio.util.trie import CompressedTrie, _CompressedNode, fasta_to_pairlist


class CompressedNodeTests(TestCase):
    """Tests for the _CompressedNode class"""

    def setUp(self):
        """Set up test data for use in compresses node unit tests"""
        self.key = "aba"
        self.values = [1, 2]
        self.node = _CompressedNode(self.key, self.values)

    def test_init(self):
        """Node init should construct the right structure"""
        # With no values should create a node with an empty list for values,
        # the provided key as key, and an empty dictionary as children
        n = _CompressedNode(self.key)
        self.assertEqual(n.values, [])
        self.assertEqual(n.key, self.key)
        self.assertEqual(n.children, {})
        # With values should create a node with the provided values list as
        # values, the provided key as key, and an empty dictionary as children
        n = _CompressedNode(self.key, self.values)
        self.assertEqual(n.values, self.values)
        self.assertEqual(n.key, self.key)
        self.assertEqual(n.children, {})

    def test_truth_value(self):
        """Non zero should check for any data on the node"""
        n = _CompressedNode("")
        self.assertFalse(bool(n))
        self.assertTrue(bool(self.node))

    def test_len(self):
        """Should return the number of values attached to the node"""
        self.assertEqual(len(self.node), 2)

    def test_size(self):
        """Should return the number of nodes attached to the node"""
        self.assertEqual(self.node.size, 1)

    def test_prefix_map(self):
        """Should return the prefix map of the node"""
        exp = {1: [2]}
        self.assertEqual(self.node.prefix_map, exp)

    def test_insert(self):
        """Correctly inserts a new key in the node"""
        n = _CompressedNode(self.key, self.values)
        n.insert("abb", [3])

        # A new node has been create with the common prefix
        self.assertEqual(n.key, "ab")
        self.assertEqual(n.values, [])
        # Tests the old node and the new one has been correctly added
        # as children
        exp_keys = set(["b", "a"])
        self.assertEqual(set(n.children.keys()), exp_keys)
        # Check that the children have the current values
        self.assertEqual(n.children["b"].key, "b")
        self.assertEqual(n.children["b"].values, [[3]])
        self.assertEqual(n.children["b"].children, {})

        self.assertEqual(n.children["a"].key, "a")
        self.assertEqual(n.children["a"].values, [1, 2])
        self.assertEqual(n.children["a"].children, {})

    def test_find(self):
        """The key could be found"""
        # Correctly retrieves the key stored in the calling node
        self.assertEqual(self.node.find("aba"), [1, 2])

        # Correctly retrieves the key stored in a node attached to calling one
        n = _CompressedNode(self.key, self.values)
        n.insert("abb", [3])
        self.assertEqual(n.find("aba"), [1, 2])
        self.assertEqual(n.find("abb"), [[3]])
        self.assertEqual(n.find("ab"), [])

        # Correctly retrieves an empty list for a non existent key
        self.assertEqual(n.find("cd"), [])


class CompressedTrieTests(TestCase):
    """Tests for the CompressedTrie class"""

    def setUp(self):
        """Set up test data for use in compressed trie unit tests"""
        self.data = [("ab",  "0"),
                     ("abababa", "1"),
                     ("abab", "2"),
                     ("baba", "3"),
                     ("ababaa", "4"),
                     ("a", "5"),
                     ("abababa", "6"),
                     ("bab", "7"),
                     ("babba", "8")]
        self.empty_trie = CompressedTrie()
        self.trie = CompressedTrie(self.data)

    def test_init(self):
        """Trie init should construct the right structure"""
        # In no pair_list is provided, it should create an empty Trie
        t = CompressedTrie()
        self.assertEqual(t._root.key, "")
        self.assertEqual(t._root.values, [])
        self.assertEqual(t._root.children, {})
        # If a pair_list is provided, it should insert all the data
        t = CompressedTrie(self.data)
        self.assertEqual(t._root.key, "")
        self.assertEqual(t._root.values, [])
        self.assertEqual(set(t._root.children.keys()), set(["a", "b"]))

    def test_non_zero(self):
        """Non zero should check for any data on the trie"""
        self.assertFalse(self.empty_trie)
        self.assertTrue(self.trie)

    def test_len(self):
        """Should return the number of values attached to the trie"""
        self.assertEqual(len(self.empty_trie), 0)
        self.assertEqual(len(self.trie), 9)

    def test_size(self):
        """Should return the number of nodes attached to the trie"""
        self.assertEqual(self.empty_trie.size, 1)
        self.assertEqual(self.trie.size, 10)

    def test_prefix_map(self):
        """Should map prefix to values"""
        exp = {"1": ["6", "2", "0", "5"],
               "8": ["7"],
               "3": [],
               "4": []}
        self.assertEqual(self.trie.prefix_map, exp)

    def test_insert(self):
        """Correctly inserts a new key into the trie"""
        t = CompressedTrie(self.data)
        t.insert("babc", "9")
        self.assertTrue("9" in t.find("babc"))

        exp = {"1": ["6", "2", "0", "5"],
               "9": ["7"],
               "3": [],
               "4": [],
               "8": []}
        self.assertEqual(t.prefix_map, exp)

    def test_find(self):
        """Correctly founds the values present on the trie"""
        for key, value in self.data:
            self.assertTrue(value in self.trie.find(key))
        self.assertEqual(self.trie.find("cac"), [])
        self.assertEqual(self.trie.find("abababa"), ["1", "6"])


class FastaToPairlistTests(TestCase):
    """Tests for the fasta_to_pairlist function"""

    def setUp(self):
        self.seqs = [("sid_0", "AC"),
                     ("sid_1", "ACAGTC"),
                     ("sid_2", "ACTA"),
                     ("sid_3", "CAGT"),
                     ("sid_4", "CATGAA"),
                     ("sid_5", "A"),
                     ("sid_6", "CATGTA"),
                     ("sid_7", "CAA"),
                     ("sid_8", "CACCA")]

    def test_fasta_to_pairlist(self):
        """Correctly returns a list of (seq, label)"""
        exp = [("AC", "sid_0"),
               ("ACAGTC", "sid_1"),
               ("ACTA", "sid_2"),
               ("CAGT", "sid_3"),
               ("CATGAA", "sid_4"),
               ("A", "sid_5"),
               ("CATGTA", "sid_6"),
               ("CAA", "sid_7"),
               ("CACCA", "sid_8")]

        for obs, exp in izip(fasta_to_pairlist(self.seqs), exp):
            self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
