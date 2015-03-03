r"""
Tree representations (:mod:`skbio.tree`)
========================================

.. currentmodule:: skbio.tree

This module provides functionality for working with trees, including
phylogenetic trees and hierarchies, and prefix trees (i.e., tries).
Functionality is provided for constructing trees, for traversing in multiple
ways, comparisons, fetching subtrees, and more. This module supports trees that
are multifurcating and nodes that have single descendants.

Classes
-------

.. autosummary::
   :toctree: generated/

    TreeNode
    CompressedTrie

Phylogenetic Reconstruction
---------------------------

.. autosummary::
   :toctree: generated/

    nj

Utility Functions
-----------------

.. autosummary::
   :toctree: generated/

    fasta_to_pairlist
    majority_rule

Exceptions
----------

.. autosummary::
   :toctree: generated/

   TreeError
   NoLengthError
   DuplicateNodeError
   MissingNodeError
   NoParentError

Examples
--------

>>> from skbio import TreeNode
>>> from io import StringIO

A new tree can be constructed from a Newick string. Newick is a common format
used to represent tree objects within a file. Newick was part of the original
PHYLIP package from Joseph Felsenstein's group (defined `here
<http://goo.gl/fIY1Iq>`_), and is based around representing nesting with
parentheses. For instance, the following string describes a 3 taxon tree, with
one internal node:

    ((A, B)C, D)root;

Where A, B, and D are tips of the tree, and C is an internal node that covers
tips A and B.

Now let's construct a simple tree and dump an ASCII representation:

>>> tree = TreeNode.read(StringIO(u"((A, B)C, D)root;"))
>>> print tree.is_root()  # is this the root of the tree?
True
>>> print tree.is_tip()  # is this node a tip?
False
>>> print tree.ascii_art()
                    /-A
          /C-------|
-root----|          \-B
         |
          \-D

There are a few common ways to traverse a tree, and depending on your use,
some methods are more appropriate than others. Wikipedia has a well written
page on tree `traversal methods <http://goo.gl/K4Ufl>`_, and will go into
further depth than what we'll cover here. We're only going to cover two of the
commonly used traversals here, preorder and postorder, but we will show
examples of two other common helper traversal methods to gather tips or
internal nodes.

The first traversal we'll cover is a preorder traversal in which you evaluate
from root to tips, looking at the left most child first. For instance:

>>> for node in tree.preorder():
...    print node.name
root
C
A
B
D

The next method we'll look at is a postorder traveral which will evaluate the
left subtree tips first before walking back up the tree:

>>> for node in tree.postorder():
...    print node.name
A
B
C
D
root

`TreeNode` provides two helper methods as well for iterating over just the tips
or for iterating over just the internal nodes.

>>> for node in tree.tips():
...    print "Node name: %s, Is a tip: %s" % (node.name, node.is_tip())
Node name: A, Is a tip: True
Node name: B, Is a tip: True
Node name: D, Is a tip: True

>>> for node in tree.non_tips():
...    print "Node name: %s, Is a tip: %s" % (node.name, node.is_tip())
Node name: C, Is a tip: False

Note, by default, `non_tips` will ignore `self` (which is the root in this
case).  You can pass the `include_self` flag to `non_tips` if you wish to
include `self`.

The `TreeNode` provides a few ways to compare trees. First, let's create two
similar trees and compare their topologies using `compare_subsets`. This
distance is the fraction of common clades present in the two trees, where a
distance of 0 means the trees contain identical clades, and a distance of 1
indicates the trees do not share any common clades:

>>> tree1 = TreeNode.read(StringIO(u"((A, B)C, (D, E)F, (G, H)I)root;"))
>>> tree2 = TreeNode.read(StringIO(u"((G, H)C, (D, E)F, (B, A)I)root;"))
>>> tree3 = TreeNode.read(StringIO(u"((D, B)C, (A, E)F, (G, H)I)root;"))
>>> print tree1.compare_subsets(tree1)  # identity case
0.0
>>> print tree1.compare_subsets(tree2)  # same tree but different clade order
0.0
>>> print tree1.compare_subsets(tree3)  # only 1 of 3 common subsets
0.666666666667

We can additionally take into account branch length when computing distances
between trees. First, we're going to construct two new trees with described
branch length, note the difference in the Newick strings:

>>> tree1 = \
...     TreeNode.read(StringIO(u"((A:0.1, B:0.2)C:0.3, D:0.4, E:0.5)root;"))
>>> tree2 = \
...     TreeNode.read(StringIO(u"((A:0.4, B:0.8)C:0.3, D:0.1, E:0.5)root;"))

In these two trees, we've added on a description of length from the node to
its parent, so for instance:

>>> for node in tree1.postorder():
...     print node.name, node.length
A 0.1
B 0.2
C 0.3
D 0.4
E 0.5
root None

Now let's compare two trees using the distances computed pairwise between tips
in the trees. The distance computed, by default, is the correlation of all
pairwise tip-to-tip distances between trees:

>>> print tree1.compare_tip_distances(tree1)  # identity case
0.0
>>> print tree1.compare_tip_distances(tree2)
0.120492524415

Prefix trees (i.e., tries) examples
-----------------------------------

Construct a Trie from a (key, value) list

>>> from skbio.tree import CompressedTrie
>>> pair_list = [("ab",  "0"),
...              ("abababa", "1"),
...              ("abab", "2"),
...              ("baba", "3"),
...              ("ababaa", "4"),
...              ("a", "5"),
...              ("abababa", "6"),
...              ("bab", "7"),
...              ("babba", "8")]
>>> t = CompressedTrie(pair_list)

Get the number of keys stored in the trie

>>> len(t)
9

Get the number of nodes in the trie

>>> t.size
10

Get the trie's prefix map

>>> t.prefix_map
{'1': ['6', '2', '0', '5'], '8': ['7'], '3': [], '4': []}

Find the value attached to a given key

>>> t.find("ababaa")
['4']

Add a new (key, value) pair to the Trie

>>> t.insert("bac", "9")
>>> t.find("bac")
['9']
>>> t.prefix_map
{'1': ['6', '2', '0', '5'], '9': [], '3': [], '4': [], '8': ['7']}

Create a new trie with a list of sequences

>>> from skbio.tree import fasta_to_pairlist

>>> seqs = [("s0", "ACA"),
...         ("s1", "ACAGTC"),
...         ("s2", "ACTA"),
...         ("s3", "CAGT"),
...         ("s4", "CATGAA"),
...         ("s5", "A"),
...         ("s6", "CATGTA"),
...         ("s7", "CACCA")]

>>> t = CompressedTrie(fasta_to_pairlist(seqs))

>>> t.prefix_map
{'s3': [], 's2': [], 's1': ['s0', 's5'], 's7': [], 's6': [], 's4': []}

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util import TestRunner

from ._tree import TreeNode
from ._trie import CompressedTrie, fasta_to_pairlist
from ._nj import nj
from ._majority_rule import majority_rule
from ._exception import (TreeError, NoLengthError, DuplicateNodeError,
                         MissingNodeError, NoParentError)

__all__ = ['TreeNode', 'CompressedTrie', 'fasta_to_pairlist', 'nj',
           'majority_rule', 'TreeError', 'NoLengthError', 'DuplicateNodeError',
           'MissingNodeError', 'NoParentError']

test = TestRunner(__file__).test
