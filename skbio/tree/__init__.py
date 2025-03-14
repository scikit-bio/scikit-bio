r"""Trees and Phylogenetics (:mod:`skbio.tree`)
===========================================

.. currentmodule:: skbio.tree

This module provides functionality for working with trees, including phylogenetic trees
and various types of hierarchies with relevance in biology or beyond. It supports trees
that are multifurcating and nodes that have single descendants. Functionality is
provided for constructing trees, for traversing in multiple ways, comparisons, fetching
subtrees, and more.

See the :ref:`tree_tutorial` section for working with trees using scikit-bio.


Tree structure and operations
-----------------------------

``TreeNode`` is the data structure for tree representation. A tree consists of an
arbitrary number of interconnected ``TreeNode`` objects representing individual nodes.
Each node has pointers to its parent and children. The class provides a large number of
methods to achieve various tasks for tree analysis and manipulation.

.. autosummary::
   :toctree: generated/
   :template: TreeNode.rst

    TreeNode


Tree Construction
-----------------

Algorithms for phylogenetic reconstruction based on distance matrices.

.. autosummary::
   :toctree: generated/

    upgma
    nj
    gme
    bme


Tree Rearrangement
------------------

Algorithms for improving an existing phylogenetic tree via rearrangement.

.. autosummary::
   :toctree: generated/

    nni


Tree Comparison
---------------

Metrics for assessing the topological or length dissimilarity among trees.

.. autosummary::
   :toctree: generated/

    rf_dists
    wrf_dists
    path_dists


Tree utilities
--------------

.. autosummary::
   :toctree: generated/

    majority_rule


Exceptions
^^^^^^^^^^

.. autosummary::

   TreeError
   NoLengthError
   DuplicateNodeError
   MissingNodeError
   NoParentError


.. _tree_tutorial:

Tutorial
--------

In scikit-bio, the :class:`TreeNode` class is the main and only data structure for tree
represention. It is simple, flexible and powerful.

>>> from skbio import TreeNode

Loading a tree
^^^^^^^^^^^^^^

A tree can be constructed from a Newick string. :wiki:`Newick <Newick_format>` is a
common file format used to represent tree structure based on nesting with parentheses.
For instance, the following string describes a 3-taxon tree, with one internal node:

    ((A,B)C,D)root;

Where A, B, and D are tips of the tree, and C is an internal node that unites tips A
and B. C and D are in turn united by the root of the tree.

If this string is stored in a file ``input.nwk``, we can read it into scikit-bio with:

.. code-block:: python

   tree = TreeNode.read("input.nwk")

Alternatively, we can directly convert a Newick string into a tree with:

>>> tree = TreeNode.read(["((A,B)C,D)root;"])
>>> tree
<TreeNode, name: root, internal node count: 1, tips count: 3>

Let's display an ASCII representation of the tree:

>>> print(tree.ascii_art())
                    /-A
          /C-------|
-root----|          \-B
         |
          \-D

Optionally, we can attach branch lengths to the nodes in a Newick string. In
phylogenetics, the branch length usually represents the amount of evolution
accumulated from one node to another. For example:

>>> tree = TreeNode.read(["((A:2,B:3)C:1,D:4)root;"])

To save a tree to a Newick file:

.. code-block:: python

   tree.write("output.nwk")

Navigating in a tree
^^^^^^^^^^^^^^^^^^^^

The variable ``tree`` is a :class:`TreeNode` object. It is merely the **root** node of
the tree, but one can refer to the entire tree with just the root.

>>> tree.name
'root'

>>> tree.is_root()
True

Starting from the root node one can navigate to its child nodes.

>>> len(tree.children)
2

The first child of root has the name "C". It is an **internal node** of the tree.

>>> node = tree.children[0]
>>> node.name
'C'

>>> node.is_tip()
False

The branch length can be accessed via attribute ``length``. Note, that it is the length
of the branch connecting the current node and its parent.

>>> node.length
1.0

scikit-bio's ``TreeNode`` is highly flexible that we can refer to a subtree (a.k.a., a
clade) within a larger tree using its root node (here "C"). There is no need to detach
the subtree and we can just work with it in-place. Most tree analyses implemented in
scikit-bio also apply to subtrees.

The second child "D" is a **tip** (a.k.a., leaf). In a phylogenetic tree, a tip often
represents a **taxon**, i.e. an organism or a group of organisms being studied.

>>> node = tree.children[1]
>>> node.name
'D'

>>> node.is_tip()
True

Each node has exactly one parent node, with the exception of the root node, whose
parent is None.

>>> node.parent.name
'root'

The ``parent`` and ``children`` attributes connect multiple nodes into a tree
structure. Additionally, one can use methods like :meth:`~TreeNode.siblings`,
:meth:`~TreeNode.neighbors`, :meth:`~TreeNode.ancestors` and :meth:`~TreeNode.root`
to navigate to node with certain relationships to the current node. For example:

>>> node.siblings()[0].name
'C'

One can quickly find a node by name using :meth:`~TreeNode.find`:

>>> node = tree.find("A")
>>> node.name
'A'

If the tree is giant, the first call of ``find`` may take a while, but all subsequent
calls will complete in a flash. That's because ``find`` will cache and reuse a lookup
table under the hood.

Traversing a tree
^^^^^^^^^^^^^^^^^

The process of visiting each and every node within a tree is called :wiki:`tree
traversal <Tree_traversal>`. There are a few common ways to traverse a tree, and
depending on your use, some methods are more appropriate than others. The wiki page
linked above will go into further depth on this topic. Here, we're only going to cover
two of the commonly used traversal methods: :meth:`~TreeNode.preorder` and
:meth:`~TreeNode.postorder`.

A preorder traversal moves from root to tips, looking at the left most child first.
For instance:

>>> for node in tree.preorder():
...    print(node.name)
root
C
A
B
D

A postorder traversal moves from tips to root,  accessing the left subtree tips first
before walking back up the tree:

>>> for node in tree.postorder():
...    print(node.name)
A
B
C
D
root

If you are not certain (or don't mind) what method should be used to traverse the tree,
simply use :meth:`~TreeNode.traverse`.

Note, by default, those traversal methods will include self, i.e., the root node from
which traversal starts. One can pass the ``include_self=False`` flag to the method call
if you wish to exclude it.

``TreeNode`` provides multiple methods for iterating over certain nodes. For examples,
:meth:`~TreeNode.tips` and :meth:`~TreeNode.non_tips` iterate over tips and internal
nodes, respectively:

>>> for node in tree.tips():
...    print(node.name)
A
B
D

>>> for node in tree.non_tips():
...    print(node.name)
C

If we just want to know the total number of nodes (or just tips) in a tree, we can do:

>>> tree.count()
5

>>> tree.count(tips=True)
3

Instead, if we want to obtain a set of taxa (tip names) from a tree, we can use the
:meth:`~TreeNode.subset` method (it is called the "subset" because one can call it at
any node and get a subset of taxa descending from it).

>>> taxa = tree.subset()
>>> sorted(taxa)
['A', 'B', 'D']

Editing a tree
^^^^^^^^^^^^^^

As discussed above, a tree consists of multiple ``TreeNode`` objects connected to each
other. To manipulate a tree, one just needs to change the connections between nodes.
``TreeNode`` provides several standard methods that let you work with nodes like Python
lists.

Create a new tip "E" and attach it to the internal node "C" as a child of it:

>>> node = TreeNode("E")
>>> tree.find("C").append(node)
>>> print(tree.ascii_art())
                    /-A
                   |
          /C-------|--B
         |         |
-root----|          \-E
         |
          \-D

Create two tips "F" and "G" and attach both of them to tip "D". Now "D" becomes an
internal node because it has children:

>>> D = tree.find("D")
>>> D.extend([TreeNode("F"), TreeNode("G")])
>>> print(tree.ascii_art())
                    /-A
                   |
          /C-------|--B
         |         |
-root----|          \-E
         |
         |          /-F
          \D-------|
                    \-G

Prune clade "D" and graft it to node "B". Once a node is attached to a new parent, the
connection to its old parent is automatically removed.

>>> B = tree.find("B")
>>> B.append(D)
>>> print(tree.ascii_art())
                    /-A
                   |
                   |                    /-F
-root---- /C-------|-B------- /D-------|
                   |                    \-G
                   |
                    \-E

Remove clade "D" from node "B". This operation detaches the two nodes, but the removed
node still exists in the memory and once can still use it.

>>> B.remove(D)
True

>>> print(tree.ascii_art())
                    /-A
                   |
-root---- /C-------|--B
                   |
                    \-E

>>> print(D.ascii_art())
          /-F
-D-------|
          \-G

These methods automatically clear lookup tables and other caches because they may be
obsolete for the modified tree. This safety measure could impact performance when
working with very large trees. If caches are not present or relevant, we can add
``uncache=False`` to these methods to skip this check. See :meth:`~TreeNode.has_caches`
for details.

``TreeNode`` provides several advanced methods for editing a tree as a batch. For
examples, :meth:`~TreeNode.bifurcate` converts an arbitrary tree into a bifurcating
tree (binary tree), :meth:`~TreeNode.shear` refines a tree to a given set of taxa,
:meth:`~TreeNode.unpack_by_func` unpacks (contracts) branches that suffice certain
criteria (e.g., bootstrap support below a threshold).

Building a tree
^^^^^^^^^^^^^^^

scikit-bio provides several common algorithms for building a tree structure based on
a distance matrix (:class:`~skbio.stats.distance.DistanceMatrix`) that stores pairwise
distances between taxa.

In the field of :wiki:`phylogenetics <Phylogenetics>`, the evolutionary relationships
among organisms are inferred using genetic data.
:wiki:`Distance-based methods <Distance_matrices_in_phylogeny>` are a common category
of methods for such purpose. The distances are usually computed from sequence alignment
(see :mod:`skbio.alignment`), but they can also be inferred using alternative methods
such as _k_-mer frequency (see :func:`~skbio.sequence.Sequence.kmer_frequencies`) and
community diversity (see :mod:`skbio.diversity.beta`), and the applications of
distance-based tree building are not limited to evolutionary biology, but can be
extended to other fields (e.g., :wiki:`hierarchical clustering
<Hierarchical_clustering>` in data science).

To begin with, let's create a distance matrix of seven taxa.

>>> from skbio import DistanceMatrix
>>> dm = DistanceMatrix([
...     [0.   , 0.005, 0.033, 0.387, 0.626, 0.635, 0.665],
...     [0.005, 0.   , 0.028, 0.387, 0.626, 0.635, 0.665],
...     [0.033, 0.028, 0.   , 0.396, 0.635, 0.64 , 0.679],
...     [0.387, 0.387, 0.396, 0.   , 0.611, 0.611, 0.66 ],
...     [0.626, 0.626, 0.635, 0.611, 0.   , 0.147, 0.758],
...     [0.635, 0.635, 0.64 , 0.611, 0.147, 0.   , 0.754],
...     [0.665, 0.665, 0.679, 0.66 , 0.758, 0.754, 0.   ],
... ], ids=[
...     'human', 'chimp', 'monkey', 'pig', 'mouse', 'rat', 'chicken'
... ])

The following code applies :wiki:`neighbor joining <Neighbor_joining>` (:func:`nj`), a
classic distance-based tree-building algorithm:

>>> from skbio.tree import nj
>>> tree = nj(dm)

Display the reconstructed phylogenetic tree:

>>> print(tree.ascii_art())
          /-human
         |
         |--chimp
---------|
         |          /-monkey
         |         |
          \--------|          /-pig
                   |         |
                    \--------|                    /-mouse
                             |          /--------|
                              \--------|          \-rat
                                       |
                                        \-chicken

Other methods include :func:`upgma` (simplest), :func:`gme` (highly scalable) and
:func:`bme`. All of them can be applied in the same way.

These four methods create a tree from scratch through a fixed set of procedures.
Additionally, starting from an existing tree, :wiki:`tree rearrangement
<Tree_rearrangement>` methods will explore the tree space further in order to
improve it. Nearest neighbor interchange (:func:`nni`) is one such method.

>>> from skbio.tree import nni
>>> tree = nni(tree, dm)

Rooting a tree
^^^^^^^^^^^^^^

The phylogenetic tree created above seems to properly reflect the relatedness among the
seven organisms. However, the "root" of the tree doesn't appear to be aligned with our
knowledge of evolution. This is because the tree is in fact **unrooted**. That is, the
placement of branches does not inform the direction of evolution.

In scikit-bio, the distinguishment between rooted and unrooted trees is implicit: Every
tree has a root node. A tree is considered as "rooted" if its root node has exactly two
children. In contrast, an "unrooted" tree may have three (the most common case), one,
or more than three children attached to its root node.

There are several approaches to root an unrooted tree to reflect the correct direction
of evolution. Midpoint rooting (:meth:`~TreeNode.root_at_midpoint`) is a simple method
that places the root at the midpoint between the two most distant tips in the tree.

>>> rooted = tree.root_at_midpoint()
>>> print(rooted.ascii_art())
          /-chicken
         |
-root----|                    /-mouse
         |          /--------|
         |         |          \-rat
          \--------|
                   |          /-pig
                    \--------|
                             |          /-monkey
                              \--------|
                                       |          /-chimp
                                        \--------|
                                                  \-human

If we already know which taxon or taxa (called an "outgroup") were separated from all
other taxa in evolution, we can perform outgroup rooting
(:meth:`~TreeNode.root_by_outgroup`):

>>> rooted = tree.root_by_outgroup(["chicken"])

Or we can directly root a tree at a given taxon, node or branch:

>>> rooted = tree.root_at("chicken")

To convert a rooted tree into unrooted, we can apply the :meth:`~TreeNode.unroot`
method.

Analyzing topology
^^^^^^^^^^^^^^^^^^

The topology of a tree describes the grouping pattern among taxa. For example, in the
rooted tree displayed above, human, chimp and monkey form a clade (a monophyletic
group), while mouse and rat form another clade. A clade can be referred to by its root
node, which is the lowest common ancestor (:wiki:`LCA <Lowest_common_ancestor>`) of all
descending taxa:

>>> primate = rooted.lca(["human", "chimp", "monkey"])
>>> rodent = rooted.lca(["mouse", "rat"])

Note that one doesn't need to list all descending taxa in order to find the LCA. Two
or more taxa that join at the LCA node on their paths to the root are sufficient for
locating the LCA.

>>> primate = rooted.lca(["human", "monkey"])

A clade regardless of its internal structure can be simply defined by the
:meth:`~TreeNode.subset` of taxa descending from it (see above). All
:meth:`~TreeNode.subsets` within a rooted tree describes all possible taxon groups.

Assume we have another rooted tree ``tree2`` of the same taxa. To assess whether
that tree also support human, chimp and monkey grouped together, we can simply do:

.. code-block:: python

   primate.subset() in tree2.subsets()

When working with unrooted trees, we should replace the notion of subset with
bipartition (:meth:`~TreeNode.bipart`), which defines the separation of all taxa into
two parts by a given branch, regardless of the direction of evolution. For example,
both the rooted and unrooted trees support the separation of primates from other taxa:

>>> primate.bipart() in tree.biparts()
True

The topology of two trees with shared taxa can be compared using the Robinson-Foulds
(RF) distance, which is the number of bipartitions (unrooted trees) or subsets
(rooted trees) that are distinct between the two trees. For example, let's define
three trees, the first two of which have the same topology despite different orders
of clades.

>>> tree1 = TreeNode.read(["((A,B),(C,D),(E,F));"])
>>> tree2 = TreeNode.read(["((F,E),(C,D),(B,A));"])
>>> tree3 = TreeNode.read(["((C,B),(A,D),(E,F));"])

Calculate RF distances using (:meth:`~TreeNode.compare_rfd`):

>>> tree1.compare_rfd(tree1)  # identity case
0.0
>>> tree1.compare_rfd(tree2)  # same tree but different clade order
0.0
>>> tree1.compare_rfd(tree3)  # only 1 of 3 common clade
4.0

We can calculate pairwise RF distances between any number of trees using a single call
of :func:`rf_dists`, which returns a distance matrix:

>>> from skbio.tree import rf_dists
>>> dm = rf_dists([tree1, tree2, tree3])
>>> print(dm)
3x3 distance matrix
IDs:
'0', '1', '2'
Data:
[[ 0.  0.  4.]
 [ 0.  0.  4.]
 [ 4.  4.  0.]]

Analyzing distances
^^^^^^^^^^^^^^^^^^^

With branch lengths annotated, a phylogenetic tree describes the evolutionary distances
among taxa. This important property is the basis for multiple downstream analyses, such
as community diversity (see :mod:`skbio.diversity`).

The :meth:`~TreeNode.distance` between two nodes in a tree (a.k.a., patristic distance)
is the total length of branches constituting the path connecting them. For example:

>>> human, mouse = tree.find("human"), tree.find("mouse")
>>> human.distance(mouse).round(3)
0.628

The maximum distance between any pair of taxa (:meth:`~TreeNode.maxdist`) and the total
branch length connecting all taxa (:meth:`~TreeNode.total_length`) can be used to
measure the biodiversity of groups of organisms.

The distances between all pairs of tips (taxa) in a tree can be efficiently calculated
using :meth:`~TreeNode.cophenet`, which returns a distance matrix.

.. code-block:: python

   dm = tree.cophenet()

The evolutionary distances reflected by two trees can be compared. First, let's create
two trees with the same taxa but different branch lengths.

>>> tree1 = TreeNode.read(["((A:0.1,B:0.2):0.3,C:0.4,D:0.5);"])
>>> tree2 = TreeNode.read(["((A:0.4,B:0.8):0.3,C:0.1,D:0.5);"])

Calculate the distance between the cophenetic distance matrics using
(:meth:`~TreeNode.compare_cophenet`):

>>> tree1.compare_cophenet(tree2).round(3)
0.12

Likewise, we can calculate this metric between multiple trees using :func:`path_dists`.


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._tree import TreeNode
from ._upgma import upgma
from ._nj import nj
from ._me import gme, bme, nni
from ._majority_rule import majority_rule
from ._compare import rf_dists, wrf_dists, path_dists
from ._exception import (
    TreeError,
    NoLengthError,
    DuplicateNodeError,
    MissingNodeError,
    NoParentError,
)

__all__ = [
    "TreeNode",
    "upgma",
    "nj",
    "gme",
    "bme",
    "nni",
    "rf_dists",
    "wrf_dists",
    "path_dists",
    "majority_rule",
    "TreeError",
    "NoLengthError",
    "DuplicateNodeError",
    "MissingNodeError",
    "NoParentError",
]
