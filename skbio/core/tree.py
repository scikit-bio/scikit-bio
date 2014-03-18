#!/usr/bin/env python

r"""
Tree representations (:mod:`skbio.core.tree`)
=============================================

.. currentmodule:: skbio.core.tree

This module provides functionality for working with trees, including
phylogenetic trees and hierarchies. Functionality is provided for
constructing the trees, for traversing in multiple ways, comparisons,
fetching subtrees, and more. This module supports trees that are multifurcating
and support nodes that have single descendents as well.

Classes
-------

.. autosummary::
   :toctree: generated/

   TreeNode

Examples
--------

>>> from skbio.core.tree import TreeNode

A new tree can be constructed from a Newick string. Newick is a common format
used to represent tree objects within a file. Newick was part of the original
PHYLIP package from Joseph Felsenstein's group (defined
`here <http://goo.gl/fIY1Iq>`_), and is based around representing nesting with
parentheses. For instance, the following string describes a 3 taxon tree, with
one internal node:

    ((A, B)C, D)root;

Where A, B, and D are tips of the tree, and C is an internal node that covers
tips A and B.

Now lets construct a simple tree and dump an ASCII representation:

>>> tree = TreeNode.from_newick("((A, B)C, D)root;")
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
commonly used traversals here, preorder and postorder, but we well show
examples of two other common helper traversal methods to gather tips or
internal nodes.

The first traversal we'll cover is a preorder traversal in which you evaluate
from root to tips, looking at the left most child first. For instance:

>>> for node in tree.preorder():
...    print node.Name
root
C
A
B
D

The next method we'll look at is a postorder traveral which will evaluate the
left subtree tips first before walking back up the tree:

>>> for node in tree.postorder():
...    print node.Name
A
B
C
D
root

TreeNode provides two helper methods as well for iterating over just the tips
or for iterating over just the internal nodes.

>>> for node in tree.tips():
...    print "Node name: %s, Is a tip: %s" % (node.Name, node.is_tip())
Node name: A, Is a tip: True
Node name: B, Is a tip: True
Node name: D, Is a tip: True

>>> for node in tree.non_tips():
...    print "Node name: %s, Is a tip: %s" % (node.Name, node.is_tip())
Node name: C, Is a tip: False

Note, by default, non_tips will ignore self (which is the root in this case).
You can pass the include_self flag to non_tips if you wish to include self.

The TreeNode provides a few ways to compare trees. First, let's create two
similar trees and compare their topologies using compare_subsets. This
distance is the fraction of common clades present in the two trees, where a
distance of 0 means the trees contain identical clades, and a distance of 1
indicates the trees do not share any common clades:

>>> tree1 = TreeNode.from_newick("((A, B)C, (D, E)F, (G, H)I)root;")
>>> tree2 = TreeNode.from_newick("((G, H)C, (D, E)F, (B, A)I)root;")
>>> tree3 = TreeNode.from_newick("((D, B)C, (A, E)F, (G, H)I)root;")
>>> print tree1.compare_subsets(tree1)  # identity case
0.0
>>> print tree1.compare_subsets(tree2)  # same tree but different clade order
0.0
>>> print tree1.compare_subsets(tree3)  # only 1 of 3 common subsets
0.666666666667

We can additionally take into account branch length when computing distances
between trees. First, we're going to construct two new trees with described
branch length, note the difference in the Newick strings:

>>> tree1 = TreeNode.from_newick("((A:0.1, B:0.2)C:0.3, D:0.4, E:0.5)root;")
>>> tree2 = TreeNode.from_newick("((A:0.4, B:0.8)C:0.3, D:0.1, E:0.5)root;")

In these two trees, we've added on a description of length from the node to
its parent, so for instance:

>>> for node in tree1.postorder():
...     print node.Name, node.Length
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
"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import re
from operator import or_
from random import shuffle
from copy import deepcopy
from itertools import combinations

import numpy as np

from skbio.maths.stats.test import correlation_t
from skbio.core.exception import (NoLengthError, DuplicateNodeError,
                                  NoParentError, MissingNodeError,
                                  TreeError, RecordError)


def _dnd_tokenizer(data):
    """Tokenizes data into a stream of punctuation, labels and lengths.

    Parameters
    ----------
    data : str
        a DND-like (e.g., newick) string

    Returns
    -------
    GeneratorType
        Yields successive DND tokens

    See also
    --------
    TreeNode.from_newick
    TreeNode.to_newick

    Examples
    --------
    >>> from skbio.core.tree import _dnd_tokenizer
    >>> for token in _dnd_tokenizer("((tip1, tip2)internal1)"):
    ...     print token
    (
    (
    tip1
    ,
    tip2
    )
    internal1
    )
    """
    dnd_token_str = '(:),;'
    dnd_tokens = set(dnd_token_str)

    in_quotes = False
    saved = []
    sa = saved.append
    for d in data:
        if d == "'":
            in_quotes = not in_quotes
        if d in dnd_tokens and not in_quotes:
            curr = ''.join(saved).strip()
            if curr:
                yield curr
            yield d
            saved = []
            sa = saved.append
        else:
            sa(d)


def distance_from_r(m1, m2):
    """Estimates distance as (1-r)/2: neg correl = max distance

    Parameters
    ----------
    m1 : SymmetricDistanceMatrix
        a distance matrix to compare
    m2 : SymmetricDistanceMatrix
        a distance matrix to compare

    Returns
    -------

    float
        The distance between m1 and m2
    """
    return (1-correlation_t(m1.flat, m2.flat)[0])/2


class TreeNode(object):
    """Representation of a node within a tree

    A `TreeNode` instance stores links to its parent and optional children
    nodes. In addition, the `TreeNode` can represent a `Length` (e.g., a
    branch length) between itself and its parent.

    Parameters
    ----------
    name : str or None
        A node can have a name. It is common for tips in particular to have
        names, for instance, in a phylogenetic tree where the tips correspond
        to species.
    length : float, int, or None
        Distances between nodes can be used to represent evolutionary
        distances, time, etc.
    parent : TreeNode or None
        Connect this node to a parent
    children : list of TreeNode or None
        Connect this node to existing children

    Attributes
    ----------
    name
    length
    parent
    children
    id

    """

    _exclude_from_copy = set(['Parent', 'Children', '_node_cache'])

    def __init__(self, Name=None, Length=None, Parent=None, Children=None):
        self.Name = Name
        self.Length = Length
        self.Parent = Parent
        self._node_cache = {}
        self.Children = []
        self.Id = None

        if Children:
            for c in Children:
                self.append(c)

    ### start operators ###
    def __repr__(self):
        """Returns summary of the tree

        Returns
        -------
        str
            A summary of this node and all descendents

        Notes
        -----
        This method returns the name of the node and a count of tips and the
        number of internal nodes in the tree

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c, d)root;")
        >>> repr(tree)
        '<TreeNode, name: root internal node count: 1, tips count: 3>'

        """
        nodes = [n for n in self.traverse(include_self=False)]
        n_tips = sum([n.is_tip() for n in nodes])
        n_nontips = len(nodes) - n_tips
        classname = self.__class__.__name__
        name = self.Name if self.Name is not None else "unnamed"

        return "<%s, name: %s internal node count: %d, tips count: %d>" % \
               (classname, name, n_nontips, n_tips)

    def __str__(self):
        """Returns string version of self, with names and distances

        Returns
        -------
        str
            Returns a Newick representation of the tree

        See Also
        --------
        TreeNode.to_newick
        TreeNode.from_newick

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c);")
        >>> str(tree)
        '((a,b)c);'

        """

        return self.to_newick(with_distances=True)

    def __iter__(self):
        """Node iter iterates over the Children."""
        return iter(self.Children)

    def __len__(self):
        return len(self.Children)

    def __getitem__(self, i):
        """Node delegates slicing to Children"""
        return self.Children[i]

    ### end operators ###

    ### start topology updates ###
    def _adopt(self, node):
        """Update parent references but does NOT update self.Children"""
        self.invalidate_node_cache()
        if node.Parent is not None:
            node.Parent.remove(node)
        node.Parent = self
        return node

    def append(self, node):
        """Appends a node to self.Children, in-place, cleaning up refs

        `append` will invalidate any node lookup caches, remove an existing
        parent on `node` if one exists, set the parent of `node` to `self`
        and add the `node` to `self`s `children`.

        Parameters
        ----------
        node : TreeNode
            An existing TreeNode object

        See Also
        --------
        TreeNode.extend

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> root = TreeNode(Name="root")
        >>> child1 = TreeNode(Name="child1")
        >>> child2 = TreeNode(Name="child2")
        >>> root.append(child1)
        >>> root.append(child2)
        >>> print root
        (child1,child2)root;

        """
        self.Children.append(self._adopt(node))

    def extend(self, nodes):
        """Append a list of nodes to self

        `extend` will invalidate any node lookup caches, remoev existing
        parents of the `nodes` if they have any, set their parents to `self
        and add the nodes to `self`s `children`.

        Parameters
        ----------
        nodes : list of TreeNode
            A list of TreeNode objects

        See Also
        --------
        TreeNode.append

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> root = TreeNode(Name="root")
        >>> root.extend([TreeNode(Name="child1"), TreeNode(Name="child2")])
        >>> print root
        (child1,child2)root;

        """
        self.Children.extend([self._adopt(n) for n in nodes])

    def pop(self, index=-1):
        """Remove a node from self

        Remove a child node by its index position. All node lookup caches
        are invalidated, and the parent reference for the popped node will be
        set to None.

        Parameters
        ----------
        index : int
            The index position in children to pop

        Returns
        -------
        TreeNode
            The popped child

        See Also
        --------
        TreeNode.remove
        TreeNode.remove_deleted

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(a,b)c;")
        >>> print tree.pop(0)
        a;

        """
        return self._remove_node(index)

    def _remove_node(self, idx):
        """The actual (and only) method that performs node removal"""
        self.invalidate_node_cache()
        node = self.Children.pop(idx)
        node.Parent = None
        return node

    def remove(self, node):
        """Remove a node from self

        Remove a `node` from `self` by identity of the node.

        Parameters
        ----------
        node : TreeNode
            The node to remove from self's children

        Returns
        -------
        bool
            True if the node was removed, False otherwise

        See Also
        --------
        TreeNode.pop
        TreeNode.remove_deleted

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(a,b)c;")
        >>> tree.remove(tree.Children[0])
        True

        """
        for (i, curr_node) in enumerate(self.Children):
            if curr_node == node:
                self._remove_node(i)
                return True
        return False

    def remove_deleted(self, func):
        """Delete nodes in which func(node) evaluates True

        Remove all descendents from self that evaluate True from `func`. This
        has the potential to drop clades.

        Parameters
        ----------
        func : a function
            A function that evaluates `True` when a node should be deleted

        See Also
        --------
        TreeNode.pop
        TreeNode.remove

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(a,b)c;")
        >>> tree.remove_deleted(lambda x: x.Name == 'b')
        >>> print tree
        (a)c;
        """
        for node in self.traverse(include_self=False):
            if func(node):
                node.Parent.remove(node)

    def prune(self):
        """Reconstructs correct topology after nodes have been removed.

        Internal nodes with only one child will be removed and new connections
        will be made to reflect change. This method is useful to call
        following node removals as it will clean up nodes with singular
        children.

        Names and properties of singular children will override the names and
        properties of their parents following the prune.

        Node lookup caches are invalidated.

        See Also
        --------
        TreeNode.shear
        TreeNode.remove
        TreeNode.pop
        Treenode.remove_deleted

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c,(d,e)f)root;")
        >>> to_delete = tree.find('b')
        >>> tree.remove_deleted(lambda x: x == to_delete)
        >>> print tree
        ((a)c,(d,e)f)root;
        >>> tree.prune()
        >>> print tree
        ((d,e)f,a)root;
        """
        # build up the list of nodes to remove so the topology is not altered
        # while traversing
        nodes_to_remove = []
        for node in self.traverse(include_self=False):
            if len(node.Children) == 1:
                nodes_to_remove.append(node)

        # clean up the single children nodes
        for node in nodes_to_remove:
            node.Parent.append(node.Children[0])
            node.Parent.remove(node)

#   def shear(self, names):
#       """Lop off tips until the tree just has the desired tip names"""
#       tcopy = self.deepcopy()
#       all_tips = set([n.Name for n in tcopy.tips()])
#       ids = set(names)
#
#       if not ids.issubset(all_tips):
#           raise ValueError("ids are not a subset of the tree!")
#
#       while len(tcopy.tips()) != len(ids):
#           for n in tcopy.tips():
#               if n.Name not in ids:
#                   n.Parent.removeNode(n)
#
#       tcopy.prune()
#       return tcopy
#
    ### end topology updates ###

    ### copy like methods ###
    def copy(self):
        """Returns a copy of self using an iterative approach

        Perform an iterative deepcopy of self. It is not assured that the copy
        of node attributes will be performed iteratively as that depends on
        the copy method of the types being copied

        Returns
        -------
        TreeNode
            A new copy of self

        See Also
        --------
        TreeNode.unrooted_deepcopy
        TreeNode.unrooted_copy

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c,(d,e)f)root;")
        >>> tree_copy = tree.copy()
        >>> tree_nodes = set([id(n) for n in tree.traverse()])
        >>> tree_copy_nodes = set([id(n) for n in tree_copy.traverse()])
        >>> print len(tree_nodes.intersection(tree_copy_nodes))
        0
        """
        def __copy_node(node_to_copy):
            """Helper method to copy a node"""
            # this is _possibly_ dangerous, we're assuming the node to copy is
            # of the same class as self, and has the same exclusion criteria.
            # however, it is potentially dangerous to mix TreeNode subclasses
            # within a tree, so...
            result = self.__class__()
            efc = self._exclude_from_copy
            for key in node_to_copy.__dict__:
                if key not in efc:
                    result.__dict__[key] = deepcopy(node_to_copy.__dict__[key])
            return result

        root = __copy_node(self)
        nodes_stack = [[root, self, len(self.Children)]]

        while nodes_stack:
            #check the top node, any children left unvisited?
            top = nodes_stack[-1]
            new_top_node, old_top_node, unvisited_children = top

            if unvisited_children:
                top[2] -= 1
                old_child = old_top_node.Children[-unvisited_children]
                new_child = __copy_node(old_child)
                new_top_node.append(new_child)
                nodes_stack.append([new_child, old_child,
                                    len(old_child.Children)])
            else:  # no unvisited children
                nodes_stack.pop()
        return root

    __copy__ = copy
    __deepcopy__ = deepcopy = copy

    def unrooted_deepcopy(self, parent=None):
        """Walks the tree unrooted-style and returns a new copy

        Perform a deepcopy of self and return a new copy of the tree as an
        unrooted copy. This is useful for defining new roots of the tree as
        the `TreeNode`.

        This method calls TreeNode.unrooted_copy which is recursive.

        Parameters
        ----------
        parent : TreeNode or None
            Used to avoid infinite loops when performing the unrooted traverse

        Returns
        -------
        TreeNode
            A new copy of the tree

        See Also
        --------
        TreeNode.copy
        TreeNode.unrooted_copy
        TreeNode.root_at

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,(b,c)d)e,(f,g)h)i;")
        >>> new_tree = tree.find('d').unrooted_deepcopy()
        >>> print new_tree
        (b,c,(a,((f,g)h)e)d)root;
        """
        root = self.root()
        root.assign_ids()

        new_tree = root.copy()
        new_tree.assign_ids()

        new_tree_self = new_tree.find_by_id(self.Id)
        return new_tree_self.unrooted_copy(parent)

    def unrooted_copy(self, parent=None):
        """Walks the tree unrooted-style and returns a copy

        Perform a copy of self and return a new copy of the tree as an
        unrooted copy. This is useful for defining new roots of the tree as
        the `TreeNode`.

        This method is recursive.

        Warning, this is _NOT_ a deepcopy

        Parameters
        ----------
        parent : TreeNode or None
            Used to avoid infinite loops when performing the unrooted traverse

        Returns
        -------
        TreeNode
            A new copy of the tree

        See Also
        --------
        TreeNode.copy
        TreeNode.unrooted_deepcopy
        TreeNode.root_at

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,(b,c)d)e,(f,g)h)i;")
        >>> new_tree = tree.find('d').unrooted_copy()
        >>> print new_tree
        (b,c,(a,((f,g)h)e)d)root;
        """
        neighbors = self.neighbors(ignore=parent)
        children = [c.unrooted_copy(parent=self) for c in neighbors]

        # we might be walking UP the tree, so:
        if parent is None:
            # base edge
            edgename = None
            length = None
        elif parent.Parent is self:
            # self's parent is becoming self's child
            edgename = parent.Name
            length = parent.Length
        else:
            if parent is not self.Parent:
                raise TreeError("Should never happen...")
            edgename = self.Name
            length = self.Length

        result = self.__class__(Name=edgename, Children=children,
                                Length=length)

        if parent is None:
            result.Name = "root"

        return result

    def subtree(self, tip_list=None):
        """Make a copy of the subtree"""
        raise NotImplementedError()

    def subset(self):
        """Returns set of names that descend from specified node

        Get the set of Names on tips that descend from this node

        Returns
        -------
        frozenset
            The set of names at the tips of the clade that descends from self

        See Also
        --------
        TreeNode.subsets
        TreeNode.compare_subsets

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,(b,c)d)e,(f,g)h)i;")
        >>> sorted(tree.subset())
        ['a', 'b', 'c', 'f', 'g']
        """
        return frozenset({i.Name for i in self.tips()})

    def subsets(self):
        """Return all sets of names that come from self and its descendents

        Compute all subsets of tip names over self, or, represent a tree as a
        set of nested sets.

        Returns
        -------
        frozenset
            A frozenset of frozensets of str

        See Also
        --------
        TreeNode.subset
        TreeNode.compare_subsets

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(((a,b)c,(d,e)f)h)root;")
        >>> for s in sorted(tree.subsets()):
        ...     print sorted(s)
        ['a', 'b']
        ['d', 'e']
        ['a', 'b', 'd', 'e']
        """
        sets = []
        for i in self.postorder(include_self=False):
            if not i.Children:
                i.__leaf_set = frozenset([i.Name])
            else:
                leaf_set = reduce(or_, [c.__leaf_set for c in i.Children])
                if len(leaf_set) > 1:
                    sets.append(leaf_set)
                i.__leaf_set = leaf_set
        return frozenset(sets)

    def root_at(self, node):
        """Return a new tree rooted at the provided node.

        This can be useful for drawing unrooted trees with an orientation that
        reflects knowledge of the true root location.

        Parameters
        ----------
        node : TreeNode or str
            The node to root at

        Returns
        -------
        TreeNode
            A new copy of the tree

        Raises
        ------
        TreeError
            Raises a `TreeError` if a tip is specified as the new root

        See Also
        --------
        TreeNode.root_at_midpoint
        TreeNode.unrooted_deepcopy

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(((a,b)c,(d,e)f)g,h)i;")
        >>> print tree.root_at('c')
        (a,b,((d,e)f,(h)g)c)root;

        """
        if isinstance(node, str):
            node = self.find(node)

        if not node.Children:
            raise TreeError("Can't use a tip (%s) as the root" %
                            repr(node.Name))
        return node.unrooted_deepcopy()

    def root_at_midpoint(self):
        """Return a new tree rooted at midpoint of the two tips farthest apart

        This method doesn't preserve the internal node naming or structure,
        but does keep tip to tip distances correct. Uses unrooted_copy() but
        operates on a full copy of the tree.

        Raises
        ------
        TreeError
            If a tip ends up being the mid point

        Returns
        -------
        TreeNode
            A tree rooted at its midpoint
        LengthError
            Midpoint rooting requires `Length` and will raise (indirectly) if
            evaluated nodes don't have `Length`

        See Also
        --------
        TreeNode.root_at
        TreeNode.unrooted_deepcopy

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("(((d:1,e:1,(g:1)f:1)c:1)b:1,h:1)a:1;")
        >>> print tree.root_at_midpoint()
        ((d:1.0,e:1.0,(g:1.0)f:1.0)c:0.5,((h:1.0)b:1.0):0.5)root;

        """
        tree = self.copy()
        max_dist, tips = tree.get_max_distance()
        half_max_dist = max_dist/2.0

        if max_dist == 0.0:  # only pathological cases with no lengths
            return tree

        tip1 = tree.find(tips[0])
        tip2 = tree.find(tips[1])
        lca = tree.lowest_common_ancestor([tip1, tip2])

        if tip1.accumulate_to_ancestor(lca) > half_max_dist:
            climb_node = tip1
        else:
            climb_node = tip2

        dist_climbed = 0.0
        while dist_climbed + climb_node.Length < half_max_dist:
            dist_climbed += climb_node.Length
            climb_node = climb_node.Parent

        # now midpt is either at on the branch to climb_node's  parent
        # or midpt is at climb_node's parent
        if dist_climbed + climb_node.Length == half_max_dist:
            # climb to midpoint spot
            climb_node = climb_node.Parent
            if climb_node.is_tip():
                raise TreeError('error trying to root tree at tip')
            else:
                return climb_node.unrooted_copy()

        else:
            # make a new node on climb_node's branch to its parent
            old_br_len = climb_node.Length

            new_root = tree.__class__()
            climb_node.Parent.append(new_root)
            new_root.append(climb_node)

            climb_node.Length = half_max_dist - dist_climbed
            new_root.Length = old_br_len - climb_node.Length

            return new_root.unrooted_copy()

    ### end copy like methods ###

    ### node checks ###

    def is_tip(self):
        """Returns True if the current node is a tip, i.e. has no children.

        Returns
        -------
        bool
            `True` if the node is a tip

        See Also
        --------
        TreeNode.is_root
        TreeNode.has_children

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c);")
        >>> print tree.is_tip()
        False
        >>> print tree.find('a').is_tip()
        True
        """
        return not self.Children

    def is_root(self):
        """Returns True if the current is a root, i.e. has no parent.

        Returns
        -------
        bool
            `True` if the node is the root

        See Also
        --------
        TreeNode.is_tip
        TreeNode.has_children

        Examples
        --------
        >>> from skbio.core.tree import TreeNode
        >>> tree = TreeNode.from_newick("((a,b)c);")
        >>> print tree.is_root()
        True
        >>> print tree.find('a').is_root()
        False
        """
        return self.Parent is None

    def has_children(self):
        """Returns True if self.Children."""
        return bool(self.Children)

    ### end node checks ###

    ### traversal methods ###

    def traverse(self, self_before=True, self_after=False, include_self=True):
        """Returns iterator over descendants. Iterative: safe for large trees.

        self_before includes each node before its descendants if True.
        self_after includes each node after its descendants if True.
        include_self includes the initial node if True.

        self_before and self_after are independent. If neither is True, only
        terminal nodes will be returned.

        Note that if self is terminal, it will only be included once even if
        self_before and self_after are both True.

        This is a depth-first traversal. Since the trees are not binary,
        preorder and postorder traversals are possible, but inorder traversals
        would depend on the data in the tree and are not handled here.
        """
        if self_before:
            if self_after:
                return self.pre_and_postorder(include_self=include_self)
            else:
                return self.preorder(include_self=include_self)
        else:
            if self_after:
                return self.postorder(include_self=include_self)
            else:
                return self.tips(include_self=include_self)

    def preorder(self, include_self=True):
        """Performs preorder iteration over tree."""
        stack = [self]
        while stack:
            curr = stack.pop()
            if include_self or (curr is not self):
                yield curr
            if curr.Children:
                stack.extend(curr.Children[::-1])

    def postorder(self, include_self=True):
        """Performs postorder iteration over tree.

        This is somewhat inelegant compared to saving the node and its index
        on the stack, but is 30% faster in the average case and 3x faster in
        the worst case (for a comb tree).

        Zongzhi Liu's slower but more compact version is:

        def postorder_zongzhi(self):
            stack = [[self, 0]]
            while stack:
                curr, child_idx = stack[-1]
                if child_idx < len(curr.Children):
                    stack[-1][1] += 1
                    stack.append([curr.Children[child_idx], 0])
                else:
                    yield stack.pop()[0]
        """
        child_index_stack = [0]
        curr = self
        curr_children = self.Children
        curr_children_len = len(curr_children)
        while 1:
            curr_index = child_index_stack[-1]
            #if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                #if the current child has children, go there
                if curr_child.Children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.Children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                #otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            #if there are no children left, return self, and move to
            #self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.Parent
                curr_children = curr.Children
                curr_children_len = len(curr_children)
                child_index_stack.pop()
                child_index_stack[-1] += 1

    def pre_and_postorder(self, include_self=True):
        """Performs iteration over tree, visiting node before and after."""
        #handle simple case first
        if not self.Children:
            if include_self:
                yield self
            raise StopIteration
        child_index_stack = [0]
        curr = self
        curr_children = self.Children
        while 1:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            #if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                #if the current child has children, go there
                if curr_child.Children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.Children
                    curr_index = 0
                #otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            #if there are no children left, return self, and move to
            #self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.Parent
                curr_children = curr.Children
                child_index_stack.pop()
                child_index_stack[-1] += 1

    def levelorder(self, include_self=True):
        """Performs levelorder iteration over tree"""
        queue = [self]
        while queue:
            curr = queue.pop(0)
            if include_self or (curr is not self):
                yield curr
            if curr.Children:
                queue.extend(curr.Children)

    def tips(self, include_self=False):
        """Iterates over tips descended from self, [] if self is a tip."""
        #bail out in easy case
        if not self.Children:
            if include_self:
                yield self
            raise StopIteration

        stack = [self]
        while stack:
            curr = stack.pop()
            if curr.Children:
                stack.extend(curr.Children[::-1])  # 20% faster than reversed
            else:
                yield curr

    def non_tips(self, include_self=False):
        """Iterates over nontips descended from self, [] if none.

        include_self, if True (default is False), will return the current
        node as part of the list of nontips if it is a nontip.
        """
        for n in self.postorder(include_self):
            if n.Children:
                yield n

    ### end traversal methods ###

    ### search methods ###

    def invalidate_node_cache(self):
        """Delete the node cache"""
        self._node_cache = {}

    def create_node_cache(self):
        """Construct an internal lookup keyed by node name, valued by node

        This method will not cache nodes in which the .Name is None. This
        method will raise DuplicateNodeError if a name conflict is discovered.
        """
        if self._node_cache:
            return

        for node in self.traverse():
            name = node.Name
            if name is None:
                continue

            if name in self._node_cache:
                raise DuplicateNodeError("%s already exists!" % name)

            self._node_cache[name] = node

    def find(self, name):
        """Find a node by name

        This method returns raises MissingNodeError if the node is not found.
        The first time this method is called, an internal cache is
        constructed to improve performance on subsequent calls.
        """
        # if what is being passed in looks like a node, just return it
        if isinstance(name, self.__class__):
            return name

        self.create_node_cache()
        node = self._node_cache.get(name, None)

        if node is None:
            raise MissingNodeError("Node %s is not in self" % name)
        else:
            return node

    def find_by_id(self, id_):
        """Find a node by id

        This method returns raises MissingNodeError if the node is not found.
        The first time this method is called, an internal cache is
        constructed to improve performance on subsequent calls.
        """
        # if this method gets used frequently, then we should cache by ID
        # as well
        root = self.root()
        root.assign_ids()

        node = None
        for n in self.traverse(include_self=True):
            if n.Id == id_:
                node = n

        if node is None:
            raise MissingNodeError("ID %d is not in self" % id_)
        else:
            return node

    ### path methods ###

    def ancestors(self):
        """Returns all ancestors back to the root. Dynamically calculated."""
        if self.is_root():
            return []

        result = []
        curr = self.Parent
        while not curr.is_root():
            result.append(curr)
            curr = curr.Parent
        result.append(curr)

        return result

    def root(self):
        """Returns root of the tree self is in. Dynamically calculated."""
        if self.is_root():
            return self

        curr = self
        while not curr.is_root():
            curr = curr.Parent
        return curr

    def siblings(self):
        """Returns all nodes that are children of the same parent as self.

        Note: excludes self from the list. Dynamically calculated.
        """
        if self.is_root():
            return []

        result = self.Parent.Children[:]
        result.remove(self)

        return result

    def neighbors(self, ignore=None):
        """Get neighbors of self"""
        nodes = [n for n in self.Children + [self.Parent] if n is not None]
        if ignore is None:
            return nodes
        else:
            return [n for n in nodes if n is not ignore]

    def lowest_common_ancestor(self, tipnames):
        """Lowest common ancestor for a list of tipnames

        This should be around O(H sqrt(n)), where H is height and n is the
        number of tips passed in.
        """
        if len(tipnames) == 1:
            return self.find(tipnames[0])

        tips = [self.find(name) for name in tipnames]

        if len(tips) == 0:
            return None

        nodes_to_scrub = []

        for t in tips:
            if t.is_root():
                # has to be the LCA...
                return t

            prev = t
            curr = t.Parent

            while curr and not hasattr(curr, 'black'):
                setattr(curr, 'black', [prev])
                nodes_to_scrub.append(curr)
                prev = curr
                curr = curr.Parent

            # increase black count, multiple children lead to here
            if curr:
                curr.black.append(prev)

        curr = self
        while len(curr.black) == 1:
            curr = curr.black[0]

        # clean up tree
        for n in nodes_to_scrub:
            delattr(n, 'black')

        return curr

    lca = lowest_common_ancestor  # for convenience

    ### end path methods ###

    ### parsers ###

    @classmethod
    def from_newick(cls, lines, unescape_name=True):
        """Returns tree from the Clustal .dnd file format and equivalent.

        Tree is made of skbio.core.tree.TreeNode objects, with branch lengths
        """
        def _new_child(old_node):
            """Returns new_node which has old_node as its parent."""
            new_node = cls()
            new_node.Parent = old_node
            if old_node is not None:
                if new_node not in old_node.Children:
                    old_node.Children.append(new_node)
            return new_node

        if isinstance(lines, str):
            data = lines
        else:
            data = ''.join(lines)

        #skip arb comment stuff if present: start at first paren
        paren_index = data.find('(')
        data = data[paren_index:]
        left_count = data.count('(')
        right_count = data.count(')')

        if left_count != right_count:
            raise RecordError("Found %s left parens but %s right parens." %
                              (left_count, right_count))

        curr_node = None
        state = 'PreColon'
        state1 = 'PreClosed'
        last_token = None

        for t in _dnd_tokenizer(data):
            if t == ':':
                #expecting branch length
                state = 'PostColon'
                #prevent state reset
                last_token = t
                continue
            if t == ')' and last_token in ',(':
                # node without name
                new_node = _new_child(curr_node)
                new_node.Name = None
                curr_node = new_node.Parent
                state1 = 'PostClosed'
                last_token = t
                continue
            if t == ')':
                #closing the current node
                curr_node = curr_node.Parent
                state1 = 'PostClosed'
                last_token = t
                continue
            if t == '(':
                #opening a new node
                curr_node = _new_child(curr_node)
            elif t == ';':  # end of data
                last_token = t
                break
            elif t == ',' and last_token in ',(':
                # node without name
                new_node = _new_child(curr_node)
                new_node.Name = None
                curr_node = new_node.Parent
            elif t == ',':
                # separator: next node adds to this node's parent
                curr_node = curr_node.Parent
            elif state == 'PreColon' and state1 == 'PreClosed':
                # data for the current node
                new_node = _new_child(curr_node)
                if unescape_name:
                    if t.startswith("'") and t.endswith("'"):
                        while t.startswith("'") and t.endswith("'"):
                            t = t[1:-1]
                    else:
                        if '_' in t:
                            t = t.replace('_', ' ')
                new_node.Name = t
                curr_node = new_node
            elif state == 'PreColon' and state1 == 'PostClosed':
                if unescape_name:
                    while t.startswith("'") and t.endswith("'"):
                        t = t[1:-1]
                curr_node.Name = t
            elif state == 'PostColon':
                # length data for the current node
                curr_node.Length = float(t)
            else:
                # can't think of a reason to get here
                raise RecordError("Incorrect PhyloNode state? %s" % t)
            state = 'PreColon'  # get here for any non-colon token
            state1 = 'PreClosed'
            last_token = t

        if curr_node is not None and curr_node.Parent is not None:
            raise RecordError("Didn't get back to root of tree.")

        if curr_node is None:  # no data -- return empty node
            return cls()
        return curr_node  # this should be the root of the tree

    ### end parsers ###

    ### formatters ###

    def to_newick(self, with_distances=False, semicolon=True,
                  escape_name=True):
        """Return the newick string for this tree.

        Arguments:
            - with_distances: whether branch lengths are included.
            - semicolon: end tree string with a semicolon
            - escape_name: if any of these characters []'"(),:;_ exist in a
                nodes name, wrap the name in single quotes

        NOTE: This method returns the Newick representation of this node
        and its descendents. This method is a modification of an implementation
        by Zongzhi Liu
        """
        result = ['(']
        nodes_stack = [[self, len(self.Children)]]
        node_count = 1

        while nodes_stack:
            node_count += 1
            #check the top node, any children left unvisited?
            top = nodes_stack[-1]
            top_node, num_unvisited_children = top
            if num_unvisited_children:  # has any child unvisited
                top[1] -= 1  # decrease the #of children unvisited
                next_child = top_node.Children[-num_unvisited_children]
                # pre-visit
                if next_child.Children:
                    result.append('(')
                nodes_stack.append([next_child, len(next_child.Children)])
            else:  # no unvisited children
                nodes_stack.pop()
                #post-visit
                if top_node.Children:
                    result[-1] = ')'

                if top_node.Name is None:
                    name = ''
                else:
                    name = str(top_node.Name)
                    if escape_name and not (name.startswith("'") and
                                            name.endswith("'")):
                        if re.search("""[]['"(),:;_]""", name):
                            name = "'%s'" % name.replace("'", "''")
                        else:
                            name = name.replace(' ', '_')
                result.append(name)

                if with_distances and top_node.Length is not None:
                    result[-1] = "%s:%s" % (result[-1], top_node.Length)

                result.append(',')

        len_result = len(result)
        if len_result == 2:  # single node no name
            if semicolon:
                return ";"
            else:
                return ''
        elif len_result == 3:  # single node with name
            if semicolon:
                return "%s;" % result[1]
            else:
                return result[1]
        else:
            if semicolon:
                result[-1] = ';'
            else:
                result.pop(-1)
            return ''.join(result)

    def _ascii_art(self, char1='-', show_internal=True, compact=False):
        LEN = 10
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        namestr = self.Name or ''  # prevents name of NoneType
        if self.Children:
            mids = []
            result = []
            for c in self.Children:
                if c is self.Children[0]:
                    char2 = '/'
                elif c is self.Children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._ascii_art(char2, show_internal, compact)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = np.int(np.trunc((lo + hi) / 2))
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p, l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + namestr + stem[len(namestr)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + namestr], 0)

    def ascii_art(self, show_internal=True, compact=False):
        """Returns a string containing an ascii drawing of the tree.

        Arguments:
        - show_internal: includes internal edge names.
        - compact: use exactly one line per tip.

        Note, this method calls a private recursive function and is not safe
        for large trees.
        """
        (lines, mid) = self._ascii_art(show_internal=show_internal,
                                       compact=compact)
        return '\n'.join(lines)

    ### end formatters ###

    ### distance methods ###

    def accumulate_to_ancestor(self, ancestor):
        """Return the sum of the distance between self and ancestor"""
        accum = 0.0
        curr = self
        while curr is not ancestor:
            if curr.is_root():
                raise NoParentError("Provided ancestor is not in the path")

            if curr.Length is None:
                raise NoLengthError("No length on node %s found!" %
                                    curr.Name or "unnamed")

            accum += curr.Length
            curr = curr.Parent

        return accum

    def distance(self, other):
        """Return the distance between self and other"""
        if self is other:
            return 0.0

        root = self.root()
        lca = root.lowest_common_ancestor([self, other])
        accum = self.accumulate_to_ancestor(lca)
        accum += other.accumulate_to_ancestor(lca)

        return accum

    ### make max distance a property?
    def _set_max_distance(self):
        """Propagate tip distance information up the tree

        This method was originally implemented by Julia Goodrich with the
        intent of being able to determine max tip to tip distances between
        nodes on large trees efficiently. The code has been modified to track
        the specific tips the distance is between
        """
        for n in self.postorder():
            if n.is_tip():
                n.MaxDistTips = [[0.0, n], [0.0, n]]
            else:
                if len(n.Children) == 1:
                    raise TreeError("No support for single descedent nodes")
                else:
                    tip_info = [(max(c.MaxDistTips), c) for c in n.Children]
                    dists = [i[0][0] for i in tip_info]
                    best_idx = np.argsort(dists)[-2:]
                    tip_a, child_a = tip_info[best_idx[0]]
                    tip_b, child_b = tip_info[best_idx[1]]
                    tip_a[0] += child_a.Length or 0.0
                    tip_b[0] += child_b.Length or 0.0
                n.MaxDistTips = [tip_a, tip_b]

    def _get_max_distance_singledesc(self):
        """returns the max distance between any pair of tips

        Also returns the tip names  that it is between as a tuple"""
        distmtx, tip_order = self.tip_tip_distances()
        idx_max = divmod(distmtx.argmax(), distmtx.shape[1])
        max_pair = (tip_order[idx_max[0]].Name, tip_order[idx_max[1]].Name)
        return distmtx[idx_max], max_pair

    def get_max_distance(self):
        """Returns the max tip tip distance between any pair of tips

        Returns (dist, tips)
        """
        if not hasattr(self, 'MaxDistTips'):
            try:
                self._set_max_distance()
            except TreeError:
                return self._get_max_distance_singledesc()

        longest = 0.0
        tips = [None, None]
        for n in self.non_tips(include_self=True):
            tip_a, tip_b = n.MaxDistTips
            dist = (tip_a[0] + tip_b[0])

            if dist > longest:
                longest = dist
                tips = [tip_a[1], tip_b[1]]
        return longest, tips

    def tip_tip_distances(self, endpoints=None, default_length=1):
        """Returns distance matrix between all pairs of tips, and a tip order.

        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """
        all_tips = list(self.tips())
        if endpoints is None:
            tip_order = all_tips
        else:
            tip_order = [self.find(n) for n in endpoints]

        ## linearize all tips in postorder
        # .__start, .__stop compose the slice in tip_order.
        for i, node in enumerate(all_tips):
            node.__start, node.__stop = i, i+1

        # the result map provides index in the result matrix
        result_map = {n.__start: i for i, n in enumerate(tip_order)}
        num_all_tips = len(all_tips)  # total number of tips
        num_tips = len(tip_order)  # total number of tips in result
        result = np.zeros((num_tips, num_tips), float)  # tip by tip matrix
        distances = np.zeros((num_all_tips), float)  # dist from tip to tip

        def update_result():
        # set tip_tip distance between tips of different child
            for child1, child2 in combinations(node.Children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    if tip1 not in result_map:
                        continue
                    t1idx = result_map[tip1]
                    for tip2 in range(child2.__start, child2.__stop):
                        if tip2 not in result_map:
                            continue
                        t2idx = result_map[tip2]
                        result[t1idx, t2idx] = distances[tip1]+distances[tip2]

        for node in self.postorder():
            if not node.Children:
                continue
            ## subtree with solved child wedges
            ### can possibly use np.zeros
            starts, stops = [], []  # to calc ._start and ._stop for curr node
            for child in node.Children:
                if child.Length is not None:
                    child_len = child.Length
                else:
                    child_len = default_length

                distances[child.__start:child.__stop] += child_len

                starts.append(child.__start)
                stops.append(child.__stop)

            node.__start, node.__stop = min(starts), max(stops)

            if len(node.Children) > 1:
                update_result()

        return result + result.T, tip_order

    ### end distance methods ###

    ### comparison methods ###

#   def compare_rfd(self, other, proportion=False):
#       """Calculates the Robinson and Foulds symmetric difference
#
#       Implementation based off of code by Julia Goodrich
#       """
#       t1names = {n.Name for n in self.tips()}
#       t2names = {n.Name for n in other.tips()}
#
#       if t1names != t2names:
#           if t1names < t2names:
#               tree2 = other.shear(t1names)
#           else:
#               tree1 = self.shear(t2names)
#
#       tree1_sets = tree1.subsets()
#       tree2_sets = tree2.subsets()
#
#       not_in_both = tree1_sets ^ tree2_sets
#       total_subsets = len(tree1_sets) + len(tree2_sets)
#
#       dist = len(not_in_both)
#
#       if proportion:
#           dist = dist/float(total_subsets)
#
#       return dist

    def compare_subsets(self, other, exclude_absent_taxa=False):
        """Returns fraction of overlapping subsets where self and other differ.

        Other is expected to be a tree object compatible with PhyloNode.

        Note: names present in only one of the two trees will count as
        mismatches: if you don't want this behavior, strip out the non-matching
        tips first.
        """
        self_sets, other_sets = self.subsets(), other.subsets()

        if exclude_absent_taxa:
            in_both = self.subset() & other.subset()
            self_sets = (i & in_both for i in self_sets)
            self_sets = frozenset({i for i in self_sets if len(i) > 1})
            other_sets = (i & in_both for i in other_sets)
            other_sets = frozenset({i for i in other_sets if len(i) > 1})

        total_subsets = len(self_sets) + len(other_sets)
        intersection_length = len(self_sets & other_sets)

        if not total_subsets:  # no common subsets after filtering, so max dist
            return 1

        return 1 - (2*intersection_length / float(total_subsets))

    def compare_tip_distances(self, other, sample=None, dist_f=distance_from_r,
                              shuffle_f=shuffle):
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).

        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching
        names and matching names is undefined in the tree where they don't
        match, and because we need to reorder the names in the two trees to
        match up the distance matrices).
        """
        self_names = {i.Name: i for i in self.tips()}
        other_names = {i.Name: i for i in other.tips()}
        common_names = frozenset(self_names) & frozenset(other_names)
        common_names = list(common_names)

        if not common_names:
            raise ValueError("No names in common between the two trees.")

        if len(common_names) <= 2:
            return 1  # the two trees must match by definition in this case

        if sample is not None:
            shuffle_f(common_names)
            common_names = common_names[:sample]

        self_nodes = [self_names[k] for k in common_names]
        other_nodes = [other_names[k] for k in common_names]

        self_matrix = self.tip_tip_distances(endpoints=self_nodes)[0]
        other_matrix = other.tip_tip_distances(endpoints=other_nodes)[0]

        return dist_f(self_matrix, other_matrix)

    ### end comparison methods ###

    def index_tree(self):
        """Returns ({node_id:node}, [node_id,first_child,last_child])

        Indexes nodes in-place as n._leaf_index.

        Algorithm is as follows:
        for each node in post-order traversal over tree:
            if the node has children:
                set an index on each child
                for each child with children:
                    add the child and its start and end tips to the result
        """
        id_index = {}
        child_index = []
        curr_index = 0

        for n in self.postorder():
            for c in n.Children:
                c._leaf_index = curr_index
                id_index[curr_index] = c
                curr_index += 1

                if c:
                    #c has children itself, so need to add to result
                    child_index.append((c._leaf_index,
                                        c.Children[0]._leaf_index,
                                        c.Children[-1]._leaf_index))

        # handle root, which should be t itself
        self._leaf_index = curr_index
        id_index[curr_index] = self

        # only want to add to the child_index if self has children...
        if self.Children:
            child_index.append((self._leaf_index,
                                self.Children[0]._leaf_index,
                                self.Children[-1]._leaf_index))

        return id_index, child_index

    def assign_ids(self):
        """Assign topologically stable unique IDs to self"""
        for idx, n in enumerate(self.traverse(include_self=True)):
            n.Id = idx
