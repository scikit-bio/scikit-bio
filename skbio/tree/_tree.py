# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, simplefilter
from operator import ne, gt, itemgetter
from copy import copy, deepcopy
from itertools import chain, combinations
from collections import defaultdict, deque

import numpy as np
import pandas as pd

from skbio._base import SkbioObject
from skbio.stats.distance import DistanceMatrix
from skbio.tree._exception import (
    NoLengthError,
    DuplicateNodeError,
    NoParentError,
    MissingNodeError,
    TreeError,
)
from skbio.util._decorator import (
    classonlymethod,
    deprecated,
    register_aliases,
    aliased,
    params_aliased,
)
from skbio.util._warning import _warn_once
from skbio.io.registry import Read, Write
from ._compare import (
    _check_dist_metric,
    _check_shuffler,
    _topo_dists,
    _path_dists,
)


# ----------------------------------------------------------------------------
# Important note: The TreeNode class has a large number of methods. They are
# organized under several categories, which are defined in this script as well
# as in `doc/source/_templates/TreeNode.rst`, which is a template file for the
# documentation. When methods are added, removed or re-organized, one needs to
# edit the template file to reflect the changes.
# ----------------------------------------------------------------------------


@register_aliases
class TreeNode(SkbioObject):
    r"""Represent a node within a tree.

    A ``TreeNode`` instance stores links from a node to its parent node and optionally
    child nodes. In addition, it can represent the length of the branch connecting
    itself and its parent, and the support of this branch.

    Parameters
    ----------
    name : str or None
        Name of the node. It is common for tips in particular to have names, for
        instance, in a phylogenetic tree where the tips correspond to taxa. Internal
        nodes and the root may also have names.
    length : float, int, or None
        Length of the branch connecting this node to its parent. Can represent
        ellapsed time, amount of mutations, or other measures of evolutionary
        distance.
    support : float, int, or None
        Support value of the branch connecting this node to its parent. Can be
        bootstrap value, posterior probability, or other measures of the confidence or
        frequency of this branch.
    parent : TreeNode or None
        Parent node to which this node is connected. A node without a parent is the
        root of the tree.
    children : list of TreeNode or None
        Child nodes to which this node is connected. A node without any children is a
        tip (leaf) of the tree.

    Notes
    -----
    A tree is a graph in which any two nodes (vertices) are connected by exactly one
    path. The ``TreeNode`` class is capable of representing various tree structures,
    including binary trees, phylogenetic trees, and other hierarchical systems such as
    taxonomies and ontologies. While the class is versatile, many of its terms and
    methods are specifically designed for phylogenetic analysis.

    In scikit-bio, trees are modeled as a collection of interconnected ``TreeNode``
    objects, each representing a single node in the tree. There is no explicit class
    for the entire tree, a clade, or a branch (edge). Instead, a tree is implicitly
    defined by its root node, from which the entire tree can be traversed. Starting
    from any node, one can navigate up to its parent and ancestors, down to its
    children and descendants, or sideways to its siblings.

    The underlying data structure of a tree composed of ``TreeNode`` objects is an
    ordered, rooted tree. However, the ``TreeNode`` class has the flexibility to handle
    unrooted and unordered trees as well, which are common in phylogenetics.

    """

    default_write_format = "newick"

    read = Read()
    write = Write()

    def __init__(
        self, name=None, length=None, support=None, parent=None, children=None
    ):
        self.name = name
        self.length = length
        self.support = support
        self.parent = parent
        self.children = []

        # TODO: `id` doesn't need to be a default attribute.
        self.id = None

        # TODO: This could skip cache clearing.
        if children is not None:
            self.extend(children)

    def __repr__(self):
        r"""Return summary of the tree.

        Returns
        -------
        str
            A summary of this node and all descendants

        Notes
        -----
        This method returns the name of the node and a count of tips and the
        number of internal nodes in the tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c, d)root;"])
        >>> repr(tree)
        '<TreeNode, name: root, internal node count: 1, tips count: 3>'

        """
        nodes = [n for n in self.traverse(include_self=False)]
        n_tips = sum([n.is_tip() for n in nodes])
        n_nontips = len(nodes) - n_tips
        classname = self.__class__.__name__
        name = self.name if self.name is not None else "unnamed"

        return "<%s, name: %s, internal node count: %d, tips count: %d>" % (
            classname,
            name,
            n_nontips,
            n_tips,
        )

    def __str__(self):
        r"""Return a Newick string of self, with names and distances."""
        return str("".join(self.write([])))

    def __iter__(self):
        r"""Iterate over the children of self."""
        return iter(self.children)

    def __len__(self):
        """Return the number of children of self."""
        return len(self.children)

    def __getitem__(self, i):
        r"""Slice the children of self."""
        return self.children[i]

    # ------------------------------------------------
    # Tree copying
    # ------------------------------------------------

    # node attributes that should not be copied
    _exclude_from_copy = {
        "name",
        "length",
        "support",
        "parent",
        "children",
        "id",
        "_tip_cache",
        "_non_tip_cache",
        "_registered_caches",
    }

    def _copy(self, deep, memo):
        """Return a copy of self."""

        # decide deep or shallow copy
        _copy = deepcopy if deep else copy
        _args = [memo] if deep else []

        # node attributes to exclude during copying
        # add any custom attributes that were registered as caches
        exclude_attrs = self._exclude_from_copy
        if hasattr((root := self.root()), "_registered_caches"):
            exclude_attrs = exclude_attrs | root._registered_caches

        # exclude dynamically generated methods
        exclude_attrs = exclude_attrs | {"_write_method"}

        # tree node class (default is TreeNode)
        # this is _possibly_ dangerous, we're assuming the node to copy is
        # of the same class as self, and has the same exclusion criteria.
        # however, it is potentially dangerous to mix TreeNode subclasses
        # within a tree, so...
        treenode = self.__class__

        def __copy_node(node, parent=None):
            """Copy a node."""

            # create a new instance by transferring built-in attributes, which can be
            # directly assigned
            res = treenode(
                name=node.name,
                length=node.length,
                support=node.support,
                parent=parent,
                children=None,
            )
            res.id = node.id

            # copy custom attributes, which may be compound objects therefore need to
            # be copied
            # this method of iteration is slightly faster than
            # `for key in node.__dict__.keys() - exclude_attrs:`
            for key in node.__dict__:
                if key not in exclude_attrs:
                    res.__dict__[key] = _copy(node.__dict__[key], *_args)
            return res

        # start with a copy of self, which will become the root (no parent)
        root = __copy_node(self)
        stack = [[root, self, len(self.children)]]
        stack_append = stack.append

        while stack:
            # check the top node, any children left unvisited?
            top = stack[-1]
            new_top_node, old_top_node, unvisited_children = top

            if unvisited_children:
                top[2] -= 1
                old_child = old_top_node.children[-unvisited_children]
                new_child = __copy_node(old_child, new_top_node)
                new_top_node.children.append(new_child)
                stack_append([new_child, old_child, len(old_child.children)])
            else:
                del stack[-1]
        return root

    def __copy__(self):
        """Return a shallow copy."""
        return self._copy(False, {})

    def __deepcopy__(self, memo):
        """Return a deep copy."""
        return self._copy(True, memo)

    def copy(self, deep=True):
        r"""Return a copy of self using an iterative approach.

        Parameters
        ----------
        deep : bool, optional
            Whether to perform a deep (True, default) or shallow (False) copy of node
            attributes.

            .. versionadded:: 0.6.2

            .. note:: The default value will be changed to False in 0.7.0.

        Returns
        -------
        TreeNode
            A new copy of self.

        .. versionchanged:: 0.6.3
            Node attribute caches will not be copied.

        See Also
        --------
        unrooted_copy

        Notes
        -----
        This method iteratively copies the current node and its descendants. That is,
        if the current node is not the root of the tree, only the subtree below the
        node, instead of the entire tree, will be copied.

        All nodes and their attributes except for caches will be copied. The copies are
        new objects rather than references to the original objects. The distinction
        between deep and shallow copies only applies to each node attribute.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> tree_copy = tree.copy()
        >>> tree_nodes = set([id(n) for n in tree.traverse()])
        >>> tree_copy_nodes = set([id(n) for n in tree_copy.traverse()])
        >>> print(len(tree_nodes.intersection(tree_copy_nodes)))
        0

        """
        return self._copy(deep, {})

    @deprecated("0.6.2", msg="Use `copy` instead.")
    def deepcopy(self):
        r"""Return a deep copy of self using an iterative approach.

        Returns
        -------
        TreeNode
            A new deep copy of self.

        See Also
        --------
        copy

        Notes
        -----
        ``deepcopy`` is equivalent to ``copy`` with ``deep=True``, which is
        currently the default behavior of the latter.

        """
        return self._copy(True, {})

    def subtree(self, tip_list=None):
        r"""Make a copy of the subtree.

        .. deprecated:: 0.6.3
            This method will be removed in version 0.7.0. It was never implemented, and
            its goal can be achieved by :meth:`copy`.

        """
        raise NotImplementedError()

    # ------------------------------------------------
    # Tree navigation
    # ------------------------------------------------

    def is_tip(self):
        r"""Check if the current node is a tip of a tree.

        Returns
        -------
        bool
            Whether the node is a tip.

        See Also
        --------
        is_root
        has_children

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c);"])
        >>> print(tree.is_tip())
        False
        >>> print(tree.find('a').is_tip())
        True

        """
        return not self.children

    def is_root(self):
        r"""Check if the current node is the root of a tree.

        Returns
        -------
        bool
            Whether the node is the root.

        See Also
        --------
        is_tip
        has_children

        Notes
        -----
        A root is defined as a node that has no ``parent``. A tree has exactly one
        root.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c);"])
        >>> print(tree.is_root())
        True
        >>> print(tree.find('a').is_root())
        False

        """
        return self.parent is None

    def has_children(self):
        r"""Check if the current node has any children.

        Returns
        -------
        bool
            Whether the node has at least one child.

        See Also
        --------
        is_tip
        is_root

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c);"])
        >>> print(tree.has_children())
        True
        >>> print(tree.find('a').has_children())
        False

        """
        return not self.is_tip()

    def root(self):
        r"""Return root of the tree which contains `self`.

        Returns
        -------
        TreeNode
            The root of the tree

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> tip_a = tree.find('a')
        >>> root = tip_a.root()
        >>> root == tree
        True

        """
        curr = self
        while not curr.is_root():
            curr = curr.parent
        return curr

    def ancestors(self, include_self=False):
        r"""Return all ancestral nodes from self back to the root.

        Parameters
        ----------
        include_self : bool, optional
            Whether to include the initial node in the path (default: False).

        Returns
        -------
        list of TreeNode
            The path from self toward the root.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> tip = tree.find('a')
        >>> [node.name for node in tip.ancestors()]
        ['c', 'g']
        >>> [node.name for node in tip.ancestors(include_self=True)]
        ['a', 'c', 'g']

        """
        curr = self
        result = [curr] if include_self else []
        result_append = result.append
        while (curr := curr.parent) is not None:
            result_append(curr)
        return result

    def siblings(self):
        r"""Return all nodes that are siblings of the current node.

        Siblings are nodes that are children of the current node's parent, except for
        the current node itself.

        Returns
        -------
        list of TreeNode
            The list of sibling nodes relative to self.

        See Also
        --------
        neighbors

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e,f)g)root;"])
        >>> tip_e = tree.find('e')
        >>> [n.name for n in tip_e.siblings()]
        ['d', 'f']

        """
        try:
            return [x for x in self.parent.children if x is not self]
        except AttributeError:
            return []

    def neighbors(self, ignore=None):
        r"""Return all nodes that are neighbors of the current node.

        Neighbors are nodes that are directly connected to the current node by one
        branch. They usually include parent and children of the current node, if
        present. One may optionally ignore one node from the result.

        Parameters
        ----------
        ignore : TreeNode, optional
            A node to ignore.

        Returns
        -------
        list of TreeNode
            The list of all nodes that are connected to self.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> node_c = tree.find('c')
        >>> [n.name for n in node_c.neighbors()]
        ['a', 'b', 'root']

        """
        if (parent := self.parent) is not None:
            nodes = self.children + [parent]
        else:
            nodes = self.children[:]
        if ignore is None:
            return nodes
        else:
            return [n for n in nodes if n is not ignore]

    @aliased("lowest_common_ancestor")
    @params_aliased([("nodes", "tipnames", "0.6.3", True)])
    def lca(self, nodes=None):
        r"""Find the lowest common ancestor of a list of nodes.

        Parameters
        ----------
        nodes : iterable of TreeNode or str
            Instances or names of the nodes of interest.

        Returns
        -------
        TreeNode
            The lowest common ancestor of the nodes.

        Raises
        ------
        MissingNodeError
            If some nodes cannot be found in the tree.

        Notes
        -----
        Both tips and internal nodes may be provided in ``nodes``. If internal node
        names are provided, it is the user's responsibility to ensure that they are
        unique in the tree.

        This method considers the entire tree rather than the subtree below self.
        Therefore, if some nodes are not descendants of self, the LCA of nodes will be
        ancestral to self.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> nodes = [tree.find('a'), tree.find('b')]
        >>> lca = tree.lca(nodes)
        >>> print(lca.name)
        c
        >>> nodes = [tree.find('a'), tree.find('e')]
        >>> lca = tree.lca(nodes)  # lca is an alias for convience
        >>> print(lca.name)
        root

        """
        if not nodes:
            raise ValueError("No node is specified.")
        nodes = [self.find(x) for x in nodes]
        if len(nodes) == 1:
            return nodes[0]

        # Keep a record of visited nodes, such that the temporary attribute assigned
        # to each node can be cleared after getting LCA.
        visited = []
        visited_append = visited.append

        # Path of the first node to root. LCA must be in this path.
        # A temporary attribute "prev" will be assigned to visited nodes. It represents
        # the previous node in the upward path.
        curr = next(nodes := iter(nodes))
        prev = None
        while curr is not None:
            visited_append(curr)
            curr._prev = prev
            prev = curr
            curr = curr.parent

        # Paths of other nodes to root.
        # The prev attribute no longer needs to record the previous node. It is
        # uniformly set as None. When the path hits a previously visited node, it will
        # stop. If the node is in the first path, its prev becomes None, indicating
        # that it has been visited more than once.
        for curr in nodes:
            while not hasattr(curr, "_prev"):
                visited_append(curr)
                curr._prev = None
                curr = curr.parent
            curr._prev = None

        # walk down from root until a node with prev as None
        curr = prev
        while (prev := curr._prev) is not None:
            curr = prev

        # clear temporary attribute
        for node in visited:
            del node._prev

        return curr

    def _path(self, other):
        r"""Return the path from self to other.

        Parameters
        ----------
        other : TreeNode
            Target node.

        Returns
        -------
        TreeNode
            LCA of self and other.
        list of TreeNode
            self (inclusive) to LCA (exclusive).
        list of TreeNode
            other (inclusive) to LCA (exclusive).

        Notes
        -----
        This algorithm is optimized for finding the LCA of two nodes. Instead,
        :meth:`lca` is optimized for finding the LCA of multiple nodes.

        """
        anc1 = self.ancestors(include_self=True)
        anc2 = other.ancestors(include_self=True)

        # find lowest common ancestor of the two by iterating down from root and
        # stopping at divergence
        # pos is lca's index from root + 1
        lca, pos = None, None
        for i, (n1, n2) in enumerate(zip(reversed(anc1), reversed(anc2))):
            if n1 is n2:
                lca = n1
            else:
                pos = i
                break
        if lca is None:
            raise TreeError("Could not find a path between self and other.")
        if pos is None:
            pos = i + 1

        return lca, anc1[: len(anc1) - pos], anc2[: len(anc2) - pos]

    def path(self, other, include_ends=False):
        r"""Return the list of nodes in the path from self to another node.

        .. versionadded:: 0.6.3

        Parameters
        ----------
        other : TreeNode
            Final node of path.
        include_ends: bool, optional
            Whether to include the initial (self) and final (other) nodes in the list.
            Default is False.

        Returns
        -------
        list
            List of TreeNode objects.

        See Also
        --------
        distance

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -root----|
                 |          /-d
                  \f-------|
                            \-e

        >>> node_1, node_2 = tree.find('a'), tree.find('d')
        >>> path = node_1.path(node_2)
        >>> print(len(path))
        3
        >>> print('-'.join(x.name for x in path))
        c-root-f
        >>> path_2 = node_1.path(node_2, include_ends=True)
        >>> print(len(path_2))
        5
        >>> print('-'.join(x.name for x in path_2))
        a-c-root-f-d

        """
        lca, self_path, other_path = self._path(other)
        path = self_path + [lca] + other_path[::-1]

        # remove initial and final nodes if desired
        if include_ends is False:
            path = path[1:-1]

        return path

    # ------------------------------------------------
    # Tree traversal
    # ------------------------------------------------

    def traverse(self, self_before=True, self_after=False, include_self=True):
        r"""Traverse over tree.

        Parameters
        ----------
        self_before : bool, optional
            Whether to include each node before its descendants (default: True).
        self_after : bool, optional
            Whether to include each node after its descendants (default: False).
        include_self : bool, optional
            Include the initial node if True (default).

        Yields
        ------
        TreeNode
            Visited node.

        See Also
        --------
        preorder
        postorder
        pre_and_postorder
        levelorder
        tips
        non_tips

        Notes
        -----
        This is a depth-first search (DFS). ``self_before`` and ``self_after``
        determine whether a node should be visited before and after traversing its
        children. They are independent. If both True, each internal node (and root)
        will be visited twice. If neither is True, only tips will be returned.

        This method is a generalization of :meth:`preorder`, :meth:`postorder`,
        :meth:`pre_and_postorder` and :meth:`tips`. The default mode
        (``self_before=True, self_after=False``) is equivalent to preorder
        traversal.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.traverse():
        ...     print(node.name)
        g
        c
        a
        b
        f
        d
        e

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
        r"""Perform preorder traversal over tree.

        Parameters
        ----------
        include_self : bool, optional
            Include the initial node if True (default).

        Yields
        ------
        TreeNode
            Visited node.

        See Also
        --------
        traverse
        postorder
        pre_and_postorder
        levelorder

        Notes
        -----
        Preorder traversal visits each node followed by traversing each of its
        children in order. It is also known as NLR (node - left - right). It is
        a depth-first search (DFS). The overall direction of traversal is from
        root to tips.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.preorder():
        ...     print(node.name)
        g
        c
        a
        b
        f
        d
        e

        """
        stack = [self] if include_self else self.children[::-1]
        stack_pop = stack.pop
        stack_extend = stack.extend
        while stack:
            yield (curr := stack_pop())
            if curr.children:
                stack_extend(curr.children[::-1])

    def postorder(self, include_self=True):
        r"""Perform postorder traversal over tree.

        Parameters
        ----------
        include_self : bool, optional
            Include the initial node if True (default).

        Yields
        ------
        TreeNode
            Visited node.

        See Also
        --------
        traverse
        preorder
        pre_and_postorder
        levelorder

        Notes
        -----
        Postorder traversal traverses all children of a node in order before
        visiting the parent node. It is also known as LRN (left - right -
        node). It is a depth-first search (DFS). The overall direction of
        traversal is from tips to root.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.postorder():
        ...     print(node.name)
        a
        b
        c
        d
        e
        f
        g

        """
        # This is somewhat inelegant compared to saving the node and its index
        # on the stack, but is 30% faster in the average case and 3x faster in
        # the worst case (for a comb tree).
        child_index_stack = [0]
        child_index_stack_append = child_index_stack.append
        child_index_stack_pop = child_index_stack.pop
        curr = self
        curr_children = self.children
        curr_children_len = len(curr_children)
        while True:
            curr_index = child_index_stack[-1]
            # if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack_append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                curr_children_len = len(curr_children)
                child_index_stack_pop()
                child_index_stack[-1] += 1

    def pre_and_postorder(self, include_self=True):
        r"""Perform traversal over tree, visiting nodes before and after.

        Parameters
        ----------
        include_self : bool, optional
            Include the initial node if True (default).

        Yields
        ------
        TreeNode
            Visited node.

        See Also
        --------
        traverse
        postorder
        preorder
        levelorder

        Notes
        -----
        Pre- and post-order traversal visits each node before and after
        traversing all children of the node. Therefore, each internal node (and
        root) is visited twice. It is a depth-first search (DFS). The overall
        direction of traversal is from root to tips then back to root.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.pre_and_postorder():
        ...     print(node.name)
        g
        c
        a
        b
        c
        f
        d
        e
        f
        g

        """
        # handle simple case first
        if not self.children:
            if include_self:
                yield self
            return
        child_index_stack = [0]
        child_index_stack_append = child_index_stack.append
        child_index_stack_pop = child_index_stack.pop
        curr = self
        curr_children = self.children
        while True:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            # if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack_append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                child_index_stack_pop()
                child_index_stack[-1] += 1

    def levelorder(self, include_self=True):
        r"""Perform level order traversal over tree.

        Parameters
        ----------
        include_self : bool, optional
            Include the initial node if True (default).

        Yields
        ------
        TreeNode
            Visited node.

        See Also
        --------
        postorder
        preorder
        pre_and_postorder
        traverse

        Notes
        -----
        Level order traversal visits all nodes at each depth from the root
        before visiting nodes at the next depth. It is a breadth-first search
        (BFS). The overall direction of traversal is from root to tips.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.levelorder():
        ...     print(node.name)
        g
        c
        f
        a
        b
        d
        e

        """
        queue = deque([self]) if include_self else deque(self.children)
        queue_popleft = queue.popleft
        queue_extend = queue.extend
        while queue:
            yield (curr := queue_popleft())
            if curr.children:
                queue_extend(curr.children)

    def tips(self, include_self=False):
        r"""Iterate over tips descended from the current node.

        Parameters
        ----------
        include_self : bool, optional
            Whether to include the initial node if it is a tip (default: False).

        Yields
        ------
        TreeNode
            Visited tip.

        See Also
        --------
        non_tips
        postorder

        Notes
        -----
        Nodes are ordered by a postorder traversal of the tree. The order is
        consistent between calls.

        If self is a tip, it won't be yieled unless `include_self` is True.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f);"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        ---------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.tips():
        ...     print(node.name)
        a
        b
        d
        e

        """
        for node in self.postorder(include_self=include_self):
            if not node.children:
                yield node

    def non_tips(self, include_self=False):
        r"""Iterate over non-tip nodes descended from the current node.

        Parameters
        ----------
        include_self : bool, optional
            Whether to include the initial node if it is not a tip (default: False).

        Yields
        ------
        TreeNode
            Visited non-tip node.

        See Also
        --------
        tips
        postorder

        Notes
        -----
        Nodes are ordered by a postorder traversal of the tree. The order is
        consistent between calls.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f);"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        ---------|
                 |          /-d
                  \f-------|
                            \-e

        >>> for node in tree.non_tips():
        ...     print(node.name)
        c
        f

        """
        for node in self.postorder(include_self):
            if node.children:
                yield node

    # ------------------------------------------------
    # Tree manipulation
    # ------------------------------------------------

    def append(self, node, uncache=True):
        r"""Add a node to self's children.

        Parameters
        ----------
        node : TreeNode
            Node to add as a child.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        extend

        Notes
        -----
        This method will add the node to the end of self's children. If the incoming
        node is within another tree, it will be disconnected from its original parent,
        if any, but its children will be preserved. Therefore, this method is able to
        move an entire clade.

        The ``uncache`` parameter applies to both donor and recipient trees.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> root = TreeNode(name="root")
        >>> child1 = TreeNode(name="child1")
        >>> child2 = TreeNode(name="child2")
        >>> root.append(child1)
        >>> root.append(child2)
        >>> print(root)
        (child1,child2)root;
        <BLANKLINE>

        """
        if uncache:
            self.clear_caches()
            node.clear_caches()

        # reconnect the node from its original parent to self
        # this code is similar to `remove`, but it does not return a value
        if node.parent is not None:
            for i, curr_node in enumerate((children := node.parent.children)):
                if curr_node is node:
                    del children[i]
                    break
        node.parent = self
        self.children.append(node)

    def extend(self, nodes, uncache=True):
        r"""Add a list of nodes to self's children.

        Parameters
        ----------
        nodes : iterable of TreeNode
            Nodes to add as children.

            .. versionchanged:: 0.6.2

                Can accept any iterable type in addition to list as input.

        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        append

        Notes
        -----
        This method will remove existing parents of the nodes if they have any, set
        their parents to self, and add the nodes to the end of self's children.

        The ``uncache`` parameter applies to both donor and recipient trees.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> root = TreeNode(name="root")
        >>> root.extend([TreeNode(name="child1"), TreeNode(name="child2")])
        >>> print(root)
        (child1,child2)root;
        <BLANKLINE>

        """
        # make a shallow copy of nodes, which is necessary for working with iterators
        # and containers that are mutable during reconnection (like `children`)
        nodes = list(nodes)
        if uncache:
            self.clear_caches()
            for node in nodes:
                node.clear_caches()

        # reconnect each node from original parent to self; see `append`
        for node in nodes:
            if node.parent is not None:
                for i, curr_node in enumerate((children := node.parent.children)):
                    if curr_node is node:
                        del children[i]
                        break
            node.parent = self
        self.children.extend(nodes)

    def insert(self, node, distance=None, branch_attrs=[], uncache=True):
        r"""Insert a node into the branch connecting self and its parent.

        .. versionadded:: 0.6.2

        Parameters
        ----------
        node : TreeNode
            Node to insert.
        distance : float, int or None, optional
            Distance between self and the insertion point. Must not exceed ``length``
            of self. If None whereas ``length`` is not None, will insert at the
            midpoint of the branch.
        branch_attrs : iterable of str, optional
            Attributes of self that should be transferred to the inserted node
            as they are considered as attributes of the branch. ``support``
            will be automatically included as it is always a branch attribute.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        Raises
        ------
        NoParentError
            If self has no parent.
        ValueError
            If distance is specified but branch has no length.
        ValueError
            If distance exceeds branch length.

        See Also
        --------
        append

        Notes
        -----
        This method will remove the existing parent of the node if any, set its parent
        as self's parent, and set self's parent as the incoming node. The node's index
        position in the parent's children is consistent with that of self prior to
        insertion.

        The ``uncache`` parameter applies to both donor and recipient trees.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:4,d:5)e;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
        -e-------|          \-b
                 |
                  \-d

        >>> tree.find("c").insert(TreeNode("x"))
        >>> print(tree.ascii_art())
                                      /-a
                  /x------- /c-------|
        -e-------|                    \-b
                 |
                  \-d
        >>> tree.find("c").length
        2.0
        >>> tree.find("x").length
        2.0

        """
        if (parent := self.parent) is None:
            raise NoParentError("Self has no parent.")
        if uncache:
            self.clear_caches()

        # detach node from original tree if applicable
        if node.parent is not None:
            node.parent.remove(node, uncache)

        # replace self with node in the parent's list of children
        node.parent = parent
        for i, curr_node in enumerate(parent.children):
            if curr_node is self:
                parent.children[i] = node

        # add self to the beginning of the node's list of children
        self.parent = node
        node.children.insert(0, self)

        # transfer branch attributes to new node
        branch_attrs = set(branch_attrs)
        branch_attrs.add("support")
        branch_attrs.discard("length")
        for attr in branch_attrs:
            setattr(node, attr, getattr(self, attr, None))

        # determine insertion point
        if distance is None:
            if self.length is None:
                node.length = None
            else:
                self.length *= 0.5
                node.length = self.length
        else:
            if self.length is None:
                raise ValueError("Distance is provided but branch has no length.")
            elif distance > self.length:
                raise ValueError("Distance cannot exceed branch length.")
            node.length = self.length - distance
            self.length = distance

    def pop(self, index=-1, uncache=True):
        r"""Remove and return a child node by index position from self.

        Parameters
        ----------
        index : int, optional
            The index position in ``children`` to pop.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode
            The popped child node.

        Raises
        ------
        IndexError
            If the index position does not exist.

        See Also
        --------
        remove
        remove_by_func

        Notes
        -----
        The parent of the popped node will be set to ``None``.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(a,b)c;"])
        >>> print(tree.pop(0))
        a;
        <BLANKLINE>

        """
        if uncache:
            self.clear_caches()
        node = self.children.pop(index)
        node.parent = None
        return node

    def remove(self, node, uncache=True):
        r"""Remove a child node by identity from self.

        Parameters
        ----------
        node : TreeNode
            The node to remove from self's children.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        Returns
        -------
        bool
            True if the node was removed. False if the node is not a child of self.

        See Also
        --------
        pop
        remove_by_func

        Notes
        -----
        The parent of the removed node will be set to None. The removed node and its
        children (if any) still exist, but are disconnected from the tree. Therefore,
        this method is useful for detaching a clade from a tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(a,b)c;"])
        >>> tree.remove(tree.children[0])
        True

        """
        if uncache:
            self.clear_caches()

        # it is necessary to perform removal by identity (`is`), instead of removal by
        # equality (e.g., `self.children.remove(node)`), therefore:
        for i, curr_node in enumerate(self.children):
            if curr_node is node:
                curr_node.parent = None
                del self.children[i]
                return True
        return False

    @aliased("remove_deleted", "0.6.3", True)
    def remove_by_func(self, func, uncache=True):
        r"""Remove nodes of a tree that meet certain criteria.

        Parameters
        ----------
        func : callable
            A function that accepts a ``TreeNode`` and returns True or False, where
            True indicates the node is to be deleted.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        pop
        remove

        Notes
        -----
        This method has the potential to drop entire clades. That is, if an internal
        node is removed, all its descendants are no longer connected to the tree, even
        if they are not explicitly removed.

        This method has the potential to leave single-child internal nodes in the tree,
        which can be further collapsed by :meth:`prune`.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(a,b)c;"])
        >>> tree.remove_by_func(lambda x: x.name == 'b')
        >>> print(tree)
        (a)c;
        <BLANKLINE>

        """
        if uncache:
            self.clear_caches()
        for node in self.traverse(include_self=False):
            if func(node):
                node.parent.remove(node, uncache=False)

    def prune(self, uncache=True):
        r"""Collapse single-child nodes in the tree.

        Internal nodes with only one child will be removed, and direct connections will
        be made from the parent to the child. The branch length of the node will be
        added to the child. The name and properties of the child will override those of
        the parent following the operation.

        Parameters
        ----------
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        shear
        pop
        remove
        remove_by_func

        Notes
        -----
        This method is useful for cleaning up single-child nodes after some nodes were
        removed from a tree.

        If called from an internal node of the tree, only the clade below the node will
        be pruned.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(((a,b)c,(d)e)g,((h,i)j)k)root;"])
        >>> print(tree.ascii_art())
                                      /-a
                            /c-------|
                  /g-------|          \-b
                 |         |
        -root----|          \e------- /-d
                 |
                 |                    /-h
                  \k------- /j-------|
                                      \-i

        >>> tree.prune()
        >>> print(tree.ascii_art())
                                      /-a
                            /c-------|
                  /g-------|          \-b
                 |         |
        -root----|          \-d
                 |
                 |          /-h
                  \j-------|
                            \-i

        """
        if uncache:
            self.clear_caches()

        # build up the list of nodes to remove so the topology is not altered
        # while traversing
        nodes_to_remove = []
        nodes_to_remove_append = nodes_to_remove.append
        for node in self.traverse(include_self=False):
            if len(node.children) == 1:
                nodes_to_remove_append(node)

        # clean up the single children nodes
        for node in nodes_to_remove:
            child = node.children[0]
            if child.length is None or node.length is None:
                child.length = child.length or node.length
            else:
                child.length += node.length
            if (parent := node.parent) is not None:
                # TODO: replace the original node's index position, rather than append
                # to the end.
                parent.append(child, uncache=False)
                parent.remove(node, uncache=False)

        # If there is a single descendent from the root, the root will adopt the
        # child's properties. We can't "delete" the root as that would be deleting
        # self.
        if len(self.children) == 1:
            child = self.children[0]
            if child.length is None or self.length is None:
                self.length = self.length or child.length
            else:
                self.length += child.length
            for key, value in child.__dict__.items():
                if key not in ("length", "parent", "children"):
                    self.__dict__[key] = value
            self.remove(child, uncache=False)
            self.extend(child.children, uncache=False)

    def shear(self, names, strict=True, prune=True, inplace=False, uncache=True):
        r"""Refine a tree such that it just has the desired tip names.

        Parameters
        ----------
        names : iterable of str
            The tip names on the tree to keep.
        strict : bool, optional
            In case some names are not found in the tree, whether to raise an error
            (True, default) or to refine the tree to the found names only (False).

            .. versionadded:: 0.6.3

        prune : bool, optional
            Whether to collapse single-child nodes after shearing by calling
            :meth:`prune` (default: True).

            .. versionadded:: 0.6.3

        inplace : bool, optional
            Whether to modify the tree in place (True) or to create a modified copy of
            the tree (False, default).

            .. versionadded:: 0.6.3

        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`. Only applicable when ``inplace`` is True.

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode, optional
            The resulting tree (if ``inplace`` is False).

        Raises
        ------
        ValueError
            If one or more names do not exist in the tree and ``strict`` is True.

        See Also
        --------
        prune
        remove
        remove_by_func

        Notes
        -----
        This method is useful for reducing a large tree to a relevant subset of taxa.

        If called from an internal node of the tree, only the clade below the node will
        be refined, and the copy of the tree (when ``inplace`` is False) will only
        include the clade.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(((a,b)c,(d,e)f)g,(h,i)j)root;"])
        >>> print(tree.ascii_art())
                                      /-a
                            /c-------|
                           |          \-b
                  /g-------|
                 |         |          /-d
                 |          \f-------|
        -root----|                    \-e
                 |
                 |          /-h
                  \j-------|
                            \-i

        >>> sheared = tree.shear(['a', 'd', 'h'])
        >>> print(sheared.ascii_art())
                            /-a
                  /g-------|
        -root----|          \-d
                 |
                  \-h

        """
        names = set(names)
        if strict and not names.issubset(self.subset()):
            raise ValueError("Names are not a subset of the tree.")

        # modify (sub)tree in place
        if inplace:
            tree = self
            if uncache:
                tree.clear_caches()

            # temporarily disconnect subtree from parent
            curr_parent = tree.parent
            tree.parent = None

        # make a copy of (sub)tree
        else:
            tree = self.copy()

        # mark desired tips and their ancestors
        marked = set()
        marked_add = marked.add
        for tip in tree.tips():
            if tip.name in names:
                marked_add(tip)

                # see also `tip.ancestors`, but the following code stops early if it
                # doesn't need to reach root
                anc = tip.parent
                while anc is not None:
                    if anc in marked:
                        break
                    marked_add(anc)
                    anc = anc.parent

        # TODO: This `list` can potentially be removed to save unnecessary removals
        # within clades that are already removed
        for node in list(tree.traverse()):
            if node not in marked:
                node.parent.remove(node, uncache=False)

        # remove single-child nodes
        if prune:
            tree.prune(uncache=False)

        # reconnect subtree to parent
        if inplace:
            tree.parent = curr_parent
        else:
            return tree

    def unpack(self, uncache=True):
        """Unpack an internal node in place.

        Parameters
        ----------
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        Notes
        -----
        This method sequentially: 1) elongates child nodes by branch length
        of self (omit if there is no branch length), 2) removes self from
        parent node, and 3) grafts child nodes to parent node.

        Raises
        ------
        ValueError
            If input node is root or tip.

        See Also
        --------
        unpack_by_func
        prune

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(['((c:2.0,d:3.0)a:1.0,(e:2.0,f:1.0)b:2.0);'])
        >>> tree.find('b').unpack()
        >>> print(tree)
        ((c:2.0,d:3.0)a:1.0,e:4.0,f:3.0);
        <BLANKLINE>

        """
        if self.is_root():
            raise TreeError("Cannot unpack root.")
        if self.is_tip():
            raise TreeError("Cannot unpack tip.")
        if uncache:
            self.clear_caches()
        parent = self.parent
        blen = self.length or 0.0
        for child in self.children:
            clen = child.length or 0.0
            child.length = clen + blen or None
        parent.remove(self, uncache=False)
        parent.extend(self.children, uncache=False)

    def unpack_by_func(self, func, uncache=True):
        """Unpack internal nodes of a tree that meet certain criteria.

        Parameters
        ----------
        func : callable
            A function that accepts a ``TreeNode`` and returns True or False, where
            True indicates the node is to be unpacked.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        unpack
        prune

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(['((c:2,d:3)a:1,(e:1,f:2)b:2);'])
        >>> tree.unpack_by_func(lambda x: x.length <= 1)
        >>> print(tree)
        ((e:1.0,f:2.0)b:2.0,c:3.0,d:4.0);
        <BLANKLINE>
        >>> tree = TreeNode.read(['(((a,b)85,(c,d)78)75,(e,(f,g)64)80);'])
        >>> tree.assign_supports()
        >>> tree.unpack_by_func(lambda x: x.support < 75)
        >>> print(tree)
        (((a,b)85,(c,d)78)75,(e,f,g)80);
        <BLANKLINE>

        """
        if uncache:
            self.clear_caches()
        nodes_to_unpack = []
        nodes_to_unpack_append = nodes_to_unpack.append
        for node in self.non_tips(include_self=False):
            if func(node):
                nodes_to_unpack_append(node)
        for node in nodes_to_unpack:
            node.unpack(uncache=False)

    def bifurcate(self, insert_length=None, include_self=True, uncache=True):
        r"""Convert the tree into a bifurcating tree.

        All nodes that have more than two children will have additional intermediate
        nodes inserted to ensure that every node has only two children.

        Parameters
        ----------
        insert_length : int, optional
            The branch length assigned to all inserted nodes.
        include_self : bool, optional
            If False, will not convert the current node. This is useful for keeping an
            unrooted tree unrooted. Default is True.

            .. versionadded:: 0.6.3

        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        prune
        is_bifurcating

        Notes
        -----
        This method does not modify single-child nodes. These nodes can be collapsed
        using :meth:`prune` prior to this method to create a strictly bifurcating tree.

        This method modifies the subtree under the current node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b,g,h)c,(d,e)f)root;"])
        >>> print(tree.ascii_art())
                            /-a
                           |
                           |--b
                  /c-------|
                 |         |--g
                 |         |
        -root----|          \-h
                 |
                 |          /-d
                  \f-------|
                            \-e

        >>> tree.bifurcate()
        >>> print(tree.ascii_art())
                            /-h
                  /c-------|
                 |         |          /-g
                 |          \--------|
                 |                   |          /-a
        -root----|                    \--------|
                 |                              \-b
                 |
                 |          /-d
                  \f-------|
                            \-e

        """
        if uncache:
            self.clear_caches()
        treenode = self.__class__
        for node in self.traverse(include_self=include_self):
            if len(node.children) > 2:
                stack = node.children
                while len(stack) > 2:
                    ind = stack.pop()
                    interm = treenode(length=insert_length, children=stack[:])
                    node.append(interm, uncache=False)
                    for child in stack:
                        node.remove(child, uncache=False)
                    node.extend([ind, interm], uncache=False)

    @params_aliased([("shuffler", "shuffle_f", "0.6.3", True)])
    def shuffle(self, k=None, names=None, shuffler=None, n=1):
        r"""Randomly shuffle tip names of the tree.

        Parameters
        ----------
        k : int, optional
            The number of tips to shuffle. If provided, this number of tips will be
            randomly selected by ``shuffle_f``, and only those names will be shuffled.
            Conflicts with ``names``.
        names : list, optional
            The specific tip names to shuffle. Conflicts with ``k``.
        shuffler : int, np.random.Generator or callable, optional
            Shuffling function, which must accept a list and modify in place. Default
            is the :meth:`shuffle <numpy.random.Generator.shuffle>` method of a NumPy
            random generator. If an integer is provided, a random generator will be
            constructed using this number as the seed.

            .. versionchanged:: 0.6.3
                Switched to NumPy's new random generator. Can accept a random seed or
                random generator instance.

        n : int, optional
            The number of iterations to perform. Must be a positive integer. Default
            is 1. If None or ``np.inf``, iterations will be infinite.

            .. versionchanged:: 0.6.3
                Can accept None.

        Yields
        ------
        TreeNode
            Tree with shuffled tip names.

        Raises
        ------
        ValueError
            If ``k`` < 2 or ``n`` < 1.
        ValueError
            If both ``k`` and ``names`` are specified.
        MissingNodeError
            If ``names`` is specified but one of the names cannot be found.

        See Also
        --------
        numpy.random.Generator.shuffle

        Notes
        -----
        This method does not create copies of the tree. Instead, tip names are shuffled
        in place in the original tree and the tree is yielded prior to the next round
        of shuffling. Tree caches will be cleared prior to shuffling.

        ``k`` and ``names`` cannot be specified at the same time. If neither ``k`` nor
        ``names`` are provided, all tips will be shuffled.

        Examples
        --------
        Shuffle the names of a 4-tip tree for 5 times:

        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b),(c,d));"])
        >>> for shuffled in tree.shuffle(shuffler=42, n=5):
        ...     print(shuffled)
        ((d,c),(b,a));
        <BLANKLINE>
        ((a,b),(d,c));
        <BLANKLINE>
        ((a,c),(d,b));
        <BLANKLINE>
        ((d,b),(a,c));
        <BLANKLINE>
        ((a,c),(d,b));
        <BLANKLINE>

        """
        if k is not None:
            if k < 2:
                raise ValueError("k must be None or >= 2.")
            if names is not None:
                raise ValueError("k and names cannot be specified at the same time.")
        if n is None:
            n = np.inf
        elif n < 1:
            raise ValueError("n must be > 0.")

        # determine shuffling function
        shuffler = _check_shuffler(shuffler)

        # determine tip names to shuffle
        if names is not None:
            tips = [self.find(x) for x in names]
        else:
            tips = list(self.tips())
            if k is not None:
                shuffler(tips)
                tips = tips[:k]
            names = [x.name for x in tips]

        # since the names are being shuffled, the caches are no longer reliable
        self.clear_caches()

        # iteratively shuffle tip names and yield tree
        counter = 0
        while counter < n:
            shuffler(names)
            for tip, name in zip(tips, names):
                tip.name = name
            yield self
            counter += 1

    # ------------------------------------------------
    # Tree rerooting
    # ------------------------------------------------

    def unroot(self, side=None, uncache=True):
        r"""Convert a rooted tree into unrooted.

        .. versionadded:: 0.6.2

        Parameters
        ----------
        side : int, optional
            Which basal node (i.e., children of root) will be elevated to root. Must be
            0 or 1. If not provided, will elevate the first basal node that is not a
            tip. The choice won't impact tree topology.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        root
        root_at

        Notes
        -----
        In scikit-bio, every tree has a root node. A tree is considered as
        "rooted" if its root node has exactly two children. In contrast, an
        "unrooted" tree may have three (the most common case), one, or more
        than three children attached to its root node. This method will not
        modify the tree if it is already unrooted.

        This method unroots a tree by trifucating its root. Specifically, it
        removes one of the two basal nodes of the tree (i.e., children of the
        root), transfers the name of the removed node to the root, and
        re-attaches the removed node's children to the root. Additionally, the
        removed node's branch length, if available, will be added to the other
        basal node's branch. The outcome appears as if the root is removed
        and the two basal nodes are directly connected.

        The choice of the basal node to be elevated affects the positioning of
        the resulting tree, but does not affect its topology from a
        phylogenetic perspective, as it is considered as unrooted.

        This method manipulates the tree in place. There is no return value.

        .. note:: In the case where the basal node has just one child, the
            resulting tree will still appear rooted as it has two basal nodes.
            To avoid this scenario, call :meth:`prune` to remove all one-child
            internal nodes.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(['(((a,b)c,(d,e)f)g,(h,i)j)k;'])
        >>> print(tree.ascii_art())
                                      /-a
                            /c-------|
                           |          \-b
                  /g-------|
                 |         |          /-d
                 |          \f-------|
        -k-------|                    \-e
                 |
                 |          /-h
                  \j-------|
                            \-i

        >>> tree.unroot()
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
                 |
                 |          /-d
        -g-------|-f-------|
                 |          \-e
                 |
                 |          /-h
                  \j-------|
                            \-i

        """
        # return original tree if already unrooted
        root = self.root()
        if len(bases := root.children) != 2:
            return root

        if uncache:
            root.clear_caches()

        # choose a basal node to elevate
        if side is None:
            side = 1 if (bases[0].is_tip() and not bases[1].is_tip()) else 0
        chosen, other = bases[side], bases[1 - side]

        # remove chosen node and re-attach its children to root
        chosen.parent = None
        for child in chosen.children:
            child.parent = root
        if side:
            root.children = [other] + chosen.children
        else:
            root.children = chosen.children + [other]

        # transfer basal node's attributes to root
        for key, value in chosen.__dict__.items():
            if key not in ("length", "support", "parent", "children"):
                root.__dict__[key] = value

        # add branch length to the other basal node
        if (L := chosen.length) is not None:
            if other.length is not None:
                other.length += L
            else:
                other.length = L

    def unrooted_copy(
        self,
        parent=None,
        branch_attrs={"name", "length", "support"},
        root_name="root",
        deep=False,
        exclude_attrs=None,
    ):
        r"""Walk the tree unrooted-style and return a copy.

        Parameters
        ----------
        parent : TreeNode or None
            Direction of walking (from parent to self). If specified, walking to the
            parent will be prohibited.

        branch_attrs : set of str, optional
            Attributes of ``TreeNode`` objects that should be considered as branch
            attributes during the operation.

            .. versionadded:: 0.6.2

            .. note:: ``name`` will be removed from the default in 0.7.0, as it is
                usually considered as an attribute of the node instead of the branch.

        root_name : str or None, optional
            Name for the new root node, if it doesn't have one.

            .. versionadded:: 0.6.2

            .. note:: This parameter will be removed in 0.7.0, and the root node will
                not be renamed.

        deep : bool, optional
            Whether to perform a shallow (False, default) or deep (True) copy of node
            attributes.

            .. versionadded:: 0.6.2

        exclude_attrs : set, optional
            Node attributes that should not be copied. If None (default), the caches
            will be excluded. This parameter keeps a memo during recursive copying for
            efficiency. It should not be customized by the user unless absolutely
            needed.

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode
            A new copy of the tree rooted at the given node.

            .. versionchanged:: 0.6.2

                Node attributes other than name and length will also be copied.

        Warnings
        --------
        The default behavior of ``unrooted_copy`` is subject to change in 0.7.0. The
        new default behavior can be achieved by specifying
        ``branch_attrs={"length", "support"}, root_name=None``.

        See Also
        --------
        copy
        unrooted_move

        Notes
        -----
        This method recursively walks a tree from a given node in an unrooted style
        (i.e., directions of branches are not assumed), and copies each node it
        visits, such that the copy of the given node becomes the root node of a new
        tree and the copies of all other nodes are re-positioned accordingly, whereas
        the topology of the new tree will be identical to the existing one.

        Nodes attributes except for caches will be copied to the new tree. Attributes
        in ``branch_attrs`` will be transferred to the node at the other end of a
        branch if the branch is flipped in the new tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> new_tree = tree.find('d').unrooted_copy()
        >>> print(new_tree)
        (b,c,(a,((f,g)h)e)d)root;
        <BLANKLINE>

        """
        # future warning
        if branch_attrs == {"name", "length", "support"} and root_name == "root":
            _warn_once(
                self.__class__.unrooted_copy,
                FutureWarning,
                "The default behavior of `unrooted_copy` is subject to change in "
                "0.7.0. The new default behavior can be achieved by specifying "
                '`branch_attrs={"length", "support"}, root_name=None`.',
            )

        # determine copy mode
        _copy = deepcopy if deep else copy

        # determine node attributes to exclude
        if exclude_attrs is None:
            exclude_attrs = self._exclude_from_copy
            if hasattr((root := self.root()), "_registered_caches"):
                exclude_attrs = exclude_attrs | root._registered_caches

        # identify neighbors (adjacent nodes) of self, excluding the incoming node
        neighbors = self.neighbors(ignore=parent)

        # recursively copy each neighbor; they will become outgoing nodes (children)
        children = [
            c.unrooted_copy(
                parent=self,
                branch_attrs=branch_attrs,
                root_name=root_name,
                deep=deep,
                exclude_attrs=exclude_attrs,
            )
            for c in neighbors
        ]

        # identify node from which branch attributes should be transferred
        # 1. starting point (becomes root)
        if parent is None:
            other = None
        # 2. walk up (parent becomes child)
        elif parent.parent is self:
            other = parent
        # 3. walk down (retain the same order)
        else:
            other = self

        # create a new node and attach children to it, see also `copy`
        attrs = {
            x: (
                (None if other is None else getattr(other, x))
                if x in branch_attrs
                else getattr(self, x)
            )
            for x in ("name", "length", "support")
        }
        result = self.__class__(**attrs, children=children)

        # transfer attributes to the new node, see also `copy`
        for key in self.__dict__:
            if key not in exclude_attrs:
                source = other if key in branch_attrs else self
                if source is not None and key in source.__dict__:
                    result.__dict__[key] = _copy(source.__dict__[key])

        # name the new root
        if root_name and parent is None and result.name is None:
            result.name = root_name

        return result

    @deprecated(
        "0.6.2",
        msg="Because it generates a redundant copy of the tree. Use `unrooted_copy` "
        "instead.",
    )
    def unrooted_deepcopy(self, parent=None):
        r"""Walk the tree unrooted-style and returns a new deepcopy.

        Parameters
        ----------
        parent : TreeNode or None
            Direction of walking (from parent to self). If specified, walking
            to the parent will be prohibited.

        Returns
        -------
        TreeNode
            A new copy of the tree rooted at the given node.

        See Also
        --------
        copy
        unrooted_copy
        root_at

        Notes
        -----
        Perform a deepcopy of self and return a new copy of the tree as an
        unrooted copy. This is useful for defining a new root of the tree.

        This method calls :meth:`unrooted_copy` which is recursive.

        """
        root = self.root()
        root.assign_ids()

        new_tree = root.copy()
        new_tree.assign_ids()

        new_tree_self = new_tree.find_by_id(self.id)
        return new_tree_self.unrooted_copy(parent, deep=True)

    def unrooted_move(
        self,
        branch_attrs={"length", "support"},
        uncache=True,
    ):
        r"""Walk the tree unrooted-style and rearrange it.

        .. versionadded:: 0.6.2

        .. versionchanged:: 0.6.3
            The underlying algorithm is now iterative instead of recursive, therefore
            won't be constrained by Python's maximum recursion limit when working with
            large trees. Parameter ``parent`` was removed as it is no longer needed.

        Parameters
        ----------
        parent : TreeNode or None
            Direction of walking (from parent to self). If specified, walking
            to the parent will be prohibited.
        branch_attrs : set of str, optional
            Attributes of ``TreeNode`` objects that should be considered as
            branch attributes during the operation.
        uncache : bool, optional
            Whether to clear caches of the tree if present (default: True). See
            :meth:`details <has_caches>`.

            .. versionadded:: 0.6.3

        See Also
        --------
        root_at
        unrooted_copy

        Notes
        -----
        This method walks a tree from a given node in an unrooted style (i.e.,
        directions of branches are not assumed). It rerranges the tree such that
        the given node becomes the root node and all other nodes are re-positioned
        accordingly, whereas the topology remains the same.

        This method manipulates the tree in place. There is no return value.
        The new tree should be referred to by the node where the operation
        started, as it has become the new root node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> new_root = tree.find('d')
        >>> new_root.unrooted_move()
        >>> print(new_root)
        (b,c,(a,((f,g)h)i)e)d;
        <BLANKLINE>

        """
        if uncache:
            self.clear_caches()

        # This algorithm uses an iterative approach to avoid the maximum recursion
        # limit imposed by Python. Two rounds of iterations are involved to 1) flip
        # the tree and to 2) reconstruct the upward branches.
        self.old_child = None
        curr = self
        parent = None

        # move up in the original tree and move parent to children
        while (old_parent := curr.parent) is not None:
            for i, child in enumerate(old_parent.children):
                if child is curr:
                    old_parent.children.pop(i)
                    break
            curr.children.append(old_parent)
            old_parent.old_child = curr
            curr.parent = parent
            parent = curr
            curr = old_parent

        # move up in the new tree and rebuild parent connection
        while curr.old_child is not None:
            child = curr.old_child
            curr.parent = child
            for attr in branch_attrs:
                setattr(curr, attr, getattr(child, attr, None))
            del curr.old_child
            curr = child

        del self.old_child
        for attr in branch_attrs:
            setattr(self, attr, None)

    def root_at(
        self,
        node=None,
        above=False,
        reset=False,
        branch_attrs=["name"],
        root_name="root",
        inplace=False,
    ):
        r"""Reroot the tree at the provided node.

        This is useful for positioning a tree with an orientation that reflects
        knowledge of the true root location.

        Parameters
        ----------
        node : TreeNode or str, optional
            The node to root at. Can either be a node object or the name of the node.
            If not provided, will root at self. If a root node provided, will return
            the original tree.

            .. versionchanged:: 0.6.2

                Becomes optional.

        above : bool, float, or int, optional
            Whether and where to insert a new root node. If False (default), the target
            node will serve as the root node. If True, a new root node will be created
            and inserted at the midpoint of the branch connecting the target node and
            its parent. If a number, the new root will be inserted at this distance
            from the target node. The number ranges between 0 and branch length.

            .. versionadded:: 0.6.2

        reset : bool, optional
            Whether to remove the original root of a rooted tree before performing the
            rerooting operation. Default is False.

            .. versionadded:: 0.6.2

            .. note:: The default value will be set as True in 0.7.0.

        branch_attrs : iterable of str, optional
            Attributes of each node that should be considered as attributes of the
            branch connecting the node to its parent. This is important for the correct
            rerooting operation. "length" and "support" will be automatically included
            as they are always branch attributes.

            .. versionadded:: 0.6.2

            .. note:: ``name`` will be removed from the default in 0.7.0, as it is
                usually considered as an attribute of the node instead of the branch.

        root_name : str or None, optional
            Name for the root node, if it doesn't already have one.

            .. versionadded:: 0.6.2

            .. note:: The default value will be set as ``None`` in 0.7.0.

        inplace : bool, optional
            Whether to reroot the tree in place (True) or to create a rerooted copy of
            the tree (False, default).

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode
            A tree rooted at the give node.

        Warnings
        --------
        The default behavior of ``root_at`` is subject to change in 0.7.0. The
        new default behavior can be achieved by specifying ``reset=True,
        branch_attrs=[], root_name=None``.

        See Also
        --------
        unrooted_copy
        unrooted_move
        unroot

        Notes
        -----
        The specified node will be come the root of the new tree.

        Tree caches (see :meth:`details <has_caches>`) will not be retained in the
        returned tree. In in-place mode, they will be cleared prior to rerooting. In
        copying mode, they will not be copied to the new tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(((a,b)c,(d,e)f)g,h)i;"])
        >>> print(tree.ascii_art())
                                      /-a
                            /c-------|
                           |          \-b
                  /g-------|
                 |         |          /-d
        -i-------|          \f-------|
                 |                    \-e
                 |
                  \-h

        Use the given node as the root node. This will typically create an
        unrooted tree (i.e., root node has three children).

        >>> t1 = tree.root_at("c", branch_attrs=[])
        >>> print(t1)
        (a,b,((d,e)f,(h)i)g)c;
        <BLANKLINE>
        >>> print(t1.ascii_art())
                  /-a
                 |
                 |--b
        -c-------|
                 |                    /-d
                 |          /f-------|
                  \g-------|          \-e
                           |
                            \i------- /-h

        Insert a new root node into the branch above the given node. This will
        create a rooted tree (i.e., root node has two children).

        >>> t2 = tree.root_at("c", above=True, branch_attrs=[])
        >>> print(t2)
        ((a,b)c,((d,e)f,(h)i)g)root;
        <BLANKLINE>
        >>> print(t2.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -root----|
                 |                    /-d
                 |          /f-------|
                  \g-------|          \-e
                           |
                            \i------- /-h

        """
        # future warning
        if reset is False and branch_attrs == ["name"] and root_name == "root":
            _warn_once(
                self.__class__.root_at,
                FutureWarning,
                "The default behavior of `root_at` is subject to change in 0.7.0. "
                "The new default behavior can be achieved by specifying "
                "`reset=True, branch_attrs=[], root_name=None`.",
            )

        # locate to-be root node
        tree = self.root()
        if node is None:
            node = self
        elif isinstance(node, str):
            node = tree.find(node)

        # return if already rooted
        if node.is_root():
            return node.copy()

        # check if tree is rooted
        if reset and len(tree.children) != 2:
            reset = False

        # Prior to rerooting, the tree may need to be manipulated to remove the
        # original root and/or to insert a new root node

        # For optimal performance (copying the tree only once), the following code
        # considers three scenarios:
        # 1. In-place mode: Just manipulate the tree if needed, then call
        #    `unrooted_move`.
        # 2. Copying mode, tree doesn't need to be manipulated: directly call
        #    `unrooted_copy`
        # 3. Copying mode, tree needs to be manipulated: Make a copy of the
        #    tree, manipulate, then call `unrooted_move`.

        to_copy = False
        if not inplace:
            if reset or above is not False:
                tree.assign_ids()
                new_tree = tree.copy()
                new_tree.assign_ids()
                node = new_tree.find_by_id(node.id)
                tree = new_tree
            else:
                to_copy = True

        # Clear caches, since root node will be different and caches are going to be
        # useless regardless.
        else:
            tree.clear_caches()

        # Remove original root. We need to make sure the node itself is not the basal
        # node that gets removed.
        if reset:
            side = None
            for i, base in enumerate(tree.children):
                if node is base:
                    side = 1 - i
                    break
            tree.unroot(side)

        # insert a new root node into the branch above
        if above is not False:
            to_insert = node.__class__()
            distance = None if above is True else above
            node.insert(to_insert, distance, branch_attrs, uncache=False)
            node = to_insert

        branch_attrs = set(branch_attrs)
        branch_attrs.update(["length", "support"])

        # rotate the tree to position the new root
        if to_copy:
            return node.unrooted_copy(branch_attrs=branch_attrs, root_name=root_name)
        else:
            node.unrooted_move(branch_attrs=branch_attrs, uncache=False)
            if root_name and node.name is None:
                node.name = root_name
            return node

    def root_at_midpoint(
        self, reset=False, branch_attrs=["name"], root_name="root", inplace=False
    ):
        r"""Reroot the tree at the midpoint of the two tips farthest apart.

        Parameters
        ----------
        reset : bool, optional
            Whether to remove the original root of a rooted tree before performing
            the rerooting operation. Default is False.

            .. versionadded:: 0.6.2

            .. note:: The default value will be set as True in 0.7.0.

        branch_attrs : iterable of str, optional
            Attributes of each node that should be considered as attributes of
            the branch connecting the node to its parent. This is important for
            the correct rerooting operation. "length" and "support" will be
            automatically included as they are always branch attributes.

            .. versionadded:: 0.6.2

            .. note:: ``name`` will be removed from the default in 0.7.0, as
                it is usually considered as an attribute of the node instead of
                the branch.

        root_name : str or None, optional
            Name for the new root node, if it doesn't have one.

            .. versionadded:: 0.6.2

            .. note:: The default value will be set as ``None`` in 0.7.0.

        inplace : bool, optional
            Whether to reroot the tree in place (True) or to create a rerooted copy of
            the tree (False, default).

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode
            A tree rooted at its midpoint.

        Raises
        ------
        TreeError
            If a tip ends up being the mid point.
        LengthError
            Midpoint rooting requires `length` and will raise (indirectly) if
            evaluated nodes don't have length.

        Warnings
        --------
        The default behavior of ``root_at_midpoint`` is subject to change in
        0.7.0. The new default behavior can be achieved by specifying
        ``reset=True, branch_attrs=[], root_name=None``.

        See Also
        --------
        root_at
        unrooted_copy

        Notes
        -----
        The midpoint rooting (MPR) method was originally described in [1]_.

        Tree caches (see :meth:`details <has_caches>`) will not be retained in the
        returned tree. In in-place mode, they will be cleared prior to rerooting. In
        copying mode, they will not be copied to the new tree.

        References
        ----------
        .. [1] Farris, J. S. (1972). Estimating phylogenetic trees from
           distance matrices. The American Naturalist, 106(951), 645-668.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5,g:1)h;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
                 |
        -h-------|          /-d
                 |-f-------|
                 |          \-e
                 |
                  \-g

        >>> t = tree.root_at_midpoint(branch_attrs=[])
        >>> print(t)
        ((d:3.0,e:4.0)f:2.0,((a:1.0,b:1.0)c:2.0,g:1.0)h:3.0)root;
        <BLANKLINE>
        >>> print(t.ascii_art())
                            /-d
                  /f-------|
                 |          \-e
        -root----|
                 |                    /-a
                 |          /c-------|
                  \h-------|          \-b
                           |
                            \-g

        """
        # future warning
        if reset is False and branch_attrs == ["name"] and root_name == "root":
            _warn_once(
                self.__class__.root_at_midpoint,
                FutureWarning,
                "The default behavior of `root_at_midpoint` is subject to change in "
                "0.7.0. The new default behavior can be achieved by specifying "
                "`reset=True, branch_attrs=[], root_name=None`.",
            )

        tree = self.root()
        if inplace:
            tree.clear_caches()
        else:
            tree = tree.copy()

        if reset:
            tree.unroot(uncache=False)

        max_dist, tips = tree.maxdist()
        half_max_dist = max_dist / 2.0

        if max_dist == 0.0:
            return tree

        tip1 = tree.find(tips[0])
        tip2 = tree.find(tips[1])
        lca = tree.lca([tip1, tip2])

        if tip1.depth(lca) > half_max_dist:
            climb_node = tip1
        else:
            climb_node = tip2

        dist_climbed = 0.0
        while dist_climbed + climb_node.length < half_max_dist:
            dist_climbed += climb_node.length
            climb_node = climb_node.parent

        # case 1: midpoint is at the climb node's parent
        # make the parent node as the new root
        if dist_climbed + climb_node.length == half_max_dist:
            new_root = climb_node.parent

        # case 2: midpoint is on the climb node's branch to its parent
        # insert a new root node into the branch
        else:
            new_root = tree.__class__()
            climb_node.insert(new_root, half_max_dist - dist_climbed, uncache=False)
            # TODO: Here, `branch_attrs` should be added to `insert`. However, this
            # will cause a backward-incompatible behavior. This change will be made
            # in version 0.7.0, along with the removal of `name` from the default of
            # `branch_attrs`.

        branch_attrs = set(branch_attrs)
        branch_attrs.update(["length", "support"])
        new_root.unrooted_move(branch_attrs=branch_attrs, uncache=False)
        if root_name and new_root.name is None:
            new_root.name = root_name
        return new_root

    def root_by_outgroup(
        self,
        outgroup,
        above=True,
        reset=True,
        branch_attrs=[],
        root_name=None,
        inplace=False,
    ):
        r"""Reroot the tree with a given set of taxa as outgroup.

        .. versionadded:: 0.6.2

        Parameters
        ----------
        outgroup : iterable of str
            Taxon set to serve as outgroup. Must be a proper subset of taxa in the
            tree. The tree will be rooted at the lowest common ancestor (LCA) of the
            outgroup.
        above : bool, float, or int, optional
            Whether and where to insert a new root node. If False, the LCA will serve
            as the root node. If True (default), a new root node will be created and
            inserted at the midpoint of the branch connecting the LCA and its parent
            (i.e., the midpoint between outgroup and ingroup). If a number between 0
            and branch length, the new root will be inserted at this distance from the
            LCA.
        reset : bool, optional
            Whether to remove the original root of a rooted tree before performing the
            rerooting operation. Default is True.
        branch_attrs : iterable of str, optional
            Attributes of each node that should be considered as attributes of the
            branch connecting the node to its parent. This is important for the correct
            rerooting operation. "length" and "support" will be automatically included
            as they are always branch attributes.
        root_name : str or None, optional
            Name for the root node, if it doesn't already have one.
        inplace : bool, optional
            Whether to reroot the tree in place (True) or to create a rerooted copy of
            the tree (False, default).

            .. versionadded:: 0.6.3

        Returns
        -------
        TreeNode
            A tree rooted by the outgroup.

        Raises
        ------
        TreeError
            Outgroup is not a proper subset of taxa in the tree.
        TreeError
            Outgroup is not monophyletic in the tree.

        Notes
        -----
        An outgroup is a subset of taxa that are usually distantly related from
        the remaining taxa (ingroup). The outgroup helps with locating the root
        of the ingroup, which are of interest in the study.

        This method reroots the tree at the lowest common ancestor (LCA) of the
        outgroup. By default, a new root will be placed at the midpoint between
        the LCA of outgroup and that of ingroup. But this behavior can be
        customized.

        This method requires the outgroup to be monophyletic, i.e., it forms a
        single clade in the tree. If the outgroup spans across the root of the
        tree, the method will reroot the tree within the ingroup such that the
        outgroup can form a clade in the rerooted tree, prior to rooting by
        outgroup.

        Tree caches (see :meth:`details <has_caches>`) will not be retained in the
        returned tree. In in-place mode, they will be cleared prior to rerooting. In
        copying mode, they will not be copied to the new tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(['((((a,b),(c,d)),(e,f)),g);'])
        >>> print(tree.ascii_art())
                                                /-a
                                      /--------|
                                     |          \-b
                            /--------|
                           |         |          /-c
                           |          \--------|
                  /--------|                    \-d
                 |         |
                 |         |          /-e
        ---------|          \--------|
                 |                    \-f
                 |
                  \-g

        >>> rooted = tree.root_by_outgroup(['a', 'b'])
        >>> print(rooted.ascii_art())
                            /-a
                  /--------|
                 |          \-b
                 |
        ---------|                    /-c
                 |          /--------|
                 |         |          \-d
                  \--------|
                           |                    /-e
                           |          /--------|
                            \--------|          \-f
                                     |
                                      \-g

        >>> rooted = tree.root_by_outgroup(['e', 'f', 'g'])
        >>> print(rooted.ascii_art())
                                      /-e
                            /--------|
                  /--------|          \-f
                 |         |
                 |          \-g
        ---------|
                 |                    /-c
                 |          /--------|
                 |         |          \-d
                  \--------|
                           |          /-b
                            \--------|
                                      \-a

        """
        outgroup = set(outgroup)

        if not outgroup < self.subset():
            raise TreeError("Outgroup is not a proper subset of taxa in the tree.")

        # locate the lowest common ancestor (LCA) of outgroup in the tree
        lca = self.lca(outgroup)

        # if LCA is root (i.e., outgroup is split across basal clades), root
        # the tree at a tip within the ingroup and locate LCA again
        if lca is self:
            for tip in self.tips():
                if tip.name not in outgroup:
                    tree = self.root_at(
                        tip, reset=reset, branch_attrs=branch_attrs, inplace=inplace
                    )
                    inplace = False  # no need to make copy again
                    break
            lca = tree.lca(outgroup)
        else:
            tree = self

        # test if outgroup is monophyletic
        if lca.count(tips=True) > len(outgroup):
            raise TreeError("Outgroup is not monophyletic in the tree.")

        # reroot the tree at LCA
        return tree.root_at(
            lca,
            above=above,
            reset=reset,
            branch_attrs=branch_attrs,
            root_name=root_name,
            inplace=inplace,
        )

    # ------------------------------------------------
    # Tree metrics
    # ------------------------------------------------

    def count(self, tips=False):
        r"""Get the count of nodes in the tree.

        Parameters
        ----------
        tips : bool, optional
            If True, only return the count of tips (default: False).

        Returns
        -------
        int
            The number of nodes.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> print(tree.count())
        9
        >>> print(tree.count(tips=True))
        5

        """
        if tips:
            return len(list(self.tips()))
        else:
            return len(list(self.traverse(include_self=True)))

    def subset(self, include_self=False):
        r"""Return a subset of taxa descending from self.

        A subset can be considered as taxa (tip names) within a clade defined by the
        current node (branch), selected from all taxa within the tree.

        Parameters
        ----------
        include_self : bool, optional
            Whether to include the current node if it is a tip (default: False).

            .. versionadded:: 0.6.3

        Returns
        -------
        frozenset of str
            The set of names at the tips of the clade that descends from self.

        See Also
        --------
        tips
        subsets
        bipart

        Notes
        -----
        This is a convenient method to return all taxa (tip names) rather than the tip
        nodes themselves. Internal node names will not be included.

        The returned value (a frozenset) is unordered and hashable, therefore can be
        used to define clades, lineages and taxon groups for efficient lookup. For
        example, one can check whether a taxon exists in the current tree or clade.

        By default, if this method is applied to a tip, an empty set will be returned,
        because a tip does not have descendants. If `include_self` is True, a single-
        element set containing the name of the tip will be returned. This behavior can
        be considered as returning taxa descending from the branch connecting self
        and its parent.

        Applying this method to the root node of a tree will return all taxa in the
        tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> print(tree.ascii_art())
                            /-a
                  /e-------|
                 |         |          /-b
                 |          \d-------|
        -i-------|                    \-c
                 |
                 |          /-f
                  \h-------|
                            \-g

        >>> sorted(tree.subset())
        ['a', 'b', 'c', 'f', 'g']

        >>> subset = tree.find('e').subset()
        >>> sorted(subset)
        ['a', 'b', 'c']

        >>> 'a' in subset
        True

        >>> 'f' in subset
        False

        """
        return frozenset({i.name for i in self.tips(include_self=include_self)})

    def subsets(
        self,
        within=None,
        include_full=False,
        include_tips=False,
        map_to_length=False,
    ):
        r"""Return all subsets of taxa defined by nodes descending from self.

        Parameters
        ----------
        within : iterable of str, optional
            A custom set of taxa to refine the result. Only taxa within it will be
            considered. If None (default), all taxa in the tree will be considered.

            .. versionadded:: 0.6.3

        include_full : bool, optional
            Whether to include a set of all taxa in the result. Default is False, as
            such a set provides no topological information.

            .. versionadded:: 0.6.3

        include_tips : bool, optional
            Whether to include subsets with only one taxon in the result. Default is
            False, as such sets provide no topological information.

            .. versionadded:: 0.6.3

        map_to_length : bool, optional
            If True, return a mapping of subsets to their branch lengths. Missing
            branch lengths will be replaced with 0. Default is False.

            .. versionadded:: 0.6.3

        Returns
        -------
        frozenset of frozenset of str, or
            All subsets of taxa defined by nodes descending from self. Returned if
            `map_to_length` is False.
        dict of {frozenset of str: float}
            Mapping of all subsets of taxa to their branch lengths. Returned if
            `map_to_length` is True.

        See Also
        --------
        subset
        compare_subsets
        biparts

        Notes
        -----
        The returned value represents the tree as a set of nested sets, each of which
        representing a clade in the tree. It is useful for assessing topological
        patterns of a tree.

        The returned value itself and each of its components (frozensets) are unordered
        and hashable, making it efficient for lookup and comparison. For example, one
        can check whether a group of taxa form a clade in the tree, regardless of its
        internal structure.

        This method can be applied to both rooted and unrooted trees. However, the
        underlying assumption is that the direction of descendance is from the current
        node to the tips below. That is, the root of the tree, even if not explicitly
        defined, should be at or above the current node. This should be considered when
        applying this method to an unrooted tree. If such an assumption is not present,
        one should consider using :meth:`biparts` instead.

        This method operates on the subtree below the current node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> print(tree.ascii_art())
                            /-a
                  /e-------|
                 |         |          /-b
                 |          \d-------|
        -i-------|                    \-c
                 |
                 |          /-f
                  \h-------|
                            \-g

        >>> subsets = tree.subsets()
        >>> for s in sorted(subsets, key=sorted):
        ...     print(sorted(s))
        ['a', 'b', 'c']
        ['b', 'c']
        ['f', 'g']

        >>> {'a', 'b', 'c'} in subsets
        True

        >>> {'a', 'b'} in subsets
        False

        """
        if not (getall := within is None):
            if not isinstance(within, (set, frozenset, dict)):
                within = frozenset(within)

        # initiate result
        subsets = []
        subsets_append = subsets.append
        if map_to_length:
            lengths = []
            lengths_append = lengths.append

        # If the current subset has been encountered during postorder traversal, it
        # must be the immediately previous subset. This happens when a single-child
        # node is encountered after refining taxa to the "within" set.
        last = None

        for node in self.postorder(include_self=True):
            # tip: create a one-taxon set
            if not node.children:
                if getall or node.name in within:
                    subset = frozenset([node.name])
                else:
                    subset = frozenset()

            # internal node: merge sets of children
            else:
                subset = frozenset()
                for child in node.children:
                    subset |= child._subset
                    del child._subset

            # add to result
            if subset and (include_tips or len(subset) > 1):
                if subset != last:
                    subsets_append(last := subset)
                    if map_to_length:
                        lengths_append(node.length or 0.0)
                elif map_to_length:
                    lengths[-1] += node.length or 0.0

            node._subset = subset

        # final clean up
        del self._subset

        # remove the full set
        if not include_full:
            subsets = subsets[:-1]
            if map_to_length:
                lengths = lengths[:-1]

        if map_to_length:
            return dict(zip(subsets, lengths))
        else:
            return frozenset(subsets)

    def bipart(self):
        r"""Return a bipartition of the tree at the current branch.

        .. versionadded:: 0.6.3

        A bipartition, partition or split of a tree is the division of all taxa (tip
        names) into two complementary subsets, separated at a given branch. In this
        context, it is the branch connecting self and its parent. One subset consists
        of all taxa descending from self and the other consists of all remaining taxa.
        The smaller subset of the two is returned.

        Returns
        -------
        frozenset of str
            The set of names at the tips on the smaller side of the current branch.

        See Also
        --------
        subset
        biparts

        Notes
        -----
        A bipartition describes the topological placement of a branch regardless of
        other branches and the root of the tree.

        The returned value is a set of tip names on the smaller side of the branch, as
        determined by the number of tips. If a tie is observed, the tip names on both
        sides are sorted lexicographically and the first set is returned.

        The returned value (a frozenset) is unordered and hashable, making it efficient
        for lookup and comparison. For example, one can check whether two branches in
        two unrooted trees with the same taxa agree with each other.

        Rerooting a tree will not change the bipartition of a branch. However, one
        should be cautious because this method applies to a node, and rerooting may
        change the branch above the current node.

        Applying this method to a root node will return an empty set. Applying this
        method to a tip will return a single-element set containing the tip name. These
        two situations produce outputs independent of the topology of the tree.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["(((a,(b,c)X)Y,d)Z,(e,f),g);"])
        >>> print(tree.ascii_art())
                                      /-a
                            /Y-------|
                           |         |          /-b
                  /Z-------|          \X-------|
                 |         |                    \-c
                 |         |
                 |          \-d
        ---------|
                 |          /-e
                 |---------|
                 |          \-f
                 |
                  \-g

        Clade has less than half taxa, return them.

        >>> sorted(tree.find('X').bipart())
        ['b', 'c']

        Clade has more than half taxa, return remaining taxa.

        >>> sorted(tree.find('Z').bipart())
        ['e', 'f', 'g']

        Clade has exactly half taxa, return the lexicographically smaller side.

        >>> sorted(tree.find('Y').bipart())
        ['a', 'b', 'c']

        A second tree with the same topology but different root position.

        >>> tree2 = TreeNode.read(["((c,b)X2,a,(((f,e),g)Y2,d));"])
        >>> print(tree2.ascii_art())
                            /-c
                  /X2------|
                 |          \-b
                 |
                 |--a
        ---------|
                 |                              /-f
                 |                    /--------|
                 |          /Y2------|          \-e
                 |         |         |
                  \--------|          \-g
                           |
                            \-d

        Although the tree has been re-positioned, the corresponding branches have the
        same bipartitions, whereas non-corresponding branches don't.

        >>> tree.find('X').bipart() == tree2.find('X2').bipart()
        True

        >>> tree.find('Y').bipart() == tree2.find('Y2').bipart()
        False

        """
        bipart = self.subset(include_self=True)
        full = self.root().subset(include_self=True)
        if (size := len(bipart)) > (th := len(full) * 0.5):
            bipart = full - bipart
        elif size == th:
            # sort the elements of each part by lexicographic order, then order the two
            # parts and pick the smaller part
            bipart, _ = sorted([bipart, full - bipart], key=sorted)
        return bipart

    def biparts(self, within=None, include_tips=False, map_to_length=False, full=None):
        r"""Return all bipartitions within the tree under self.

        .. versionadded:: 0.6.3

        Parameters
        ----------
        within : iterable of str, optional
            A custom set of taxa to refine the result. Only taxa within it will be
            considered. If None (default), all taxa in the tree will be considered.
        include_tips : bool, optional
            Whether to include bipartitions with only one taxon at either side.
            Default is False, as such bipartitions provide no topological
            information.
        map_to_length : bool, optional
            If True, return a mapping of subsets to their branch lengths. Missing
            branch lengths will be replaced with 0. Default is False.
        full : frozenset of str, optional
            Pre-computed full set of taxa of the current tree. Providing this parameter
            can save one tree traversal from computing.

        Returns
        -------
        frozenset of frozenset of str, or
            All sets of names at the tips on the smaller side of each branch. Returned
            if `map_to_length` is False.
        dict of {frozenset of str: float}
            Mapping of All sets of smaller-side tip names to branch lengths. Returned
            if `map_to_length` is True.

        See Also
        --------
        bipart
        subsets

        Notes
        -----
        The returned value represents the tree as a set of nested sets, each of which
        representing the position of a branch in the tree. It is useful for assessing
        topological patterns of a tree.

        The returned value itself and each of its components (frozensets) are unordered
        and hashable, making it efficient for lookup and comparison. For example, one
        can check whether the topologies of two trees are consistent, regardless of
        their root positions.

        This method can be applied to both rooted and unrooted trees. However, a rooted
        tree implies the direction of descendance, which may violate the purpose of
        bipartitioning a tree on arbitrary branches. If this is a concern, one should
        consider using :meth:`subsets` instead.

        This method operates on the subtree below the current node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)),(d,e),f);"])
        >>> print(tree.ascii_art())
                            /-a
                  /--------|
                 |         |          /-b
                 |          \--------|
                 |                    \-c
        ---------|
                 |          /-d
                 |---------|
                 |          \-e
                 |
                  \-f

        Return all bipartitions of an unrooted tree.

        >>> biparts = tree.biparts()
        >>> for s in sorted(biparts, key=sorted):
        ...     print(sorted(s))
        ['a', 'b', 'c']
        ['b', 'c']
        ['d', 'e']

        A second tree with the same topology but different root position.

        >>> tree2 = TreeNode.read(["(a,((b,c),((d,e),f)));"])
        >>> print(tree2.ascii_art())
                  /-a
                 |
        ---------|                    /-b
                 |          /--------|
                 |         |          \-c
                  \--------|
                           |                    /-d
                           |          /--------|
                            \--------|          \-e
                                     |
                                      \-f

        Although the tree has been re-positioned, the bipartitions remain the same.

        >>> biparts == tree2.biparts()
        True

        """
        # identify full set (universe)
        if full is None:
            full = self.subset()
        if not (getall := within is None):
            if not isinstance(within, (set, frozenset)):
                within = frozenset(within)
            full &= within
        th = len(full) * 0.5

        # initiate result
        if map_to_length:
            biparts = {}
            biparts_get = biparts.get
        else:
            biparts = []
            biparts_append = biparts.append

        for node in self.postorder(include_self=False):
            # tip: create a one-taxon set
            if not node.children:
                if getall or node.name in full:
                    bipart = frozenset([node.name])
                else:
                    bipart = frozenset()
                flip = False

            # internal node: merge sets of children
            # `_bipart` of a node is either the taxa below it, or, if the former has
            # reached half of the full set, it "flips" to the other half that is above
            # the node, and `_flip` will be set to True.
            # Taxa below should be united, whereas taxa above should be intersected.
            # If at least one child is already flipped, the current node should also be
            # flipped. Otherwise, the set will be compared to the half to determine the
            # flipping status.
            else:
                aboves, belows = [], []
                for child in node.children:
                    if child._flip:
                        aboves.append(child._bipart)
                    else:
                        belows.append(child._bipart)
                    del child._bipart
                    del child._flip

                if aboves:
                    bipart = frozenset.intersection(*aboves).difference(*belows)
                    flip = True
                else:
                    bipart = frozenset().union(*belows)
                    flip = False
                    if (size := len(bipart)) >= th:
                        other = full - bipart
                        if size > th or sorted(bipart) > sorted(other):
                            bipart = other
                            flip = True

            # add to result
            if bipart and (include_tips or len(bipart) > 1):
                if map_to_length:
                    biparts[bipart] = biparts_get(bipart, 0.0) + (node.length or 0.0)
                else:
                    biparts_append(bipart)

            node._bipart = bipart
            node._flip = flip

        # final clean up
        for child in self.children:
            del child._bipart
            del child._flip

        return biparts if map_to_length else frozenset(biparts)

    def _extract_support(self):
        """Extract the support value from a node label, if available.

        Returns
        -------
        tuple of
            int, float or None
                The support value extracted from the node label.
            str or None
                The node label with the support value stripped.

        """
        support, label = None, None
        if self.name:
            # separate support value from node name by the first colon
            left, _, right = self.name.partition(":")
            try:
                support = int(left)
            except ValueError:
                try:
                    support = float(left)
                except ValueError:
                    pass
            # strip support value from node name
            label = right or None if support is not None else self.name
        return support, label

    def _node_label(self):
        """Generate a node label.

        The label will be in the format of "support:name" if both exist,
        or "support" or "name" if either exists.

        Returns
        -------
        str
            Generated node label.

        """
        lblst = []
        if self.support is not None:  # prevents support of NoneType
            lblst.append(str(self.support))
        if self.name:  # prevents name of NoneType
            lblst.append(self.name)
        return ":".join(lblst)

    def assign_supports(self):
        """Extract support values from internal node labels of a tree.

        Notes
        -----
        A "support value" measures the confidence or frequency of the incoming
        branch (the branch from parent to self) of an internal node in a tree.
        Roots and tips do not have support values. To extract a support value
        from a node label, this method reads from left and stops at the first
        ":" (if any), and attempts to convert it to a number.

        For examples: "(a,b)1.0", "(a,b)1.0:2.5", and "(a,b)'1.0:species_A'".
        In these cases the support values are all 1.0.

        For examples: "(a,b):1.0" and "(a,b)species_A". In these cases there
        are no support values.

        If a support value is successfully extracted, it will be stripped from
        the node label and assigned to the `support` property.

        .. note::
            Mathematically, "support value" is a property of a branch, not a
            node, although they are usually attached to nodes in tree file
            formats [1]_.

        References
        ----------
        .. [1] Czech, Lucas, Jaime Huerta-Cepas, and Alexandros Stamatakis. "A
           Critical Review on the Use of Support Values in Tree Viewers and
           Bioinformatics Toolkits." Molecular biology and evolution 34.6
           (2017): 1535-1542.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> newick = "((a,b)95,(c,d):1.1,(e,f)'80:speciesA':1.0);"
        >>> tree = TreeNode.read([newick])
        >>> tree.assign_supports()
        >>> tree.lca(['a', 'b']).support
        95
        >>> tree.lca(['c', 'd']).support is None
        True
        >>> tree.lca(['e', 'f']).support
        80
        >>> tree.lca(['e', 'f']).name
        'speciesA'

        """
        for node in self.traverse():
            if node.is_root() or node.is_tip():
                node.support = None
            else:
                node.support, node.name = node._extract_support()

    def is_bifurcating(self, strict=False, include_self=True):
        r"""Check if the tree is bifurcating.

        .. versionadded:: 0.6.3

        Parameters
        ----------
        strict : bool, optional
            Whether to consider single-child nodes as violations of bifurcation.
            Default is False.
        include_self : bool, optional
            If False, will not check the current node. This is useful for checking an
            unrooted tree, whose root node may have three children. Default is True.

        See Also
        --------
        bifurcate
        prune

        Notes
        -----
        In a bifurcating tree (a.k.a. binary tree), every node has at most two
        children. The property of bifurcation is necessary for a wide range of tree
        analyses. In contrast, if a node has three or more children, it is considered
        as multifurcating, or polytomy in phylogenetics.

        In strict mode, every internal node (including root) has to have exactly two
        children in order for the tree to be bifurcating. Single-child nodes are
        considered as violations. These nodes can be collapsed by :meth:`prune`.

        This method operates on the subtree below the current node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b,c),(d,e))root;"])
        >>> tree.is_bifurcating()
        False

        """
        test = ne if strict else gt
        for node in self.traverse(include_self=include_self):
            if (children := node.children) and test(len(children), 2):
                return False
        return True

    def observed_node_counts(self, tip_counts):
        """Return counts of node observations from counts of tip observations.

        Parameters
        ----------
        tip_counts : dict of ints
            Counts of observations of tips. Keys correspond to tip names in
            ``self``, and counts are unsigned ints.

        Returns
        -------
        dict
            Counts of observations of nodes. Keys correspond to node names
            (internal nodes or tips), and counts are unsigned ints.

        Raises
        ------
        ValueError
            If a count less than one is observed.
        MissingNodeError
            If a count is provided for a tip not in the tree, or for an
            internal node.

        """
        result = defaultdict(int)
        for tip_name, count in tip_counts.items():
            if count < 1:
                raise ValueError("All tip counts must be greater than zero.")
            else:
                t = self.find(tip_name)
                if not t.is_tip():
                    raise MissingNodeError(
                        "Counts can only be for tips in the tree. %s is an "
                        "internal node." % t.name
                    )
                result[t] += count
                for internal_node in t.ancestors():
                    result[internal_node] += count
        return result

    @aliased("accumulate_to_ancestor", "0.6.3")
    def depth(
        self, ancestor=None, include_root=False, use_length=True, missing_as_zero=False
    ):
        r"""Calculate the depth of the current node.

        The **depth** of a node is the sum of branch lengths from it to the root of the
        tree.

        Parameters
        ----------
        ancestor : TreeNode, optional
            An ancestral node of self. If provided, the distance from self to this node
            instead of the root node will be calculated.

            .. versionchanged:: 0.6.3
                Becomes optional.

        include_root : bool, optional
            If True, the distance will include the length of the root node, or the
            given ancestral node if ``ancestor`` is provided. Default is False.

            .. versionadded:: 0.6.3

        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the number
            of branches (False) from self to root.

            .. versionadded:: 0.6.3

        missing_as_zero : bool, optional
            When a node without an associated branch length is encountered, raise an
            error (False, default) or use 0 (True). Applicable when ``use_length`` is
            True.

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The depth of self.

        Raises
        ------
        NoParentError
            If the given ancestral node is not an ancestor of self.
        NoLengthError
            If nodes without branch length are encountered, but ``missing_as_zero`` is
            False.

        See Also
        --------
        height
        distance

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        >>> tree.find('a').depth()
        4.0
        >>> tree.find('a').depth(tree.find('c'))
        1.0

        """
        curr = self
        path = [curr]
        path_append = path.append
        if ancestor is None:
            while (curr := curr.parent) is not None:
                path_append(curr)
        else:
            try:
                while curr is not ancestor:
                    path_append(curr := curr.parent)
            except AttributeError:
                raise NoParentError("Provided ancestor is not ancestral to self.")
        if not include_root:
            path = path[:-1]
        if not use_length:
            return float(len(path))
        if missing_as_zero:
            return sum(x.length or 0.0 for x in path)
        try:
            return sum(x.length for x in path)
        except TypeError:
            raise NoLengthError("Nodes without branch length are encountered.")

    def height(self, include_self=False, use_length=True, missing_as_zero=False):
        r"""Calculate the height of the current node.

        .. versionadded:: 0.6.3

        The **height** of a node is the maximum sum of branch lengths from it to any of
        its descending tips.

        Parameters
        ----------
        include_self : bool, optional
            If True, the height will include the length of the current node. Default
            is False.
        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the number
            of branches (False) from self to the most distant tip.
        missing_as_zero : bool, optional
            When a node without an associated branch length is encountered, raise an
            error (False, default) or use 0 (True). Applicable when ``use_length`` is
            True.

        Returns
        -------
        float
            The height of self.
        TreeNode
            The most distant descending tip from self.

        Raises
        ------
        NoLengthError
            If nodes without branch length are encountered, but ``missing_as_zero`` is
            False.

        See Also
        --------
        depth
        distance

        Notes
        -----
        When a tie is observed among multiple tips, only one of them will be returned.
        The choice is stable. This often happens when ``use_length=False``.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        >>> dist, tip = tree.find('c').height()
        >>> dist
        2.0
        >>> tip.name
        'b'

        """
        errmsg = "Nodes without branch length are encountered."
        maxkey = itemgetter(0)
        for node in self.postorder(include_self=True):
            if not node.children:
                node._height = (0.0, node)
            else:
                heights = []
                for child in node.children:
                    H, tip = child._height
                    del child._height
                    if not use_length:
                        H += 1.0
                    elif (L := child.length) is not None:
                        H += L
                    elif not missing_as_zero:
                        raise NoLengthError(errmsg)
                    heights.append((H, tip))
                node._height = max(heights, key=maxkey)
        H, tip = self._height
        del self._height
        if include_self:
            if not use_length:
                H += 1.0
            elif (L := self.length) is not None:
                H += L
            elif not missing_as_zero:
                raise NoLengthError(errmsg)
        return H, tip

    @aliased("descending_branch_length", "0.6.3")
    @params_aliased([("nodes", "tip_subset", "0.6.3", True)])
    def total_length(self, nodes=None, include_stem=False, include_self=False):
        r"""Calculate the total length of branches descending from self.

        Parameters
        ----------
        nodes : iterable of TreeNode or str, optional
            Instances or names of a subset of descending nodes to refine the result.
            If provided, the total length of branches connecting these nodes will be
            returned. Otherwise, the total branch length of the tree will be returned.

            .. versionchanged:: 0.6.3
                Can accept TreeNode instances in addition to names.
                Can accept internal nodes in addition to tips.

        include_stem : bool, optional
            Whether to include the path from the lowest common ancestor (LCA) of the
            subset of nodes to self. Applicable when ``nodes`` is specified. Default is
            False.

            .. versionadded:: 0.6.3

        include_self : bool, optional
            Whether to include the length of self. When ``nodes`` is provided and
            ``include_stem`` is False, it is instead the LCA of the subset of nodes.
            Default is False.

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The total descending branch length.

        Raises
        ------
        MissingNodeError
            If some nodes are not found in the tree or are not descendants of self.

        Notes
        -----
        The metric can be considered as the total amount of evolutionary change across
        all lineages in the tree.

        This metric is closely related to phylogenetic diversity (PD) in community
        ecology. When ``include_stem`` is True, it is equivalent to Faith's PD (see
        :func:`~skbio.diversity.alpha.faith_pd`). However, this method is optimized
        to handle a single set of nodes, whereas the referred function is optimized
        to simultaneously calculate for multiple taxon sets (i.e., communities).

        Missing branch lengths will be replaced with 0.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read([
        ...     "(((A:.1,B:1.2)C:.6,(D:.9,E:.6)F:.9)G:2.4,(H:.4,I:.5)J:1.3)K;"])
        >>> print(tree.ascii_art())
                                      /-A
                            /C-------|
                           |          \-B
                  /G-------|
                 |         |          /-D
                 |          \F-------|
        -K-------|                    \-E
                 |
                 |          /-H
                  \J-------|
                            \-I

        Calculate the total branch length of the tree.

        >>> L = tree.total_length()
        >>> print(round(L, 1))
        8.9

        Calculate the total branch length connecting three taxa.

        >>> L = tree.total_length(['A', 'E', 'H'])
        >>> print(round(L, 1))
        6.3

        """
        ## shortcut for the entire subtree
        if not nodes:
            return sum(
                n.length or 0.0 for n in self.postorder(include_self=include_self)
            )

        nodes = [self.find(x) for x in nodes]

        # Identify all nodes that need to be visited during the navigation from all
        # tips to the root. This algorithm resembles that of `lca`. However, we will
        # separate the visited nodes of the first path and all other paths. Also, we
        # don't need to record the previous node. All we need is whether each node is
        # unique in all paths.
        first_path = []
        first_path_append = first_path.append
        curr = next(nodes := iter(nodes))
        while curr is not None:
            first_path_append(curr)
            curr._unique = True
            curr = curr.parent

        other_paths = []
        other_paths_append = other_paths.append
        for curr in nodes:
            while not hasattr(curr, "_unique"):
                other_paths_append(curr)
                curr._unique = True
                curr = curr.parent
            curr._unique = False

        # Iterate the first path in reverse order (from root to starting node) and find
        # the indices of self and LCA.
        i_self, i_lca = None, 0
        for i in reversed(range(len(first_path))):
            if (node := first_path[i]) is self:
                i_self = i
            if node._unique is False:
                i_lca = i
                break

        # clear temporary attribute
        for node in first_path:
            del node._unique
        for node in other_paths:
            del node._unique

        # If all nodes are descendants of self, LCA must also be self or one of its
        # descendants, and self must be identified when iterating the first path.
        if i_self is None:
            raise MissingNodeError("Some nodes are not descendants of self.")

        # Identify the range of nodes to be included in calculation depending on the
        # parameter setting
        stop = (i_self if include_stem else i_lca) + include_self

        # sum up branch lengths
        return (
            sum(n.length or 0.0 for n in chain(first_path[:stop], other_paths)) or 0.0
        )

    def distance(self, other, use_length=True, missing_as_zero=False):
        r"""Calculate the distance between self and another node.

        Parameters
        ----------
        other : TreeNode
            The node to compute a distance to.
        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the number
            of branches (False) connecting self and other.

            .. versionadded:: 0.6.3

        missing_as_zero : bool, optional
            When a node without an associated branch length is encountered, raise an
            error (False, default) or use 0 (True). Applicable when ``use_length`` is
            True.

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The distance between two nodes.

        Raises
        ------
        NoLengthError
            If nodes without branch length are encountered, but ``missing_as_zero`` is
            False.

        See Also
        --------
        path
        cophenet
        compare_cophenet
        maxdist

        Notes
        -----
        The distance between two nodes is the length of the path (branches) connecting
        them. It is also known as the patristic distance [1]_.

        When ``use_length=False``, it is the number of branches in the path.

        This method can be used to compute the distance between two given nodes.
        However, it is not optimized for computing all pairwise tip distances. Use
        :meth:`cophenet` instead for that purpose.

        References
        ----------
        .. [1] Fourment, M., & Gibbs, M. J. (2006). PATRISTIC: a program for
           calculating patristic distances and graphically comparing the components of
           genetic change. BMC evolutionary biology, 6, 1-5.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        >>> tip_a = tree.find('a')
        >>> tip_d = tree.find('d')
        >>> tip_a.distance(tip_d)
        14.0
        >>> tip_a.distance(tip_d, use_length=False)
        4.0

        """
        _, self_path, other_path = self._path(other)
        if not use_length:
            return float(len(self_path) + len(other_path))
        if missing_as_zero:
            return sum(x.length or 0.0 for x in chain(self_path, other_path))
        try:
            return sum(x.length for x in chain(self_path, other_path))
        except TypeError:
            raise NoLengthError("Nodes without branch length are encountered.")

    @aliased("get_max_distance", "0.6.3")
    def maxdist(self, use_length=True):
        r"""Return the maximum distance between any pair of tips in the tree.

        This measure is also referred to as the **diameter** of a tree.

        Parameters
        ----------
        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the number
            of branches (False) connecting each pair of tips.

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The distance between the two most distant tips in the tree.
        tuple of (TreeNode, TreeNode)
            The two most distant tips in the tree.

        See Also
        --------
        distance
        cophenet
        scipy.cluster.hierarchy.maxdists

        Notes
        -----
        This method identifies the two furthest apart tips in a tree, as measured by
        the sum of branch lengths (i.e., patristic distance) connecting them. Missing
        branch lengths will be replaced with 0. When ``use_length=False``, the number
        of branches connecting two tips will be considered instead.

        When a tie is observed among more than one pair of tips, only one pair will be
        returned. The choice is stable. This often happens when ``use_length=False``.

        This method operates on the subtree below the current node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        >>> dist, tips = tree.maxdist()
        >>> dist
        16.0
        >>> [n.name for n in tips]
        ['e', 'b']

        """
        # The code performs a post-order traversal and appends two pieces of
        # information to each node:
        #   a: The maximum distance from the node to any descending tip.
        #   b: The maximum distance between any two descending tips.
        # The information is updated at each internal node:
        #   a becomes the maximum of any child's (a + length).
        #   b becomes the larger of the maximum of any b and the sum of the two
        # largest (a + length). The latter represents the new plausible maximum
        # distance that crosses the node.
        maxkey = itemgetter(0)
        for node in self.postorder():
            # initialize maximum at tip: (up_dist, up_tip, in_dist, in_tip1, in_tip2)
            if not node.children:
                node._maxdist = (0, node, 0, None, None)

            # internal node: update the maximum
            elif len(children := node.children) > 1:
                ups, ins = [], []
                for child in children:
                    up_dist, up_tip, in_dist, in_tip1, in_tip2 = child._maxdist
                    del child._maxdist
                    ups.append(
                        (up_dist + (child.length or 0.0 if use_length else 1), up_tip)
                    )
                    ins.append((in_dist, in_tip1, in_tip2))

                # compare the previous maximum with the distance between the two
                # longest descendants from any two child clades
                ups.sort(key=maxkey, reverse=True)
                (up_dist, up_tip), (up_dist2, up_tip2) = ups[:2]
                in_dist, in_tip1, in_tip2 = max(ins, key=maxkey)
                if (x_dist := up_dist + up_dist2) > in_dist:
                    node._maxdist = (up_dist, up_tip, x_dist, up_tip, up_tip2)
                else:
                    node._maxdist = (up_dist, up_tip, in_dist, in_tip1, in_tip2)

            # internal node with only one child: inherit the maximum
            else:
                (child,) = children
                up_dist, up_tip, in_dist, in_tip1, in_tip2 = child._maxdist
                del child._maxdist
                up_dist += child.length or 0.0 if use_length else 1
                node._maxdist = (up_dist, up_tip, in_dist, in_tip1, in_tip2)

        max_dist, max_tip1, max_tip2 = self._maxdist[2:]
        del self._maxdist
        if not use_length:
            max_dist = float(max_dist)
        return max_dist, (max_tip1, max_tip2)

    @aliased("tip_tip_distances", "0.6.3")
    def cophenet(self, endpoints=None, use_length=True):
        r"""Return a distance matrix between each pair of tips in the tree.

        Parameters
        ----------
        endpoints : list of TreeNode or str, optional
            Tips or their names (i.e., taxa) to be included in the calculation. The
            returned distance matrix will use this order. If not specified, all tips
            will be included.
        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the number
            of branches (False) connecting each pair of tips.

            .. versionadded:: 0.6.3

        Returns
        -------
        DistanceMatrix
            The cophenetic distance matrix.

        Raises
        ------
        MissingNodeError
            If any of the specified ``endpoints`` are not found in the tree.
        DuplicateNodeError
            If the specified ``endpoints`` have duplicates.
        ValueError
            If any of the specified ``endpoints`` are not tips.

        See Also
        --------
        distance
        compare_cophenet
        scipy.cluster.hierarchy.cophenet

        Notes
        -----
        The cophenetic distance [1]_ between a pair of tips is essentially the sum of
        branch lengths connecting them (i.e., patristic distance [2]_, see
        :meth:`distance`). It measures the divergence between two taxa in evolution.

        This method calculates the cophenetic distances between all pairs of tips in a
        tree and returns a distance matrix. Missing branch lengths will be replaced with
        0's. If ``use_length`` is False, the method instead calculates the number of
        branches connecting each pair of tips. This method operates on the subtree below
        the current node.

        In hierarchical clustering, the cophenetic distance is commonly used to measure
        the dissimilarity between two objects before they are joined in a dendrogram.
        In that context, it is also defined as the height of the lowest common ancestor
        (LCA) from the surface of the tree. However, phylogenetic trees are usually
        non-ultrametric (e.g., :func:`~skbio.tree.nj`), and the two child clades of a
        node may have different heights. Therefore, the cophenetic distance is instead
        defined as the patristic distance between the two tips. For ultrametric trees
        (e.g., :func:`~skbio.tree.upgma`), this method's result should match SciPy's
        :func:`~scipy.cluster.hierarchy.cophenet`.

        One should also distinguish cophenetic distance from a related metric:
        cophenetic value [1]_, which is the patristic distance between the LCA of two
        tips and the root of the tree. It quantifies the shared evolutionary history
        between two taxa, as in contrast to the cophenetic distance.

        References
        ----------
        .. [1] Sokal, R. R., & Rohlf, F. J. (1962). The comparison of dendrograms by
           objective methods. Taxon, 33-40.

        .. [2] Fourment, M., & Gibbs, M. J. (2006). PATRISTIC: a program for
           calculating patristic distances and graphically comparing the components of
           genetic change. BMC evolutionary biology, 6, 1-5.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])

        Calculate cophenetic distances as the sum of branch lengths (i.e., patristic
        distance).

        >>> mat = tree.cophenet()
        >>> print(mat)
        4x4 distance matrix
        IDs:
        'a', 'b', 'd', 'e'
        Data:
        [[  0.   3.  14.  15.]
         [  3.   0.  15.  16.]
         [ 14.  15.   0.   9.]
         [ 15.  16.   9.   0.]]

        Calculate cophenetic distances as the number of branches.

        >>> mat = tree.cophenet(use_length=False)
        >>> print(mat)
        4x4 distance matrix
        IDs:
        'a', 'b', 'd', 'e'
        Data:
        [[ 0.  2.  4.  4.]
         [ 2.  0.  4.  4.]
         [ 4.  4.  0.  2.]
         [ 4.  4.  2.  0.]]

        """
        taxa = []
        taxa_append = taxa.append

        # Include all tips.
        # `tips()` performs a postorder traversal, which guarantees the continuity
        # of tip indices within each node. A `_range` attribute is assigned to each
        # node, representing the range of tip indices.
        if not endpoints:
            for i, tip in enumerate(self.tips()):
                tip._range = (i, i + 1)
                taxa_append(tip.name)
            num_tips = len(taxa)

            # A tree could have duplicate taxa so this check is desired.
            if len(set(taxa)) < num_tips:
                raise DuplicateNodeError(f"Tree contains duplicate tip names.")

        # Include only selected tips in order.
        # Only selected tips are indexed, but the continuity of tip indices (see above)
        # is still ensured.
        else:
            idxmap = {}
            for i, tip in enumerate(endpoints):
                # The `find` call will raise if there are duplicate taxa in the tree.
                tip = self.find(tip)
                if tip.children:
                    raise ValueError(f"Node with name '{tip.name}' is not a tip.")
                taxa_append(name := tip.name)
                if name in idxmap:
                    raise DuplicateNodeError(f"Duplicate tip name '{name}' found.")
                idxmap[name] = i
            num_tips = len(taxa)

            # Create an index array to store the order of indices of original tips.
            order = np.empty(num_tips, dtype=int)
            i = 0
            for tip in self.tips():
                if (name := tip.name) in idxmap:
                    tip._range = (i, i + 1)
                    order[idxmap[name]] = i
                    i += 1

        # Initiate the resulting distance matrix.
        result = np.zeros((num_tips, num_tips))

        # An intermediate vector storing the accumulative distance from each tip to
        # the current node.
        depths = np.zeros(num_tips)

        # Traverse internal nodes.
        # This method involves two postorder traversals. Theoretically, one can perform
        # only one traversal, and store tip and internal node references into two lists
        # for use. However, this method isn't more efficient according to benchmarks.
        for node in self.postorder():
            if not node.children:
                continue

            # Record tip ranges of each child clade, and increment the tip depths.
            ranges = []
            for child in node.children:
                if not hasattr(child, "_range"):
                    continue
                ranges.append(range_ := slice(*child._range))
                depths[range_] += (child.length or 0.0) if use_length else 1
                del child._range

            # Calculate tip-to-tip distances between each pair of child clades, and
            # save the results to both upper and lower triangles of the resulting
            # distance matrix.
            # This is significantly faster than saving to only one triangle and doing
            # doing `result += result.T` after the iteration.
            for range1, range2 in combinations(ranges, 2):
                dists = depths[range1][:, np.newaxis] + depths[range2]
                result[range1, range2] = dists
                result[range2, range1] = dists.T

            # Due to the continuity of tip indices (see above), it is guaranteed that
            # the first child is the smallest and the last child is the largest.
            if ranges:
                node._range = (ranges[0].start, ranges[-1].stop)

        if hasattr(self, "_range"):
            del self._range

        # Reorder the distance matrix to reflect the given order of endpoints.
        if endpoints:
            result = result[order][:, order]

        # Skip validation as all items to validate are guaranteed.
        return DistanceMatrix(result, taxa, validate=False)

    @params_aliased([("shared_only", "exclude_absent_taxa", "0.6.3", True)])
    def compare_subsets(self, other, shared_only=False, proportion=True):
        r"""Calculate the difference of subsets between two trees.

        Parameters
        ----------
        other : TreeNode
            The other tree to compare with.
        shared_only : bool, optional
            Only consider taxa shared with the other tree. Default is False.
        proportion : bool, optional
            Whether to return count (False) or proportion (True, default) of different
            subsets.

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The count or proportion of subsets that differ between the trees.

        See Also
        --------
        subsets
        compare_rfd
        compare_biparts

        Notes
        -----
        This metric is equivalent to the Robinson-Foulds distance on rooted trees.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree1 = TreeNode.read(["((a,b),(c,d));"])
        >>> tree2 = TreeNode.read(["(((a,b),c),d);"])
        >>> tree1.compare_subsets(tree2)
        0.5

        """
        return _topo_dists((self, other), True, shared_only, proportion)[0]

    def compare_biparts(self, other, proportion=True):
        r"""Calculate the difference of bipartitions between two trees.

        .. versionadded:: 0.6.3

        Parameters
        ----------
        other : TreeNode
            The other tree to compare with.
        proportion : bool, optional
            Whether to return count (False) or proportion (True, default) of different
            bipartitions.

        Returns
        -------
        float
            The count or proportion of bipartitions that differ between the trees.

        See Also
        --------
        biparts
        compare_rfd
        compare_subsets

        Notes
        -----
        This metric is equivalent to the Robinson-Foulds distance on unrooted trees.

        Only taxa shared between the two trees are considered.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree1 = TreeNode.read(["((a,b),(c,d));"])
        >>> tree2 = TreeNode.read(["(((a,b),c),d);"])
        >>> tree1.compare_biparts(tree2)
        0.0

        """
        return _topo_dists((self, other), False, True, proportion)[0]

    def compare_rfd(self, other, proportion=False, rooted=None):
        r"""Calculate Robinson-Foulds distance between two trees.

        Parameters
        ----------
        other : TreeNode
            The other tree to compare with.
        proportion : bool, optional
            Whether to return the RF distance as count (False, default) or proportion
            (True).
        rooted : bool, optional
            Whether to consider the trees as rooted or unrooted. If None (default),
            this will be determined based on whether self is rooted. However, one
            can override it by explicitly specifying True (rooted) or False (unrooted).

            .. versionadded:: 0.6.3

        Returns
        -------
        float
            The Robinson-Foulds distance as count or proportion between the trees.

        .. versionchanged:: 0.6.3
            When the tree is unrooted, the calculation is based on bipartitions instead
            of subsets.

        Notes
        -----
        The Robinson-Foulds (RF) distance, a.k.a. symmetric difference, is a measure of
        topological dissimilarity between two trees. It was originally described in
        [1]_. It is calculated as the number of bipartitions that differ between two
        unrooted trees. It is equivalent to :meth:`compare_biparts`.

        .. math::

           \text{RF}(T_1, T_2) = |S_1 \triangle S_2| = |(S_1 \setminus S_2) \cup (S_2
           \setminus S_1)|

        where :math:`S_1` and :math:`S_2` are the sets of bipartitions of trees
        :math:`T_1` and :math:`T_2`, respectively.

        For rooted trees, the RF distance is calculated as the number of unshared
        clades (subsets of taxa) [2]_. It is equivalent to :meth:`compare_subsets`.

        This method automatically determines whether to use the unrooted or rooted RF
        distance based on whether self is rooted or not. Specifically, if self has two
        children (see :meth:`details <unroot>`), or has a parent (i.e., it is a subtree
        within a larger tree), it will be considered as rooted. Otherwise it will be
        considered as unrooted.

        One can override this automatic decision by setting the ``rooted`` parameter,
        which is recommended for explicity.

        By specifying ``proportion=True``, a unit distance will be returned, ranging
        from 0 (identical) to 1 (completely different).

        This method operates on the subtrees below the given nodes. Only taxa shared
        between the two trees are considered. Taxa unique to either tree are excluded
        from the calculation.

        See Also
        --------
        compare_wrfd
        compare_subsets
        compare_biparts
        skbio.tree.rf_dists

        References
        ----------
        .. [1] Robinson, D. F., & Foulds, L. R. (1981). Comparison of phylogenetic
           trees. Mathematical biosciences, 53(1-2), 131-147.

        .. [2] Bogdanowicz, D., & Giaro, K. (2013). On a matching distance between
           rooted phylogenetic trees. International Journal of Applied Mathematics
           and Computer Science, 23(3), 669-684.

        Examples
        --------
        Calculate the RF distance between two unrooted trees with the same taxa but
        different topologies. Each tree has three non-trivial bipartitions, as defined
        by individual internal branches, among which one pair (abc|def) is shared
        whereas the other two of each tree are unique (ab|cdef, abcf|de, bc|adef,
        abcd|ef). Therefore the RF distance is 2 + 2 = 4.

        >>> from skbio import TreeNode
        >>> tree1 = TreeNode.read(["((a,b),c,((d,e),f));"])
        >>> print(tree1.ascii_art())
                            /-a
                  /--------|
                 |          \-b
                 |
        ---------|--c
                 |
                 |                    /-d
                 |          /--------|
                  \--------|          \-e
                           |
                            \-f

        >>> tree2 = TreeNode.read(["((a,(b,c)),d,(e,f));"])
        >>> print(tree2.ascii_art())
                            /-a
                  /--------|
                 |         |          /-b
                 |          \--------|
                 |                    \-c
        ---------|
                 |--d
                 |
                 |          /-e
                  \--------|
                            \-f

        >>> tree1.compare_rfd(tree2)
        4.0

        """
        if rooted is None:
            rooted = self.parent is not None or len(self.children) == 2
        return _topo_dists((self, other), rooted, proportion=proportion)[0]

    def compare_wrfd(self, other, metric="cityblock", rooted=None, include_tips=True):
        r"""Calculate weighted Robinson-Foulds distance or variants between two trees.

        .. versionadded:: 0.6.3

        Parameters
        ----------
        other : TreeNode
            The other tree to compare with.
        metric : str or callable, optional
            The distance metric to use. Can be a preset, a distance function name under
            :mod:`scipy.spatial.distance`, or a custom function that takes two vectors
            and returns a number. Some notable options are:

            - "cityblock" (default): City block (Manhattan) distance. The result
              matches the original weighted Robinson-Foulds distance [1]_.
            - "euclidean": Euclidean distance. The result matches the
              Kuhner-Felsenstein (KF) distance, a.k.a. branch score (Bs) distance [2]_.
            - "correlation": 1 - Pearson's correlation coefficient (:math:`r`). Ranges
              between 0 (maximum similarity) and 2 (maximum dissimilarity). Independent
              of tree scale.
            - "unitcorr": :math:`(1 - r) / 2`, which returns a unit correlation
              distance (range: [0, 1]).

        rooted : bool, optional
            Whether to consider the trees as rooted or unrooted. If None (default),
            this will be determined based on whether self is rooted. However, one
            can override it by explicitly setting True (rooted) or False (unrooted).
            See :meth:`compare_rfd` for details.
        include_tips : bool, optional
            Whether to include single-taxon biparitions (terminal branches) in the
            calculation. Default is True, such that all branches in the trees are
            considered. Set this as False if terminal branch lengths are absent or
            irrelevant.

        Returns
        -------
        float
            The weighted Robinson-Foulds distance or variants between the trees.

        Notes
        -----
        The Robinson-Foulds (RF) distance may be weighted by the branch lengths of
        bipartitions to account for evolutionary distances in addition to branching
        patterns.

        The default behavior of this method calculates the original weighted RF (wRF)
        distance [1]_, which is the sum of differences of branch lengths of matching
        bipartitions. Bipartitions unique to one tree are given a length of 0 in the
        other tree during calculation.

        .. math::

           \text{wRF}(T_1, T_2) = \sum_{s \in S_1 \cup S_2} |l_1(s) - l_2(s)|

        where :math:`S_1` and :math:`S_2` are the sets of bipartitions of trees
        :math:`T_1` and :math:`T_2`, respectively. :math:`l_1` and :math:`l_2` are the
        branch lengths of bipartition :math:`s` in :math:`T_1` and :math:`T_2`,
        respectively (or 0 if :math:`s` is unique to the other tree).

        When ``metric="euclidean"``, it calculates the Kuhner-Felsenstein (KF)
        distance, a.k.a., branch score (Bs) distance [2]_, which replaces absolute
        difference with squared difference in the equation.

        .. math::

           \text{KF}(T_1, T_2) = \sqrt{\sum_{s \in S_1 \cup S_2} (l_1(s) - l_2(s))^2}

        This method operates on the subtrees below the given nodes. Only taxa shared
        between the two trees are considered. Taxa unique to either tree are excluded
        from the calculation.

        See Also
        --------
        compare_rfd
        compare_cophenet
        skbio.tree.wrf_dists

        References
        ----------
        .. [1] Robinson, D. F., & Foulds, L. R. (1979) Comparison of weighted labelled
           trees. In Combinatorial Mathematics VI: Proceedings of the Sixth Australian
           Conference on Combinatorial Mathematics, Armidale, Australia (pp. 119-126).

        .. [2] Kuhner, M. K., & Felsenstein, J. (1994). A simulation comparison of
           phylogeny algorithms under equal and unequal evolutionary rates. Molecular
           biology and evolution, 11(3), 459-468.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree1 = TreeNode.read(["((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);"])
        >>> print(tree1.ascii_art())
                            /-a
                  /--------|
                 |          \-b
                 |
        ---------|--c
                 |
                 |                    /-d
                 |          /--------|
                  \--------|          \-e
                           |
                            \-f

        >>> tree2 = TreeNode.read(["((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);"])
        >>> print(tree2.ascii_art())
                            /-a
                  /--------|
                 |         |          /-b
                 |          \--------|
                 |                    \-c
        ---------|
                 |--d
                 |
                 |          /-e
                  \--------|
                            \-f

        Calculate the weighted RF (wRF) distance between two unrooted trees with branch
        lengths.

        >>> tree1.compare_wrfd(tree2)
        16.0

        Calculated the wRF distance while considering trees as rooted (therefore based
        on subsets instead of bipartitions).

        >>> tree1.compare_wrfd(tree2, rooted=True)
        18.0

        Calculate the Kuhner-Felsenstein (KF) distance.

        >>> d = tree1.compare_wrfd(tree2, metric="euclidean")
        >>> print(round(d, 5))
        6.16441

        Calculate the KF distance without considering terminal branches.

        >>> d = tree1.compare_wrfd(tree2, metric="euclidean", include_tips=False)
        >>> print(round(d, 5))
        3.74166

        """
        if rooted is None:
            rooted = self.parent is not None or len(self.children) == 2
        metric = _check_dist_metric(metric)
        return _topo_dists(
            (self, other),
            rooted=rooted,
            include_tips=include_tips,
            weighted=True,
            metric=metric,
        )[0]

    @aliased("compare_tip_distances", "0.6.3")
    @params_aliased(
        [
            ("shuffler", "shuffle_f", "0.6.3", True),
            ("metric", "dist_f", "0.6.3", True),
        ]
    )
    def compare_cophenet(
        self,
        other,
        sample=None,
        metric="unitcorr",
        shuffler=None,
        use_length=True,
        ignore_self=False,
    ):
        r"""Calculate the distance between two trees based on cophenetic distances.

        Parameters
        ----------
        other : TreeNode
            The other tree to compare with.
        sample : int, optional
            Randomly subsample this number of tips in common between the trees to
            compare. This is useful when comparing very large trees.
        metric : str or callable, optional
            The distance metric to use. Can be a preset, a distance function name under
            :mod:`scipy.spatial.distance`, or a custom function that takes two vectors
            and returns a number. Some notable options are:

            - "cityblock": City block (Manhattan) distance.
            - "euclidean": Euclidean distance. The result matches the path-length
              distance [1]_, or the path distance [2]_ if ``use_length`` is False.
            - "correlation": 1 - Pearson's correlation coefficient (:math:`r`). Ranges
              between 0 (maximum similarity) and 2 (maximum dissimilarity). Independent
              of tree scale.
            - "unitcorr" (default): :math:`(1 - r) / 2`, which returns a unit
              correlation distance (range: [0, 1]).

            .. versionchanged:: 0.6.3
                Accepts a function on two vectors instead of two `DistanceMatrix`
                instances. The default value "unitcorr" is consistent with the previous
                default behavior.

        shuffler : int, np.random.Generator or callable, optional
            The shuffling function to use if ``sample`` is specified. Default is
            :meth:`~numpy.random.Generator.shuffle`. If an integer is provided, a
            random generator will be constructed using this number as the seed.

            .. versionchanged:: 0.6.3
                Switched to NumPy's new random generator. Can accept a random seed or
                random generator instance.

        use_length : bool, optional
            Whether to calculate the sum of branch lengths (True, default) or the
            number of branches (False) connecting each pair of tips.

            .. versionadded:: 0.6.3

        ignore_self : bool, optional
            Whether to ignore the distance between each tip and itself (which must be
            0). Default is False.

            .. versionadded:: 0.6.3

            .. note:: The default value will be set as True in 0.7.0.

        Returns
        -------
        float
            The distance between the trees.

        .. versionchanged:: 0.6.3
            Improved customizability to allow calculation of published metrics, such
            as path distance and path-length distance, while preserving the previous
            default behavior.

            Edge cases are now handled by the specified distance metric rather than
            being treated separately.

        Raises
        ------
        ValueError
            If there are no common tips between the trees.

        See Also
        --------
        cophenet
        compare_rfd
        compare_wrfd
        skbio.tree.path_dists

        Notes
        -----
        This method calculates the dissimilarity between the cophenetic distance [1]_
        (i.e., tip-to-tip distance) matrices of two trees. Tips are identified by
        their names (i.e., taxa). Only tips shared between the trees are considered.
        Tips unique to either tree are excluded from the calculation.

        The default behavior returns a unit correlation distance (range: [0, 1]),
        measuring the dissimilarity between the relative evolutionary distances among
        taxa, regardless of the tree scale (i.e., multiply all branch lengths in one
        tree by a factor and the result remains the same). This measure is closely
        related to **cophenetic correlation**, which measures the similarity (instead
        of dissimilarity) between two cophenetic distance matrices, or between a
        cophenetic distance matrix and the original distance matrix among taxa on
        which hierarchical clustering was performed.

        When the metric is Euclidean and lengths are used, it returns the **path-length
        distance** [2]_, which is the square root of the sum of squared differences of
        path lengths among all pairs of taxa.

        .. math::

           d(T_1, T_2) = \sqrt{\sum (d_1(i,j) - d_2(i,j))^2}

        where :math:`d_1` and :math:`d_2` are the sums of branch lengths connecting a
        pair of tips :math:`i` and :math:`j` in trees :math:`T_1` and :math:`T_2`,
        respectively.

        When the metric is Euclidean and lengths are not used, it returns the **path
        distance** [3]_, which insteads considers the number of edges in the path.

        References
        ----------
        .. [1] Sokal, R. R., & Rohlf, F. J. (1962). The comparison of dendrograms by
           objective methods. Taxon, 33-40.

        .. [2] Lapointe, F. J., & Cucumel, G. (1997). The average consensus procedure:
           combination of weighted trees containing identical or overlapping sets of
           taxa. Systematic Biology, 46(2), 306-312.

        .. [3] Steel, M. A., & Penny, D. (1993). Distributions of tree comparison
           metrics—some new results. Systematic Biology, 42(2), 126-141.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree1 = TreeNode.read(["((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);"])
        >>> print(tree1.ascii_art())
                            /-a
                  /--------|
                 |          \-b
                 |
        ---------|--c
                 |
                 |                    /-d
                 |          /--------|
                  \--------|          \-e
                           |
                            \-f

        >>> tree2 = TreeNode.read(["((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);"])
        >>> print(tree2.ascii_art())
                            /-a
                  /--------|
                 |         |          /-b
                 |          \--------|
                 |                    \-c
        ---------|
                 |--d
                 |
                 |          /-e
                  \--------|
                            \-f

        Calculate the unit correlation distance between the two trees.

        >>> d = tree1.compare_cophenet(tree2, ignore_self=True)
        >>> print(round(d, 5))
        0.14131

        Calculate the path-length distance between the two trees.

        >>> d = tree1.compare_cophenet(tree2, metric="euclidean",
        ...                                 ignore_self=True)
        >>> print(round(d, 5))
        13.71131

        Calculate the path distance between the two trees.

        >>> tree1.compare_cophenet(
        ...     tree2, metric="euclidean", use_length=False, ignore_self=True)
        4.0

        """
        # future warning
        if ignore_self is False:
            _warn_once(
                self.__class__.compare_cophenet,
                FutureWarning,
                "The default behavior of `compare_cophenet` is subject to change in "
                "0.7.0. The new default behavior can be achieved by specifying "
                "`ignore_self=True`.",
            )

        metric = _check_dist_metric(metric)
        if sample is not None:
            shuffler = _check_shuffler(shuffler)

        return _path_dists(
            trees=(self, other),
            sample=sample,
            metric=metric,
            shuffler=shuffler,
            use_length=use_length,
            ignore_self=ignore_self,
        )[0]

    # ------------------------------------------------
    # Tree indexing and searching
    # ------------------------------------------------

    def has_caches(self):
        r"""Check if the current tree has caches.

        .. versionadded:: 0.6.3

        Returns
        -------
        set or None
            Names of present node attribute caches, or None if none is present.
        bool
            Presence (True) or absence (False) of node lookup caches.

        See Also
        --------
        clear_caches
        cache_attr
        find

        Notes
        -----
        Caches are optional but can significantly accelerate certain analyses of the
        tree. Two types of caches may be created:

        1. **Node attributes**, which may be created by calling :meth:`cache_attr` and
           assigned to individual nodes within the tree. The names of these attributes
           are optionally registered at the root.

        2. **Node lookup table**, which is automatically created during the first node
           search (e.g., by calling :meth:`find`) and reused in subsequent searches.
           This table is attached to the root of the tree.

        This method checks if a node lookup table and any registered node attributes
        are present at the root node of the tree.

        The returned set of node attribute names is a reference instead of a copy. One
        may edit the set in place to explicitly enable/disable names. Use this feature
        with caution.

        When the tree is manipulated, caches typically become obsolete and are
        automatically cleared. If the caches are not present or not relevant to the
        analysis, you may set ``uncache=False`` when performing individual operations
        to suppress clearing. This can improve the performance of these operations.

        You may explicitly call :meth:`clear_caches` to clear caches of a tree.

        """
        tree = self.root()
        attrs = getattr(tree, "_registered_caches", None)
        lookup = hasattr(tree, "_tip_cache") and hasattr(tree, "_non_tip_cache")
        return attrs, lookup

    @aliased("invalidate_caches", "0.6.3", True)
    def clear_caches(self, attr=True, lookup=True):
        r"""Delete node attribute and lookup caches of a tree.

        Parameters
        ----------
        attr : bool or str, optional
            Whether to delete attribute caches created by ``cache_attr`` (default:
            True). One may instead provide an attribute name such that only this
            attribute will be deleted.

            .. versionchanged:: 0.6.3

                Can provide a specific attribute name.

        lookup : bool, optional
            Whether to delete lookup caches created during name searching (default:
            True).

            .. versionadded:: 0.6.3

        See Also
        --------
        has_caches
        cache_attr
        find

        Notes
        -----
        This method may be called from any node within a tree. The caches, which were
        attached to the root node of the tree, will be deleted.

        This method will silently skip if the specified caches do not exist.

        """
        tree = self.root()

        # delete attribute caches
        if attr and hasattr(tree, "_registered_caches"):
            attrs = tree._registered_caches

            # delete a single attribute
            if isinstance(attr, str):
                if attr not in attrs:
                    return
                for node in tree.traverse():
                    if hasattr(node, attr):
                        delattr(node, attr)
                if len(attrs) == 1:
                    del tree._registered_caches
                else:
                    attrs.remove(attr)

            # delete all attributes
            else:
                for node in tree.traverse():
                    for attr in attrs:
                        if hasattr(node, attr):
                            delattr(node, attr)
                del tree._registered_caches

        # delete lookup caches
        if lookup:
            for key in ("_tip_cache", "_non_tip_cache"):
                if hasattr(tree, key):
                    delattr(tree, key)

    def cache_attr(self, func, cache_attrname, cache_type=list, register=True):
        r"""Cache attributes on nodes of the tree through a postorder traversal.

        Parameters
        ----------
        func : callable
            Function to calculate the attribute of the current node. The result will be
            combined with the attributes of the previous nodes, if applicable.

        cache_attrname : str
            Name of the attribute to be attached to each node.

        cache_type : {list, tuple, set, frozenset}, callable, or None
            The type of the cache. Can be any of the four iterable types: list
            (default), tuple, set, or frozenset. In these cases, combination of
            attributes of the node's children and itself will be automated.

            Or a custom function that takes two arguments: list of attributes of its
            children, and attribute calculated from itself by ``func``, and returns the
            combined attribute of the node.

            Or None, in which case combination of attributes of children and self
            will not take place, unless explicitly customized within ``func``.

            .. versionchanged:: 0.6.3

                Tuple, custom function and None were added to the options.

        register : bool, optional
            Whether to register the attribute name as a cache of the tree, such that
            the attributes will be deleted from all nodes when the tree is manipulated
            or the ``clear_caches`` method is explicitly invoked. Default is True.

            .. versionadded:: 0.6.3

        Raises
        ------
        TypeError
            If ``cache_type`` is invalid.

        See Also
        --------
        has_caches
        clear_caches

        Notes
        -----
        This method provides an efficient interface to assign a custom attribute to
        every node of a tree through one postorder travesal. It is particularly useful
        if one needs to frequently look up attributes that would normally require one
        traversal of the tree per lookup. The assigned attributes may be automatically
        deleted when the tree is manipulated.

        Examples
        --------
        This method facilitates evaluation for various useful node properties. Some
        representative examples are provided below.

        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a:1.2,b:1.6)c:0.3,(d:0.8,e:1.0)f:0.6)g;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -g-------|
                 |          /-d
                  \f-------|
                            \-e

        Cache a list of all descending tip names on each node. This faciliates the
        assignment of taxon set under each clade in the tree. It resembles but is more
        efficient than calling :meth:`subset` multiple times.

        >>> f = lambda n: [n.name] if n.is_tip() else []
        >>> tree.cache_attr(f, 'tip_names')
        >>> for node in tree.traverse(include_self=True):
        ...     print(f"Node: {node.name}, tips: {node.tip_names}")
        Node: g, tips: ['a', 'b', 'd', 'e']
        Node: c, tips: ['a', 'b']
        Node: a, tips: ['a']
        Node: b, tips: ['b']
        Node: f, tips: ['d', 'e']
        Node: d, tips: ['d']
        Node: e, tips: ['e']

        Cache the number of nodes per clade. The function ``sum`` is used in place of
        cache type such that the count will be accumulated. This resembles but is more
        efficient than calling :meth:`count` multiple times.

        >>> f = lambda n: 1
        >>> tree.cache_attr(f, 'node_count', sum)
        >>> tree.node_count
        7

        Cache the sum of branch lengths per clade. This resembles but is more efficient
        than calling :meth:`total_length` multiple times.

        >>> f = lambda n: n.length or 0.0
        >>> tree.cache_attr(f, 'clade_size', sum)
        >>> tree.clade_size
        5.5

        Cache the accumulative distances from all tips to the common ancestor of each
        clade. This is more efficient than calling :meth:`depth` multiple times. One
        can further apply calculations like mean and standard deviation to the results.

        >>> import numpy as np
        >>> dist_f = lambda n: np.array(n.length or 0.0, ndmin=1)
        >>> comb_f = lambda prev, curr: np.concatenate(prev) + curr if prev else curr
        >>> tree.cache_attr(dist_f, 'accu_dists', comb_f)
        >>> tree.accu_dists
        array([ 1.5,  1.9,  1.4,  1.6])

        """
        if cache_type in (set, frozenset):

            def combine_f(prev, curr):
                return cache_type().union(*prev + [curr])

        elif cache_type in (list, tuple):

            def combine_f(prev, curr):
                return cache_type(chain.from_iterable(prev + [curr]))

        elif callable(cache_type) or cache_type is None:
            combine_f = cache_type
        else:
            raise TypeError("Cache type is invalid.")

        # register attribute name as a cache
        if register:
            tree = self.root()
            if not hasattr(tree, "_registered_caches"):
                tree._registered_caches = set()
            tree._registered_caches.add(cache_attrname)

        # traverse tree and assign attributes
        if combine_f is not None:
            for node in self.postorder(include_self=True):
                prev_attrs = [getattr(c, cache_attrname) for c in node.children]
                curr_attr = func(node)
                setattr(node, cache_attrname, combine_f(prev_attrs, curr_attr))
        else:
            for node in self.postorder(include_self=True):
                setattr(node, cache_attrname, func(node))

    def assign_ids(self):
        """Assign topologically stable unique IDs to all nodes of the tree.

        See Also
        --------
        find_by_id
        postorder

        Notes
        -----
        This method assigns unique IDs to all nodes of the tree via a postorder
        traversal. The IDs are incremental integers starting from 0. The order is
        topologically stable. Following the call, all nodes in the tree will have
        their ``id`` attribute set.

        """
        curr_index = 0
        for n in self.postorder():
            for c in n.children:
                c.id = curr_index
                curr_index += 1

        self.id = curr_index

    def index_tree(self):
        r"""Index a tree for rapid lookups within a tree array.

        Indexes nodes in-place as ``n._leaf_index``.

        Returns
        -------
        dict
            A mapping {node_id: TreeNode}
        ndarray of int
            This arrays describes the IDs of every internal node, and the ID
            range of the immediate descendents. The first column in the array
            corresponds to node_id. The second column is the left most
            descendent's ID. The third column is the right most descendent's
            ID.

        """
        self.assign_ids()

        id_index = {}
        child_index = []

        for n in self.postorder():
            for c in n.children:
                id_index[c.id] = c

                if c:
                    # c has children itself, so need to add to result
                    child_index.append((c.id, c.children[0].id, c.children[-1].id))

        # handle root, which should be t itself
        id_index[self.id] = self

        # only want to add to the child_index if self has children...
        if self.children:
            child_index.append((self.id, self.children[0].id, self.children[-1].id))
        child_index = np.asarray(child_index, dtype=np.int64)
        child_index = np.atleast_2d(child_index)

        return id_index, child_index

    def create_caches(self):
        r"""Construct an internal lookup table to facilitate searching by name.

        Raises
        ------
        DuplicateNodeError
            If there are duplicate tip names.

        See Also
        --------
        has_caches
        clear_caches
        find

        Notes
        -----
        This method is automatically called during the first search in a tree (methods
        :meth:`find` and :meth:`find_all`). After that, subsequent searches will
        utilize the constructed lookup table, until it is deleted.

        This method may be called from any node within a tree. The lookup table will be
        attached to the root node of the tree.

        This method will not cache nodes whose name is ``None``. This method will
        raise an error if a name conflict in the tips is discovered, but will not raise
        if on internal nodes. This is because, in practice, the tips of a tree are
        required to be unique while no such requirement holds for internal nodes.

        """
        tree = self.root()
        if hasattr(tree, "_tip_cache") and hasattr(tree, "_non_tip_cache"):
            return

        tip_cache, non_tip_cache = {}, {}
        non_tip_cache_setdefault = non_tip_cache.setdefault
        for node in tree.postorder():
            if (name := node.name) is None:
                continue
            if node.is_tip():
                if name in tip_cache:
                    raise DuplicateNodeError(f"Duplicate tip name '{name}' found.")
                tip_cache[name] = node
            else:
                non_tip_cache_setdefault(name, []).append(node)

        tree._tip_cache = tip_cache
        tree._non_tip_cache = non_tip_cache

    def find(self, name):
        r"""Find a node by name.

        Parameters
        ----------
        name : TreeNode or str
            The name of the node to find. If a ``TreeNode`` object is provided,
            will find this particular node in the tree.

        Returns
        -------
        TreeNode
            The found node.

        Raises
        ------
        MissingNodeError
            If the node to be searched for is not found in the current tree.

        See Also
        --------
        find_all
        find_by_id
        find_by_func

        Notes
        -----
        This method will first attempt to find the node in the tips. If it cannot find
        a corresponding tip, it will then search through the internal nodes of the
        tree. In practice, phylogenetic trees and other common trees in biology do not
        have unique internal node names. As a result, this find method will only return
        the first occurrence of an internal node encountered on a postorder traversal
        of the tree.

        The first call of ``find`` will cache a node lookup table in the tree on the
        assumption that additional calls to ``find`` will be made. See
        :meth:`details <has_caches>`.

        This method searches within the entire tree where self is located, regardless
        if self is the root node.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f);"])
        >>> node = tree.find('c')
        >>> node.name
        'c'

        """
        tree = self.root()

        # create lookup table if not already
        tree.create_caches()

        # if input is a node, get its name
        name_is_node = isinstance(name, tree.__class__)
        name_ = name.name if name_is_node else name

        # look up name in tips
        node = tree._tip_cache.get(name_, None)
        if node is not None:
            if not name_is_node or node is name:
                return node

        # look up name in non-tips
        nodes = tree._non_tip_cache.get(name_, None)
        if nodes is not None:
            if name_is_node:
                for node in nodes:
                    if node is name:
                        return node
            else:
                return nodes[0]

        raise MissingNodeError(f"Node '{name_}' is not found in the tree.")

    def find_all(self, name):
        r"""Find all nodes that match a given name.

        Parameters
        ----------
        name : TreeNode or str
            The name or node to find. If a ``TreeNode`` object is provided, all nodes
            with the same name will be returned.

        Returns
        -------
        list of TreeNode
            The found nodes.

        Raises
        ------
        MissingNodeError
            If the node to be searched for is not found.

        See Also
        --------
        find
        find_by_id
        find_by_func

        Notes
        -----
        All internal nodes (including root) and tips with the given name will be
        returned, with the former placed before the latter in the returned list.

        The first call to ``find_all`` will cache a node lookup table in the tree on
        the assumption that additional calls to ``find_all`` will be made. See
        :meth:`details <has_caches>`.

        This method searches within the entire tree where self is located, regardless
        if self is the root node.

        Examples
        --------
        >>> from skbio.tree import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)d,(f,g)c);"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
                 |
                 |          /-d
        ---------|-d-------|
                 |          \-e
                 |
                 |          /-f
                  \c-------|
                            \-g

        >>> for node in tree.find_all('c'):
        ...     print(node.name, node.children[0].name, node.children[1].name)
        c a b
        c f g
        >>> for node in tree.find_all('d'):
        ...     print(node.name, str(node))
        d (d,e)d;
        <BLANKLINE>
        d d;
        <BLANKLINE>

        """
        tree = self.root()
        if isinstance(name, tree.__class__):
            name = name.name
        tree.create_caches()
        tip = tree._tip_cache.get(name, None)
        nodes = tree._non_tip_cache.get(name, [])
        if tip is not None:
            nodes.append(tip)
        if not nodes:
            raise MissingNodeError(f"Node '{name}' is not found.")
        else:
            return nodes

    def find_by_id(self, node_id):
        r"""Find a node by ID.

        Parameters
        ----------
        node_id : int
            The ID of a node in the tree.

        Returns
        -------
        TreeNode
            The node with the matching ID.

        Raises
        ------
        MissingNodeError
            If the ID cannot be found.

        See Also
        --------
        assign_ids
        find

        Notes
        -----
        This method searches within the subtree under the current node. But the IDs
        are assigned from the root of the entire tree.

        This method does not cache ID associations. A full traversal of the
        tree is performed to find a node by an ID on every call.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f);"])
        >>> print(tree.find_by_id(2).name)
        d

        """
        self.root().assign_ids()
        for node in self.traverse(include_self=True):
            if node.id == node_id:
                return node
        raise MissingNodeError(f"ID {node_id} is not in self.")

    def find_by_func(self, func):
        r"""Find all nodes in a tree that meet certain criteria.

        Parameters
        ----------
        func : callable
            A function that accepts a ``TreeNode`` and returns True or False, where
            True indicates the node is to be yielded.

        Yields
        ------
        TreeNode
            The found node.

        See Also
        --------
        find
        find_all
        find_by_id

        Notes
        -----
        This search method is based on the current subtree, not the root.

        This method does not cache search results.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f);"])
        >>> func = lambda x: x.parent == tree.find('c')
        >>> [n.name for n in tree.find_by_func(func)]
        ['a', 'b']

        """
        for node in self.traverse(include_self=True):
            if func(node):
                yield node

    # ------------------------------------------------
    # Tree visualization
    # ------------------------------------------------

    def _ascii_art(self, char1="-", show_internal=True, compact=False):
        LEN = 10
        PAD = " " * LEN
        PA = " " * (LEN - 1)
        namestr = self._node_label()
        if self.children:
            mids = []
            result = []
            for c in self.children:
                if c is self.children[0]:
                    char2 = "/"
                elif c is self.children[-1]:
                    char2 = "\\"
                else:
                    char2 = "-"
                (clines, mid) = c._ascii_art(char2, show_internal, compact)
                mids.append(mid + len(result))
                result.extend(clines)
                if not compact:
                    result.append("")
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = (
                [PAD] * (lo + 1) + [PA + "|"] * (hi - lo - 1) + [PAD] * (end - hi)
            )
            mid = int(np.trunc((lo + hi) / 2))
            prefixes[mid] = char1 + "-" * (LEN - 2) + prefixes[mid][-1]
            result = [p + L for (p, L) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + namestr + stem[len(namestr) + 1 :]
            return (result, mid)
        else:
            return ([char1 + "-" + namestr], 0)

    def ascii_art(self, show_internal=True, compact=False):
        r"""Return a string containing an ascii drawing of the tree.

        Note, this method calls a private recursive function and is not safe
        for large trees.

        Parameters
        ----------
        show_internal : bool
            includes internal edge names
        compact : bool
            use exactly one line per tip

        Returns
        -------
        str
            an ASCII formatted version of the tree

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,b)c,(d,e)f)root;"])
        >>> print(tree.ascii_art())
                            /-a
                  /c-------|
                 |          \-b
        -root----|
                 |          /-d
                  \f-------|
                            \-e

        """
        (lines, mid) = self._ascii_art(show_internal=show_internal, compact=compact)
        return "\n".join(lines)

    # ------------------------------------------------
    # Format conversion
    # ------------------------------------------------

    def _balanced_distance_to_tip(self):
        """Return the distance to tip from this node.

        The distance to every tip from this node must be equal for this to
        return a correct result.

        Returns
        -------
        float
            The distance to tip of a length-balanced tree.

        """
        node = self
        distance = 0.0
        while node.has_children():
            distance += node.children[0].length
            node = node.children[0]
        return distance

    @classonlymethod
    def from_linkage_matrix(cls, linkage_matrix, id_list):
        r"""Return tree from SciPy linkage matrix.

        Parameters
        ----------
        linkage_matrix : ndarray
            A linkage matrix generated by ``scipy.cluster.hierarchy.linkage``.
        id_list : list
            Corresponding IDs of the indices in the linkage matrix.

        Returns
        -------
        TreeNode
            An unrooted bifurcated tree.

        See Also
        --------
        scipy.cluster.hierarchy.linkage

        """
        tip_width = len(id_list)
        cluster_count = len(linkage_matrix)
        lookup_len = cluster_count + tip_width
        node_lookup = np.empty(lookup_len, dtype=cls)

        for i, name in enumerate(id_list):
            node_lookup[i] = cls(name=name)

        for i in range(tip_width, lookup_len):
            node_lookup[i] = cls()

        newest_cluster_index = cluster_count + 1
        for link in linkage_matrix:
            child_a = node_lookup[int(link[0])]
            child_b = node_lookup[int(link[1])]

            path_length = link[2] / 2
            child_a.length = path_length - child_a._balanced_distance_to_tip()
            child_b.length = path_length - child_b._balanced_distance_to_tip()

            new_cluster = node_lookup[newest_cluster_index]
            new_cluster.append(child_a, uncache=False)
            new_cluster.append(child_b, uncache=False)

            newest_cluster_index += 1

        return node_lookup[-1]

    @classonlymethod
    def from_taxonomy(cls, lineage_map):
        r"""Construct a tree from a taxonomy.

        Parameters
        ----------
        lineage_map : dict, iterable of tuples, or pd.DataFrame
            Mapping of taxon IDs to lineages (iterables of taxonomic units
            from high to low in ranking).

        Returns
        -------
        TreeNode
            The constructed taxonomy.

        See Also
        --------
        from_taxdump

        Examples
        --------
        >>> from skbio.tree import TreeNode
        >>> lineages = [
        ...     ('1', ['Bacteria', 'Firmicutes', 'Clostridia']),
        ...     ('2', ['Bacteria', 'Firmicutes', 'Bacilli']),
        ...     ('3', ['Bacteria', 'Bacteroidetes', 'Sphingobacteria']),
        ...     ('4', ['Archaea', 'Euryarchaeota', 'Thermoplasmata']),
        ...     ('5', ['Archaea', 'Euryarchaeota', 'Thermoplasmata']),
        ...     ('6', ['Archaea', 'Euryarchaeota', 'Halobacteria']),
        ...     ('7', ['Archaea', 'Euryarchaeota', 'Halobacteria']),
        ...     ('8', ['Bacteria', 'Bacteroidetes', 'Sphingobacteria']),
        ...     ('9', ['Bacteria', 'Bacteroidetes', 'Cytophagia'])]

        >>> tree = TreeNode.from_taxonomy(lineages)
        >>> print(tree.ascii_art())
                                      /Clostridia-1
                            /Firmicutes
                           |          \Bacilli- /-2
                  /Bacteria|
                 |         |                    /-3
                 |         |          /Sphingobacteria
                 |          \Bacteroidetes      \-8
                 |                   |
        ---------|                    \Cytophagia-9
                 |
                 |                              /-4
                 |                    /Thermoplasmata
                 |                   |          \-5
                  \Archaea- /Euryarchaeota
                                     |          /-6
                                      \Halobacteria
                                                \-7

        """
        root = cls(name=None)
        root._lookup = {}

        if isinstance(lineage_map, dict):
            lineage_map = lineage_map.items()
        elif isinstance(lineage_map, pd.DataFrame):
            lineage_map = ((idx, row.tolist()) for idx, row in lineage_map.iterrows())

        for id_, lineage in lineage_map:
            cur_node = root

            # for each name, see if we've seen it, if not, add that puppy on
            for name in lineage:
                if name in cur_node._lookup:
                    cur_node = cur_node._lookup[name]
                else:
                    new_node = cls(name=name)
                    new_node._lookup = {}
                    cur_node._lookup[name] = new_node
                    cur_node.append(new_node, uncache=False)
                    cur_node = new_node

            cur_node.append(cls(name=id_), uncache=False)

        # scrub the lookups
        for node in root.non_tips(include_self=True):
            del node._lookup

        return root

    def to_taxonomy(self, allow_empty=False, filter_f=None):
        """Return a taxonomy representation of self.

        Parameters
        ----------
        allow_empty : bool, optional
            Allow gaps the taxonomy (e.g., internal nodes without names).
        filter_f : function, optional
            Specify a filtering function that returns True if the lineage is
            to be returned. This function must accept a ``TreeNode`` as its
            first parameter, and a ``list`` that represents the lineage as the
            second parameter.

        Yields
        ------
        tuple
            ``(tip, [lineage])`` where ``tip`` corresponds to a tip in the tree
            and ``[lineage]`` is the expanded names from root to tip. ``None``
            and empty strings are omitted from the lineage.

        Notes
        -----
        If ``allow_empty`` is True and the root node does not have a name, that name
        will not be included. This is because it is common to have multiple domains
        represented in the taxonomy, which would result in a root node that does not
        have a name and does not make sense to represent in the output.

        Examples
        --------
        >>> from skbio.tree import TreeNode
        >>> lineages = {'1': ['Bacteria', 'Firmicutes', 'Clostridia'],
        ...             '2': ['Bacteria', 'Firmicutes', 'Bacilli'],
        ...             '3': ['Bacteria', 'Bacteroidetes', 'Sphingobacteria'],
        ...             '4': ['Archaea', 'Euryarchaeota', 'Thermoplasmata'],
        ...             '5': ['Archaea', 'Euryarchaeota', 'Thermoplasmata'],
        ...             '6': ['Archaea', 'Euryarchaeota', 'Halobacteria'],
        ...             '7': ['Archaea', 'Euryarchaeota', 'Halobacteria'],
        ...             '8': ['Bacteria', 'Bacteroidetes', 'Sphingobacteria'],
        ...             '9': ['Bacteria', 'Bacteroidetes', 'Cytophagia']}
        >>> tree = TreeNode.from_taxonomy(lineages.items())
        >>> lineages = sorted([(n.name, l) for n, l in tree.to_taxonomy()])
        >>> for name, lineage in lineages:
        ...     print(name, '; '.join(lineage))
        1 Bacteria; Firmicutes; Clostridia
        2 Bacteria; Firmicutes; Bacilli
        3 Bacteria; Bacteroidetes; Sphingobacteria
        4 Archaea; Euryarchaeota; Thermoplasmata
        5 Archaea; Euryarchaeota; Thermoplasmata
        6 Archaea; Euryarchaeota; Halobacteria
        7 Archaea; Euryarchaeota; Halobacteria
        8 Bacteria; Bacteroidetes; Sphingobacteria
        9 Bacteria; Bacteroidetes; Cytophagia

        """
        if filter_f is None:

            def filter_f(a, b):
                return True

        self.assign_ids()
        seen = set()
        seen_add = seen.add
        lineage = []
        lineage_pop = lineage.pop
        lineage_append = lineage.append

        # visit internal nodes while traversing out to the tips, and on the
        # way back up
        for node in self.traverse(self_before=True, self_after=True):
            if node.is_tip():
                if filter_f(node, lineage):
                    yield (node, lineage[:])
            else:
                if allow_empty:
                    if node.is_root() and not node.name:
                        continue
                else:
                    if not node.name:
                        continue

                if node.id in seen:
                    lineage_pop()
                else:
                    lineage_append(node.name)
                    seen_add(node.id)

    @classonlymethod
    def from_taxdump(cls, nodes, names=None):
        r"""Construct a tree from the NCBI taxonomy database.

        Parameters
        ----------
        nodes : pd.DataFrame
            Taxon hierarchy.
        names : pd.DataFrame or dict, optional
            Taxon names.

        Returns
        -------
        TreeNode
            The constructed tree.

        Notes
        -----
        ``nodes`` and ``names`` correspond to "nodes.dmp" and "names.dmp" of
        the NCBI taxonomy database. The should be read into data frames using
        ``skbio.io.read`` prior to this operation. Alternatively, ``names``
        may be provided as a dictionary. If ``names`` is omitted, taxonomy IDs
        be used as taxon names.

        Raises
        ------
        ValueError
            If there is no top-level node.
        ValueError
            If there are more than one top-level node.

        See Also
        --------
        from_taxonomy
        skbio.io.format.taxdump

        Examples
        --------
        >>> import pandas as pd
        >>> from skbio.tree import TreeNode
        >>> nodes = pd.DataFrame([
        ...             [1, 1, 'no rank'],
        ...             [2, 1, 'domain'],
        ...             [3, 1, 'domain'],
        ...             [4, 2, 'phylum'],
        ...             [5, 2, 'phylum']], columns=[
        ...     'tax_id', 'parent_tax_id', 'rank']).set_index('tax_id')
        >>> names = {1: 'root', 2: 'Bacteria', 3: 'Archaea',
        ...          4: 'Firmicutes', 5: 'Bacteroidetes'}
        >>> tree = TreeNode.from_taxdump(nodes, names)
        >>> print(tree.ascii_art())
                            /-Firmicutes
                  /Bacteria|
        -root----|          \-Bacteroidetes
                 |
                  \-Archaea

        """
        # identify top level of hierarchy
        tops = nodes[nodes["parent_tax_id"] == nodes.index]

        # validate root uniqueness
        n_top = tops.shape[0]
        if n_top == 0:
            raise ValueError("There is no top-level node.")
        elif n_top > 1:
            raise ValueError("There are more than one top-level node.")

        # get root taxid
        root_id = tops.index[0]

        # get parent-to-child(ren) map
        to_children = {
            p: g.index.tolist()
            for p, g in nodes[nodes.index != root_id].groupby("parent_tax_id")
        }

        # get rank map
        ranks = nodes["rank"].to_dict()

        # get taxon-to-name map
        # if not provided, use tax_id as name
        if names is None:
            names = {x: str(x) for x in nodes.index}

        # use "scientific name" as name
        elif isinstance(names, pd.DataFrame):
            names = names[names["name_class"] == "scientific name"][
                "name_txt"
            ].to_dict()

        # initiate tree
        tree = cls(names[root_id])
        tree.id = root_id
        tree.rank = ranks[root_id]

        # helper for extending tree
        def _extend_tree(node):
            self_id = node.id
            if self_id not in to_children:
                return
            children = []
            for id_ in to_children[self_id]:
                child = TreeNode(names[id_])
                child.id = id_
                child.rank = ranks[id_]
                _extend_tree(child)
                children.append(child)
            node.extend(children, uncache=False)

        # extend tree
        _extend_tree(tree)
        return tree

    def to_array(self, attrs=None, nan_length_value=None):
        """Return an array representation of self.

        Parameters
        ----------
        attrs : list of tuple or None
            The attributes and types to return. The expected form is
            [(attribute_name, type)]. If `None`, then `name`, `length`, and
            `id` are returned.
        nan_length_value : float, optional
            If provided, replaces any `nan` in the branch length vector
            (i.e., ``result['length']``) with this value. `nan` branch lengths
            can arise from an edge not having a length (common for the root
            node parent edge), which can making summing problematic.

        Returns
        -------
        dict of array
            {id_index: {id: TreeNode},
             child_index: ((node_id, left_child_id, right_child_id)),
             attr_1: array(...),
             ...
             attr_N: array(...)}

        Notes
        -----
        Attribute arrays are in index order such that TreeNode.id can be used
        as a lookup into the array.

        Examples
        --------
        >>> from skbio import TreeNode
        >>> t = TreeNode.read(['(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7);'])
        >>> res = t.to_array()
        >>> sorted(res.keys())
        ['child_index', 'id', 'id_index', 'length', 'name']
        >>> res['child_index'] # doctest: +ELLIPSIS
        array([[4, 0, 2],
               [5, 3, 3],
               [6, 4, 5],
               [7, 6, 6]]...
        >>> for k, v in res['id_index'].items():
        ...     print(k, v)
        ...
        0 a:1.0;
        <BLANKLINE>
        1 b:2.0;
        <BLANKLINE>
        2 c:3.0;
        <BLANKLINE>
        3 d:5.0;
        <BLANKLINE>
        4 (a:1.0,b:2.0,c:3.0)x:4.0;
        <BLANKLINE>
        5 (d:5.0)y:6.0;
        <BLANKLINE>
        6 ((a:1.0,b:2.0,c:3.0)x:4.0,(d:5.0)y:6.0)z:7.0;
        <BLANKLINE>
        7 (((a:1.0,b:2.0,c:3.0)x:4.0,(d:5.0)y:6.0)z:7.0);
        <BLANKLINE>
        >>> res['id']
        array([0, 1, 2, 3, 4, 5, 6, 7])
        >>> res['name']
        array(['a', 'b', 'c', 'd', 'x', 'y', 'z', None], dtype=object)

        """
        if attrs is None:
            attrs = [("name", object), ("length", float), ("id", int)]
        else:
            for attr, dtype in attrs:
                if not hasattr(self, attr):
                    raise AttributeError("Invalid attribute '%s'." % attr)

        id_index, child_index = self.index_tree()
        n = self.id + 1  # assign_ids starts at 0
        tmp = [np.zeros(n, dtype=dtype) for attr, dtype in attrs]

        for node in self.traverse(include_self=True):
            n_id = node.id
            for idx, (attr, dtype) in enumerate(attrs):
                tmp[idx][n_id] = getattr(node, attr)

        results = {"id_index": id_index, "child_index": child_index}
        results.update({attr: arr for (attr, dtype), arr in zip(attrs, tmp)})
        if nan_length_value is not None:
            length_v = results["length"]
            length_v[np.isnan(length_v)] = nan_length_value
        return results
