# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Derived from improved-octo-waddle (https://github.com/biocore/improved-octo-waddle)
# originally authored by Daniel McDonald, distributed under the Modified BSD License.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

### NOTE: some doctext strings are copied and pasted from manuscript
### http://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf

import math
from operator import attrgetter

import numpy as np

from skbio._base import SkbioObject
from skbio._config import _resolve_engine
from skbio.io.descriptors import Read, Write
from skbio.stats.distance import DistanceMatrix
from skbio.tree._exception import DuplicateNodeError, MissingNodeError

from . import _bp_cy, _bp_numba
from ._bp_cy import _BPKernel
from ._bp_numba import NUMBA_AVAILABLE


# Navigation methods implemented by the compiled kernel. ``BPTree`` binds each
# of these onto the instance at construction, so ``bp.close(i)`` resolves in the
# instance ``__dict__`` straight to the kernel's compiled method (~5 ns over a
# direct call) instead of running the class-level forwarding method (~35 ns).
# The class-level methods remain the documented API and serve unbound calls.
_KERNEL_METHODS = (
    "name",
    "length",
    "edge",
    "edge_from_number",
    "close",
    "rmq",
    "rMq",
    "depth",
    "root",
    "parent",
    "is_tip",
    "first_child",
    "last_child",
    "next_sibling",
    "previous_sibling",
    "preorder_rank",
    "preorder_select",
    "postorder_rank",
    "postorder_select",
    "is_ancestor",
    "count",
    "level_ancestor",
    "level_next",
    "lca",
    "deepest_node",
    "height",
)
_get_kernel_methods = attrgetter(*_KERNEL_METHODS)

# Compute engines of the batch operations. What ``engine="fast"`` resolves to
# depends on the build: see ``BPTree._fast_engine``.
_ENGINES = ("cython", "numba")


def _rmm_geometry(n):
    """Block size and height of the range min-max tree for ``n`` parentheses.

    The block size is ``ceil(ln(n) * ln(ln(n)))``, but at least 2. The formula
    gives less than that for ``n = 2`` (it is not even positive) and for
    ``n = 4``, where it gives 1, and the backward search does not work with
    one parenthesis per block: in ``(())`` it missed position 0, so
    ``open(2)``, and every ``lca`` through it, answered 0 instead of 1. The
    ``ceil(n / b)`` blocks are the leaves of a complete binary tree of height
    ``ceil(log2(n / b))``.
    """
    b = max(2, math.ceil(math.log(n) * math.log(math.log(n))))
    n_tip = -(-n // b)
    height = (n_tip - 1).bit_length()  # exact ceil(log2(n_tip))
    return b, height


def _build_index(B):
    """Build the navigation index of a balanced-parentheses array.

    See ``_bp_cy.build_index``.
    """
    b, height = _rmm_geometry(B.shape[0])
    return _bp_cy.build_index(B, b, height)


def _check_array(arr, dtype, name, size=None):
    """Validate a 1-D NumPy array of an exact dtype, made C-contiguous."""
    if not isinstance(arr, np.ndarray):
        raise TypeError(f"{name} must be a numpy.ndarray, not {type(arr).__name__}.")
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1-D, not {arr.ndim}-D.")
    if arr.dtype != dtype:
        raise ValueError(f"{name} must have dtype {np.dtype(dtype)}, not {arr.dtype}.")
    if size is not None and arr.shape[0] != size:
        raise ValueError(
            f"{name} must have one entry per parenthesis ({size}), not {arr.shape[0]}."
        )
    return np.ascontiguousarray(arr)


def _readonly(arr):
    arr.flags.writeable = False
    return arr


class BPTree(SkbioObject):
    """A balanced parentheses succinct data structure tree representation.

    The basis for this implementation is the data structure described by
    Cordova and Navarro [1]. In some instances, some docstring text was copied
    verbatim from the manuscript. This does not implement the bucket-based
    trees.

    A node in this data structure is represented by 2 bits, an open parenthesis
    and a close parenthesis. The implementation uses a NumPy uint8 type where
    an open parenthesis is a 1 and a close is a 0. In general, operations on
    this tree are best suited for passing in the opening parenthesis index, so
    for instance, if you'd like to use BPTree.is_tip to determine if a node is a
    leaf, the operation is defined only for using the opening parenthesis. At
    this time, there is some ambiguity over what methods can handle a closing
    parenthesis.

    Node attributes, such as names, are stored external to this data structure.

    The motivator for this data structure is pure performance both in space and
    time. As such, there is minimal sanity checking. It is advised to use this
    structure with care, and ideally within a framework which can assure
    sanity.

    Parameters
    ----------
    B : numpy.ndarray of uint8
        The parentheses bit array encoding the tree topology, where an open
        parenthesis is 1 and a close parenthesis is 0. A bool array is also
        accepted, and is viewed as uint8.
    lengths : numpy.ndarray of float64, optional
        Branch length per parenthesis (read at opening parentheses). Defaults
        to 0.
    names : numpy.ndarray of object, optional
        Node name per parenthesis (read at opening parentheses). Defaults to
        None.
    edges : numpy.ndarray of int32, optional
        Edge number per parenthesis (read at opening parentheses). Defaults to
        0, without an edge-number lookup.

    Attributes
    ----------
    data : numpy.ndarray of uint8
        The parentheses bit array encoding the tree topology, where an open
        parenthesis is 1 and a close parenthesis is 0.

    Notes
    -----
    The tree's arrays, including the navigation index built from ``data``, are
    held by this Python object. Operations on them are computed by a compiled
    (Cython) engine: per-node navigation methods (e.g., :meth:`parent`,
    :meth:`lca`) are bound directly to the engine when the tree is created, so
    they run at compiled speed.

    **The** ``fast`` **engine.** The batch operations and :meth:`cophenet`
    run the same algorithm, with bit-identical results, on either engine, so
    ``engine="fast"`` is chosen from how scikit-bio was built, rather than
    fixed:

    1. ``"cython"`` if scikit-bio was built with OpenMP, as it is with GCC on
       Linux: its multithreaded kernels are up to twice as fast as Numba's on
       ``close_batch``, ``parent_batch`` and ``level_ancestor_batch``, equal on
       :meth:`cophenet`, and compile nothing at the first call, where Numba
       compiles each kernel (from a fraction of a second to a couple of
       seconds).
    2. Otherwise ``"numba"`` if Numba is installed, as the Cython kernels then
       run on a single thread.
    3. Otherwise ``"cython"``.

    The exception, not worth a rule of its own: on the CPU, Numba is up to a
    quarter faster on :meth:`lca_batch`.

    References
    ----------
    [1] http://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    """

    default_write_format = "newick"

    # ``read`` and ``write`` are provided by the ``skbio.io`` registry via the
    # descriptor protocol, as for ``TreeNode`` and other SkbioObjects.
    read = Read()
    write = Write()

    def __init__(self, B, lengths=None, names=None, edges=None):
        # a bool array has the same one-byte 0/1 layout, which the compiled
        # engine has always accepted: view it as uint8, without a copy
        if isinstance(B, np.ndarray) and B.dtype == np.bool_:
            B = B.view(np.uint8)
        B = _check_array(B, np.uint8, "B")
        size = B.shape[0]
        if size == 0:
            raise ValueError("The topology array is empty.")

        # the tree is only valid if it is balanced (equal opens and closes)
        if int(B.sum()) * 2 != size:
            raise ValueError(
                "The topology array is unbalanced; it must contain an equal "
                "number of opening (1) and closing (0) parentheses."
            )

        if names is None:
            names = np.full(size, None, dtype=object)
        else:
            names = _check_array(names, object, "names", size)

        if lengths is None:
            lengths = np.zeros(size, dtype=np.float64)
        else:
            lengths = _check_array(lengths, np.float64, "lengths", size)

        if edges is None:
            edges = np.full(size, 0, dtype=np.int32)
            edge_lookup = None
        else:
            edges = _check_array(edges, np.int32, "edges", size)
            edge_lookup = self._edge_lookup_for(B, edges)

        index = _build_index(B)
        self._data = B
        self._size = size
        self._names = names
        self._lengths = lengths
        self._edges = edges
        self._edge_lookup = edge_lookup
        self._e_index = _readonly(index["e_index"])
        self._k_index_0 = _readonly(index["k_index_0"])
        self._k_index_1 = _readonly(index["k_index_1"])
        self._m = _readonly(index["m"])
        self._M = _readonly(index["M"])
        self._r = _readonly(index["r"])
        self._b = index["b"]
        self._height = index["height"]

        self._kernel = _BPKernel(
            B,
            self._e_index,
            self._k_index_0,
            self._k_index_1,
            self._m,
            self._M,
            self._r,
            self._b,
            self._height,
            names,
            lengths,
            edges,
            edge_lookup,
        )
        self.__dict__.update(zip(_KERNEL_METHODS, _get_kernel_methods(self._kernel)))

    @property
    def data(self):
        """The parentheses bit array (1 = open, 0 = close)."""
        return self._data

    @staticmethod
    def _edge_lookup_for(B, edges):
        opens = edges[B == 1]
        if opens.size and (opens.min() < 0 or opens.max() >= B.shape[0]):
            raise ValueError("Edge numbers must be in [0, %d)." % B.shape[0])
        return _bp_cy.edge_lookup(B, edges)

    # ------------------------------------------------------------------
    # Construction, conversion and serialization
    # ------------------------------------------------------------------

    def to_npz(self, file):
        """Save the tree to a compressed NumPy ``.npz`` archive.

        The parentheses bit array, node names, and branch lengths are stored;
        edge numbers are not. This is a lightweight binary dump, distinct from
        the registry-based :meth:`write` (which defaults to the ``newick``
        format).

        Parameters
        ----------
        file : str or file-like object
            Path or open file handle to write to.

        See Also
        --------
        from_npz
        write

        """
        np.savez_compressed(
            file, names=self._names, lengths=self._lengths, B=self._data
        )

    @classmethod
    def from_npz(cls, file):
        """Load a tree from a NumPy ``.npz`` archive written by ``to_npz``.

        Parameters
        ----------
        file : str or file-like object
            Path or open file handle to read from.

        Returns
        -------
        BPTree
            The reconstructed tree, with node names and branch lengths.

        Warnings
        --------
        This method calls :func:`numpy.load` with ``allow_pickle=True`` in
        order to restore the object-dtype ``names`` array. Loading a pickled
        array can execute arbitrary code, so only read ``.npz`` archives from
        trusted sources.

        See Also
        --------
        to_npz
        read

        """
        # names is an object array pickled by ``to_npz``, so unpickling must be
        # allowed to restore it
        data = np.load(file, allow_pickle=True)
        return cls(data["B"], names=data["names"], lengths=data["lengths"])

    @classmethod
    def from_treenode(cls, tree):
        """Construct a BPTree from a :class:`~skbio.tree.TreeNode`.

        Parameters
        ----------
        tree : skbio.tree.TreeNode
            The tree to convert.

        Returns
        -------
        BPTree
            The tree represented in balanced-parentheses form.

        See Also
        --------
        skbio.tree.TreeNode.from_bptree

        """
        topo, names, lengths, edges = _bp_cy.from_treenode_arrays(tree)
        return cls(topo, names=names, lengths=lengths, edges=edges)

    def to_array(self):
        """Return an array representation of the tree.

        This mirrors :meth:`skbio.tree.TreeNode.to_array`.

        Returns
        -------
        dict
            Dictionary with keys ``'child_index'``, ``'length'``,
            ``'id_index'`` and ``'name'``.

        See Also
        --------
        skbio.tree.TreeNode.to_array

        """
        return _bp_cy.to_array(self._kernel)

    def _to_node_arrays(self):
        """Return preorder per-node arrays for :meth:`TreeNode.from_bptree`.

        The balanced-parentheses traversal is performed in a single compiled
        pass so that :meth:`skbio.tree.TreeNode.from_bptree` pays no per-node
        Python/C call overhead; the pure-Python side is then left with only the
        ``TreeNode`` object assembly.

        Returns
        -------
        tuple of numpy.ndarray
            ``(name, length, edge, parent)``, each indexed by preorder position.
            ``name`` is dtype ``object``, ``length`` is ``float64`` and ``edge``
            is ``int32`` (mirroring :meth:`name`, :meth:`length`, :meth:`edge`).
            ``parent`` (dtype ``intp``) holds the preorder index of each node's
            parent, or ``-1`` for the root.

        See Also
        --------
        skbio.tree.TreeNode.from_bptree

        """
        return _bp_cy.to_node_arrays(self._kernel)

    def __reduce__(self):
        return (BPTree, (self._data, self._lengths, self._names))

    # ------------------------------------------------------------------
    # Node attributes
    # ------------------------------------------------------------------

    def set_names(self, names):
        """Replace the node names.

        Parameters
        ----------
        names : numpy.ndarray of object
            Node name per parenthesis (read at opening parentheses).

        """
        names = _check_array(names, object, "names", self._size)
        self._names = self._kernel._names = names

    def set_lengths(self, lengths):
        """Replace the branch lengths.

        Parameters
        ----------
        lengths : numpy.ndarray of float64
            Branch length per parenthesis (read at opening parentheses).

        """
        lengths = _check_array(lengths, np.float64, "lengths", self._size)
        self._lengths = self._kernel._lengths = lengths

    def set_edges(self, edges):
        """Replace the edge numbers and rebuild the edge-number lookup.

        Parameters
        ----------
        edges : numpy.ndarray of int32
            Edge number per parenthesis (read at opening parentheses).

        """
        edges = _check_array(edges, np.int32, "edges", self._size)
        edge_lookup = self._edge_lookup_for(self._data, edges)
        self._edges = self._kernel._edges = edges
        self._edge_lookup = self._kernel._edge_lookup = edge_lookup

    def name(self, i):
        """Name of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        str or None
            The name of node ``i``.
        """
        return self._kernel.name(i)

    def length(self, i):
        """Branch length of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        float
            The length of the branch leading to node ``i``.
        """
        return self._kernel.length(i)

    def edge(self, i):
        """Edge number of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            The edge number of the branch leading to node ``i``.
        """
        return self._kernel.edge(i)

    def edge_from_number(self, n):
        """Index of the node carrying an edge number.

        Parameters
        ----------
        n : int
            The edge number to look up.

        Returns
        -------
        int
            Index of the node whose edge number is ``n``.
        """
        return self._kernel.edge_from_number(n)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def __len__(self):
        """The number of nodes in the tree."""
        return self._size // 2

    def __str__(self):
        """Return a concise summary of the tree.

        Implements the abstract ``__str__`` required of scikit-bio objects.
        """
        return self.__repr__()

    def __repr__(self):
        """Returns summary of the tree.

        Returns
        -------
        str
            A summary of this tree

        Notes
        -----
        This method returns the name of the node and a count of tips and the
        number of internal nodes in the tree.
        """
        total_nodes = len(self)
        tip_count = self.count(tips=True)

        return "<BPTree, name: %s, internal node count: %d, tips count: %d>" % (
            self.name(0),
            total_nodes - tip_count,
            tip_count,
        )

    # ------------------------------------------------------------------
    # Navigation (computed by the compiled kernel)
    # ------------------------------------------------------------------

    def rmq(self, i, j):
        """The leftmost minimum excess in i -> j.

        Parameters
        ----------
        i : int
            Start position (inclusive).
        j : int
            End position (inclusive).

        Returns
        -------
        int
            The leftmost position in ``[i, j]`` with the minimum excess.
        """
        return self._kernel.rmq(i, j)

    def rMq(self, i, j):
        """The leftmost maximum excess in i -> j.

        Parameters
        ----------
        i : int
            Start position (inclusive).
        j : int
            End position (inclusive).

        Returns
        -------
        int
            The leftmost position in ``[i, j]`` with the maximum excess.
        """
        return self._kernel.rMq(i, j)

    def close(self, i):
        """The position of the closing parenthesis that matches B[i].

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Position of the matching closing parenthesis.
        """
        return self._kernel.close(i)

    def depth(self, i):
        """The depth of given node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Depth of node relative to the root of tree.
        """
        return self._kernel.depth(i)

    def root(self):
        """The index of the root node of the tree."""
        return self._kernel.root()

    def parent(self, i):
        """The parent of node.

        Parameters
        ----------
        i : int
            Index of node to evaluate.

        Returns
        -------
        int
            Index of parent node. Returns -1 if node does not have a parent.
        """
        return self._kernel.parent(i)

    def is_tip(self, i):
        """Whether the node is a tip of a tree.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        bool
            Whether the node is a tip of a tree or not.
        """
        return self._kernel.is_tip(i)

    def first_child(self, i):
        """Index of the first (leftmost) child of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the first child of node ``i``, or 0 if ``i`` is a tip
            (0 is the root, which can never be a child).

        See Also
        --------
        last_child
        skbio.tree.TreeNode

        Notes
        -----
        Returns an integer index into the parentheses bit array, not a node
        object. Corresponds to accessing ``children[0]`` on a
        :class:`~skbio.tree.TreeNode`.
        """
        return self._kernel.first_child(i)

    def last_child(self, i):
        """Index of the last (rightmost) child of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the last child of node ``i``, or 0 if ``i`` is a tip
            (0 is the root, which can never be a child).

        See Also
        --------
        first_child
        skbio.tree.TreeNode

        Notes
        -----
        Returns an integer index into the parentheses bit array, not a node
        object. Corresponds to accessing ``children[-1]`` on a
        :class:`~skbio.tree.TreeNode`.
        """
        return self._kernel.last_child(i)

    def mincount(self, i, j):
        """Number of occurrences of the minimum in excess(i), ..., excess(j)."""
        excess, counts = np.unique(self._e_index[i : j + 1], return_counts=True)
        return counts[excess.argmin()]

    def minselect(self, i, j, q):
        """Position of the qth minimum in excess(i), ..., excess(j).

        Parameters
        ----------
        i, j : int
            The range of positions, ``i <= j``.
        q : int
            Which occurrence of the minimum, counting from 1.

        Returns
        -------
        int or None
            The position of the ``q``-th occurrence of the minimum excess in
            the range, or None if there is no such occurrence: fewer than ``q``
            of them, or ``q < 1``.
        """
        if q < 1:
            # ranks count from 1; a lower q would index from the end
            return None
        counts = self._e_index[i : j + 1]
        index = counts == counts.min()

        if index.sum() < q:
            return None
        else:
            return i + index.nonzero()[0][q - 1]

    def next_sibling(self, i):
        """Index of the next (right) sibling of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the next sibling of node ``i``, or 0 if ``i`` has no
            next sibling (0 is the root, which can never be a sibling).

        See Also
        --------
        previous_sibling
        skbio.tree.TreeNode.siblings

        Notes
        -----
        Returns the integer index of a single sibling, unlike
        :meth:`~skbio.tree.TreeNode.siblings`, which returns a list of all
        sibling nodes.
        """
        return self._kernel.next_sibling(i)

    def previous_sibling(self, i):
        """Index of the previous (left) sibling of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the previous sibling of node ``i``, or 0 if ``i`` has no
            previous sibling (0 is the root, which can never be a sibling).

        See Also
        --------
        next_sibling
        skbio.tree.TreeNode.siblings

        Notes
        -----
        Returns the integer index of a single sibling, unlike
        :meth:`~skbio.tree.TreeNode.siblings`, which returns a list of all
        sibling nodes.
        """
        return self._kernel.previous_sibling(i)

    def preorder_rank(self, i):
        """Preorder rank of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            The position of node ``i`` in a preorder traversal of the tree.

        See Also
        --------
        preorder_select
        skbio.tree.TreeNode.preorder

        Notes
        -----
        Returns the node's integer position in preorder, not a node object.
        This differs from :meth:`~skbio.tree.TreeNode.preorder`, which yields
        the nodes of the tree in preorder. The inverse of
        :meth:`preorder_select`.
        """
        return self._kernel.preorder_rank(i)

    def preorder_select(self, k):
        """Index of the node with a given preorder rank.

        Parameters
        ----------
        k : int
            Preorder rank to look up.

        Returns
        -------
        int
            Index of the node whose preorder rank is ``k``.

        See Also
        --------
        preorder_rank
        skbio.tree.TreeNode.preorder

        Notes
        -----
        The inverse of :meth:`preorder_rank`. Returns an integer index into
        the parentheses bit array, not a node object.
        :meth:`~skbio.tree.TreeNode.preorder` yields the nodes in this order.
        """
        return self._kernel.preorder_select(k)

    def postorder_rank(self, i):
        """Postorder rank of a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            The position of node ``i`` in a postorder traversal of the tree.

        See Also
        --------
        postorder_select
        skbio.tree.TreeNode.postorder

        Notes
        -----
        Returns the node's integer position in postorder, not a node object.
        This differs from :meth:`~skbio.tree.TreeNode.postorder`, which yields
        the nodes of the tree in postorder. The inverse of
        :meth:`postorder_select`.
        """
        return self._kernel.postorder_rank(i)

    def postorder_select(self, k):
        """Index of the node with a given postorder rank.

        Parameters
        ----------
        k : int
            Postorder rank to look up.

        Returns
        -------
        int
            Index of the node whose postorder rank is ``k``.

        See Also
        --------
        postorder_rank
        skbio.tree.TreeNode.postorder

        Notes
        -----
        The inverse of :meth:`postorder_rank`. Returns an integer index into
        the parentheses bit array, not a node object.
        :meth:`~skbio.tree.TreeNode.postorder` yields the nodes in this order.
        """
        return self._kernel.postorder_select(k)

    def is_ancestor(self, i, j):
        """Whether a node is an ancestor of another node.

        Parameters
        ----------
        i : int
            A node index
        j : int
            A node index

        Note
        ----
        False is returned if i == j. A node cannot be an ancestor of itself.

        Returns
        -------
        bool
            True if i is an ancestor of j, False otherwise.
        """
        return self._kernel.is_ancestor(i, j)

    def count(self, i=0, tips=False):
        """Get the count of nodes in the subtree rooted at a node.

        Parameters
        ----------
        i : int, optional
            Index of the node whose subtree is evaluated. Defaults to the root
            (``0``), i.e., the whole tree.
        tips : bool, optional
            If True, only count the tips (leaves) in the subtree (default:
            False).

        Returns
        -------
        int
            The number of nodes (or tips, if ``tips`` is True) in the subtree
            rooted at node ``i``, including ``i`` itself.

        See Also
        --------
        skbio.tree.TreeNode.count

        Notes
        -----
        Returns a node count, not a subtree.

        """
        return self._kernel.count(i, tips)

    def level_ancestor(self, i, d):
        """Index of the ancestor a given number of levels above a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.
        d : int
            Number of levels to ascend toward the root.

        Returns
        -------
        int
            Index of the ancestor ``d`` levels above node ``i``, or -1 if
            ``d`` is not positive.

        See Also
        --------
        level_next
        skbio.tree.TreeNode.ancestors

        Notes
        -----
        Returns an integer index into the parentheses bit array, not a node
        object. :meth:`~skbio.tree.TreeNode.ancestors` returns the full list
        of ancestor nodes from a node toward the root.
        """
        return self._kernel.level_ancestor(i, d)

    def level_next(self, i):
        """Index of the next node at the same depth.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the next node at the same depth as node ``i``, or -1 if
            there is no such node.

        See Also
        --------
        level_ancestor
        skbio.tree.TreeNode.levelorder

        Notes
        -----
        Returns an integer index into the parentheses bit array, not a node
        object. :meth:`~skbio.tree.TreeNode.levelorder` traverses all nodes
        depth by depth.
        """
        return self._kernel.level_next(i)

    def lca(self, i, j):
        """The lowest common ancestor of two nodes.

        Parameters
        ----------
        i : int
            A node index to evaluate
        j : int
            A node index to evaluate

        Returns
        -------
        int
           The index of the lowest common ancestor. A node is its own lowest
           common ancestor.
        """
        return self._kernel.lca(i, j)

    def deepest_node(self, i):
        """Index of the deepest node descending from a node.

        Parameters
        ----------
        i : int
            Index of the node to evaluate.

        Returns
        -------
        int
            Index of the deepest (most distant) tip descending from node ``i``.

        See Also
        --------
        height
        skbio.tree.TreeNode.height

        Notes
        -----
        Returns an integer index into the parentheses bit array, not a node
        object. This is the tip that :meth:`~skbio.tree.TreeNode.height`
        returns as the second element of its ``(height, tip)`` result.
        """
        return self._kernel.deepest_node(i)

    def height(self, i):
        """The height of node i with respect to its deepest descendent

        Parameters
        ----------
        i : int
            The node to evaluate

        Notes
        -----
        Height is in terms of number of edges, not in terms of branch length

        Returns
        -------
        int
            The number of edges between node i and its deepest node
        """
        return self._kernel.height(i)

    # ------------------------------------------------------------------
    # Batch operations (computed by the selected engine)
    # ------------------------------------------------------------------

    def _fast_engine(self):
        """What ``engine="fast"`` resolves to for this tree.

        See the Notes of :class:`BPTree` for the rules and their measurements.
        """
        if _bp_cy.OPENMP or not NUMBA_AVAILABLE:
            return "cython"  # multithreaded Cython, or nothing else installed
        return "numba"  # Cython runs serially here; Numba runs in parallel

    def _run(self, kernel, engine, *args):
        """Run a batch kernel of the selected compute engine on this tree."""
        engine = _resolve_engine(engine, _ENGINES, fast=self._fast_engine())
        if engine == "numba":
            return getattr(_bp_numba, kernel)(_bp_numba.bp_arrays(self), *args)
        return getattr(_bp_cy, kernel)(self._kernel, *args)

    def _positions(self, *arrays):
        """Validate node positions: broadcast, flatten, check the range."""
        arrays = np.broadcast_arrays(*(np.asarray(a) for a in arrays))
        shape = arrays[0].shape
        out = []
        for arr in arrays:
            # an empty list is float64 but holds no positions to check
            if arr.size and arr.dtype.kind not in "iu":
                raise TypeError("Node positions must be integers, not %s." % arr.dtype)
            arr = np.ascontiguousarray(arr.ravel(), dtype=np.intp)
            if arr.size and (arr.min() < 0 or arr.max() >= self._size):
                raise IndexError("Node positions must be in [0, %d)." % self._size)
            out.append(arr)
        return shape, out

    def close_batch(self, i, engine=None):
        """Matching closing parenthesis of each position in a batch.

        The batch counterpart of :meth:`close`: one independent query per
        element of ``i``, computed by the selected engine.

        Parameters
        ----------
        i : array_like of int
            Positions to evaluate.
        engine : {'cython', 'numba', 'fast'}, optional
            The :ref:`compute engine <compute_engines>`. Defaults to the global
            ``compute_engine`` option. ``'fast'`` is chosen per tree: see the
            Notes of :class:`BPTree`.

        Returns
        -------
        numpy.ndarray of intp
            ``close(i)`` for each position, in the shape of ``i``.

        See Also
        --------
        close

        """
        shape, (i,) = self._positions(i)
        return self._run("close_batch", engine, i).reshape(shape)

    def parent_batch(self, i, engine=None):
        """Parent of each node in a batch.

        The batch counterpart of :meth:`parent`: one independent query per
        element of ``i``, computed by the selected engine.

        Parameters
        ----------
        i : array_like of int
            Node positions to evaluate.
        engine : {'cython', 'numba', 'fast'}, optional
            The :ref:`compute engine <compute_engines>`. Defaults to the global
            ``compute_engine`` option. ``'fast'`` is chosen per tree: see the
            Notes of :class:`BPTree`.

        Returns
        -------
        numpy.ndarray of intp
            ``parent(i)`` for each position (-1 for the root), in the shape of
            ``i``.

        See Also
        --------
        parent

        """
        shape, (i,) = self._positions(i)
        return self._run("parent_batch", engine, i).reshape(shape)

    def lca_batch(self, i, j, engine=None):
        """Lowest common ancestor of each pair of nodes in a batch.

        The batch counterpart of :meth:`lca`: one independent query per pair
        ``(i[k], j[k])``, computed by the selected engine. This is not the
        lowest common ancestor of a set of nodes (as
        :meth:`skbio.tree.TreeNode.lca` computes).

        Parameters
        ----------
        i, j : array_like of int
            Node positions of each pair, broadcast against each other. The
            order within a pair does not matter.
        engine : {'cython', 'numba', 'fast'}, optional
            The :ref:`compute engine <compute_engines>`. Defaults to the global
            ``compute_engine`` option. ``'fast'`` is chosen per tree: see the
            Notes of :class:`BPTree`.

        Returns
        -------
        numpy.ndarray of intp
            The position of the lowest common ancestor of each pair, in the
            broadcast shape of ``i`` and ``j``.

        See Also
        --------
        lca

        Notes
        -----
        Unlike :meth:`lca`, which requires ``i <= j``, each pair is ordered
        before the query.

        """
        shape, (i, j) = self._positions(i, j)
        return self._run("lca_batch", engine, i, j).reshape(shape)

    def level_ancestor_batch(self, i, d, engine=None):
        """Ancestor a given number of levels above each node in a batch.

        The batch counterpart of :meth:`level_ancestor`: one independent query
        per element of ``i`` (and ``d``), computed by the selected engine.

        Parameters
        ----------
        i : array_like of int
            Node positions to evaluate.
        d : int or array_like of int
            Number of levels to ascend, broadcast against ``i``.
        engine : {'cython', 'numba', 'fast'}, optional
            The :ref:`compute engine <compute_engines>`. Defaults to the global
            ``compute_engine`` option. ``'fast'`` is chosen per tree: see the
            Notes of :class:`BPTree`.

        Returns
        -------
        numpy.ndarray of intp
            ``level_ancestor(i, d)`` for each pair, in the broadcast shape of
            ``i`` and ``d``.

        See Also
        --------
        level_ancestor

        """
        i, d = np.broadcast_arrays(np.asarray(i), np.asarray(d))
        if d.size and d.dtype.kind not in "iu":
            raise TypeError("Levels must be integers, not %s." % d.dtype)
        shape, (i,) = self._positions(i)
        d = np.ascontiguousarray(d.ravel(), dtype=np.intp)
        return self._run("level_ancestor_batch", engine, i, d).reshape(shape)

    def cophenet(self, endpoints=None, use_length=True, engine=None):
        r"""Return a distance matrix between each pair of tips in the tree.

        The balanced-parentheses counterpart of
        :meth:`skbio.tree.TreeNode.cophenet`.

        Parameters
        ----------
        endpoints : iterable of str, optional
            Names of the tips to be included in the calculation. The returned
            distance matrix will use this order. If not specified, all named
            tips will be included, in the order they appear in the tree.
        use_length : bool, optional
            Whether to return the sum of branch lengths (True, default) or the
            number of branches (False) connecting each pair of tips.
        engine : {'cython', 'numba', 'fast'}, optional
            The :ref:`compute engine <compute_engines>`. Defaults to the global
            ``compute_engine`` option. ``'fast'`` is chosen per tree: see the
            Notes of :class:`BPTree`.

        Returns
        -------
        DistanceMatrix
            The cophenetic distance matrix.

        Raises
        ------
        MissingNodeError
            If any of the specified ``endpoints`` are not found in the tree.
        DuplicateNodeError
            If the specified ``endpoints`` have duplicates, or if the tree has
            duplicate tip names when ``endpoints`` is not specified.
        ValueError
            If any of the specified ``endpoints`` are not tips.

        See Also
        --------
        lca_batch
        skbio.tree.TreeNode.cophenet

        Notes
        -----
        The distance between tips :math:`a` and :math:`b` with lowest common
        ancestor :math:`c` is :math:`d(a) + d(b) - 2d(c)`, where :math:`d` is
        the sum of branch lengths (or the number of branches) from the root.
        The tip pairs are independent and are computed in parallel by the
        selected engine. Missing branch lengths are 0.

        Examples
        --------
        >>> from skbio.tree import BPTree
        >>> tree = BPTree.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        >>> print(tree.cophenet())
        4x4 distance matrix
        IDs:
        'a', 'b', 'd', 'e'
        Data:
        [[  0.   3.  14.  15.]
         [  3.   0.  15.  16.]
         [ 14.  15.   0.   9.]
         [ 15.  16.   9.   0.]]

        """
        B = self._data
        is_tip = np.zeros(self._size, dtype=bool)
        is_tip[:-1] = (B[:-1] == 1) & (B[1:] == 0)
        tips = np.flatnonzero(is_tip)

        if not endpoints:
            names = self._names[tips]
            named = np.array([name is not None for name in names], dtype=bool)
            tips = tips[named]
            taxa = names[named].tolist()
            if len(set(taxa)) < len(taxa):
                raise DuplicateNodeError("Tree contains duplicate tip names.")
        else:
            # name lookup as TreeNode.find: tips first, then the first internal
            # node (in preorder) with the name
            lookup = {}
            for pos in tips.tolist():
                lookup.setdefault(self._names[pos], pos)
            for pos in np.flatnonzero(B & ~is_tip).tolist():
                lookup.setdefault(self._names[pos], pos)
            taxa, positions = [], []
            for name in endpoints:
                if name is None:
                    raise MissingNodeError("Cannot find a node without a name.")
                if name in taxa:
                    raise DuplicateNodeError(f"Duplicate tip name '{name}' found.")
                pos = lookup.get(name)
                if pos is None:
                    raise MissingNodeError(f"Node '{name}' is not found in the tree.")
                if not is_tip[pos]:
                    raise ValueError(f"Node with name '{name}' is not a tip.")
                taxa.append(name)
                positions.append(pos)
            tips = np.array(positions, dtype=np.intp)

        tips = np.ascontiguousarray(tips, dtype=np.intp)
        if tips.size < 2:
            # no pairs; an empty condensed vector would expand to 1 x 1
            return DistanceMatrix(
                np.zeros((tips.size, tips.size)), taxa, validate=False
            )
        engine = _resolve_engine(engine, _ENGINES, fast=self._fast_engine())

        # the kernel takes the tips in tree order, and the output row (slot) of
        # each: its place in the requested order
        order = np.argsort(tips, kind="stable")
        sorted_tips = tips[order]
        slot = np.ascontiguousarray(order, dtype=np.intp)

        # per node: its parent, and the end of its tips among the sorted tips
        # (the number of them before its closing parenthesis)
        opens = np.flatnonzero(B)
        parent = np.full(self._size, -1, dtype=np.intp)
        parent[opens] = self._run("parent_batch", engine, opens)
        selected = np.zeros(self._size, dtype=np.intp)
        selected[sorted_tips] = 1
        end = np.zeros(self._size, dtype=np.intp)
        end[opens] = np.cumsum(selected)[self._run("close_batch", engine, opens)]
        del selected

        if use_length:
            dist = self._run("root_distances", engine, self._lengths)
        else:
            dist = self._e_index.astype(np.float64)
        kernels = _bp_numba if engine == "numba" else _bp_cy
        matrix = kernels.tip_distances(sorted_tips, slot, parent, end, dist)
        return DistanceMatrix(matrix, taxa, validate=False)

    # ------------------------------------------------------------------
    # Whole-tree operations
    # ------------------------------------------------------------------

    def shear(self, tips):
        """Remove all nodes from the tree except tips and ancestors of tips.

        Parameters
        ----------
        tips : set of str
            The set of tip names to retain

        Returns
        -------
        BPTree
            A new BPTree corresponding to only the described tips and their
            ancestors.
        """
        if not isinstance(tips, set):
            raise TypeError("tips must be a set, not %s." % type(tips).__name__)
        mask, count = _bp_cy.shear_mask(self._kernel, tips)
        if count == 0:
            raise ValueError("No requested tips found")
        return self._from_mask(mask, self._lengths)

    def collapse(self):
        """Collapse single-child internal nodes.

        Every internal node with exactly one child is removed from the tree,
        and the removed node's branch length is added to that of its single
        child so that root-to-tip path lengths are preserved. The root and all
        tips are always retained, as are internal nodes with two or more
        children.

        Returns
        -------
        BPTree
            A new tree with all single-child internal nodes removed. Node names
            and the merged branch lengths are carried over; edge numbers are
            not retained.

        Notes
        -----
        A new ``BPTree`` is returned; the original tree is not modified.

        """
        mask, lengths = _bp_cy.collapse_mask(self._kernel)
        return self._from_mask(mask, lengths)

    def _from_mask(self, mask, lengths):
        """A new tree of the positions set in ``mask``."""
        keep = mask.view(bool)
        return BPTree(self._data[keep], names=self._names[keep], lengths=lengths[keep])
