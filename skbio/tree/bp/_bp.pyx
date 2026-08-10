# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
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

from libc.math cimport ceil, log as ln, pow, log2
import numpy as np
cimport numpy as cnp
cimport cython

from ._bp_binary_tree cimport *
from ._ba cimport *

cnp.import_array()

cdef extern from "limits.h":
    int INT_MAX


DOUBLE = np.float64
SIZE = np.intp
BOOL = np.uint8
INT32 = np.int32


cdef inline int min(int a, int b) nogil:
    if a > b:
        return b
    else:
        return a


cdef inline int max(int a, int b) nogil:
    if a > b:
        return a
    else:
        return b


cdef class mM:
    def __cinit__(self, BOOL_t[:] B, int B_size):
        self.m_idx = 0
        self.M_idx = 1

        self.rmm(B, B_size)

    cdef void rmm(self, BOOL_t[:] B, int B_size) nogil:
        """Construct the rmM tree based off of Navarro and Sadakane

        http://www.dcc.uchile.cl/~gnavarro/ps/talg12.pdf
        """
        cdef int i, j, lvl, pos  # for loop support
        cdef int offset  # tip offset in binary tree for a given parenthesis
        cdef int lower_limit  # the lower limit of the bucket a parenthesis is in
        cdef int upper_limit  # the upper limit of the bucket a parenthesis is in
        cdef int min_ = 0 # m, absolute minimum for a blokc
        cdef int max_ = 0 # M, absolute maximum for a block
        cdef int excess = 0 # e, absolute excess
        cdef int vbar
        cdef int r = 0

        # build tip info
        self.b = <int>ceil(ln(<double> B_size) * ln(ln(<double> B_size)))

        # determine the number of nodes and height of the binary tree
        self.n_tip = <int>ceil(B_size / <double> self.b)
        self.height = <int>ceil(log2(self.n_tip))
        self.n_internal = <int>(pow(2, self.height)) - 1
        self.n_total = self.n_tip + self.n_internal

        with gil:
            # creation of a memoryview directly or via numpy requires the GIL:
            # http://stackoverflow.com/a/22238012
            self.mM = np.zeros((self.n_total, 2), dtype=SIZE)
            self.r = np.zeros(self.n_total, dtype=SIZE)

        # annoying, cannot do step in range if step is not known at runtime
        # see https://github.com/cython/cython/pull/520
        # for i in range(0, B_size, b):
        # as a result, doing a custom range using a while loop
        # compute for tips of rmM tree
        i = 0
        while i < B_size:
            offset = i // self.b
            lower_limit = i
            upper_limit = min(i + self.b, B_size)
            min_ = INT_MAX
            max_ = 0

            self.r[offset + self.n_internal] = r
            for j in range(lower_limit, upper_limit):
                # G function, a +-1 method where if B[j] == 1 we +1, and if
                # B[j] == 0 we -1
                excess += -1 + (2 * B[j])
                r += B[j]

                if excess < min_:
                    min_ = excess

                if excess > max_:
                    max_ = excess

                # at the left bound of the bucket

            self.mM[offset + self.n_internal, self.m_idx] = min_
            self.mM[offset + self.n_internal, self.M_idx] = max_

            i += self.b

        # compute for internal nodes of rmM tree in reverse level order starting
        # at the level above the tips
        for lvl in range(self.height - 1, -1, -1):
            num_curr_nodes = <int>pow(2, lvl)

            # for each node in the level
            for pos in range(num_curr_nodes):
                # obtain the node, and the index to its children
                node = bt_node_from_left(pos, lvl)
                lchild = bt_left_child(node)
                rchild = bt_right_child(node)

                if lchild >= self.n_total:
                    continue

                elif rchild >= self.n_total:
                    self.mM[node, self.m_idx] = self.mM[lchild, self.m_idx]
                    self.mM[node, self.M_idx] = self.mM[lchild, self.M_idx]
                else:
                    self.mM[node, self.m_idx] = min(self.mM[lchild, self.m_idx],
                                                    self.mM[rchild, self.m_idx])
                    self.mM[node, self.M_idx] = max(self.mM[lchild, self.M_idx],
                                                    self.mM[rchild, self.M_idx])

                self.r[node] = self.r[lchild]


@cython.final
cdef class BPTree:
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

    Attributes
    ----------
    data : numpy.ndarray of uint8
        The parentheses bit array encoding the tree topology, where an open
        parenthesis is 1 and a close parenthesis is 0.

    References
    ----------
    [1] http://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    """

    def __cinit__(self, cnp.ndarray[BOOL_t, ndim=1] B,
                  cnp.ndarray[DOUBLE_t, ndim=1] lengths=None,
                  cnp.ndarray[object, ndim=1] names=None,
                  cnp.ndarray[INT32_t, ndim=1] edges=None):
        cdef SIZE_t i
        cdef SIZE_t size
        cdef SIZE_t[:] _e_index
        cdef SIZE_t[:] _k_index_0
        cdef SIZE_t[:] _k_index_1
        cdef SIZE_t[:] _r_index_0
        cdef SIZE_t[:] _r_index_1
        cdef cnp.ndarray[object, ndim=1] _names
        cdef cnp.ndarray[DOUBLE_t, ndim=1] _lengths
        cdef cnp.ndarray[INT32_t, ndim=1] _edges
        cdef cnp.ndarray[SIZE_t, ndim=1] _edge_lookup

        # the tree is only valid if it is balanced (equal opens and closes)
        if B.sum() * 2 != B.size:
            raise ValueError(
                "The topology array is unbalanced; it must contain an equal "
                "number of opening (1) and closing (0) parentheses."
            )
        self.data = B
        self._b_ptr = &B[0]
        self.size = B.size

        self._rmm = mM(B, B.size)

        if names is not None:
            self._names = names
        else:
            self._names = np.full(self.data.size, None, dtype=object)

        if lengths is not None:
            self._lengths = lengths
        else:
            self._lengths = np.zeros(self.data.size, dtype=DOUBLE)

        if edges is not None:
            self._set_edges(edges)
        else:
            self._edges = np.full(self.data.size, 0, dtype=INT32)
            self._edge_lookup = None

        # precursor for select index cache
        _r_index_0 = np.cumsum((1 - B), dtype=SIZE)
        _r_index_1 = np.cumsum(B, dtype=SIZE)

        # construct a select index. These operations are performed frequently,
        # and easy to cache at a relatively minor memory expense. It cannot be
        # assumed that open and close will be same length so can't stack
        _k_index_0 = np.unique(_r_index_0,
                               return_index=True)[1].astype(SIZE)
        self._k_index_0 = _k_index_0
        _k_index_1 = np.unique(_r_index_1,
                               return_index=True)[1].astype(SIZE)
        self._k_index_1 = _k_index_1

        # construct an excess index. These operations are performed a lot, and
        # similarly can to rank and select, can be cached at a minimal expense.
        _e_index = np.empty(B.size, dtype=SIZE)
        for i in range(B.size):
            _e_index[i] = self._excess(i)
        self._e_index = _e_index

    default_write_format = "newick"

    def write(self, object file, format=None, **kwargs):
        """Write the tree to a file via the ``skbio.io`` registry.

        Parameters
        ----------
        file : str or file-like object
            Path or open file handle to write to.
        format : str, optional
            The file format to write. Defaults to ``"newick"``.
        kwargs : dict, optional
            Format-specific parameters passed to the writer.

        See Also
        --------
        read
        skbio.io.registry.write
        skbio.io.format.newick

        """
        # imported lazily to avoid an import cycle at module load time
        import skbio.io

        return skbio.io.write(
            self, into=file, format=format or "newick", **kwargs
        )

    @staticmethod
    def read(object file, format=None, **kwargs):
        """Read a tree from a file via the ``skbio.io`` registry.

        Parameters
        ----------
        file : str or file-like object
            Path or open file handle to read from.
        format : str, optional
            The file format to read. If omitted, the format is inferred.
        kwargs : dict, optional
            Format-specific parameters passed to the reader.

        Returns
        -------
        BPTree
            The parsed tree.

        See Also
        --------
        write
        skbio.io.registry.read
        skbio.io.format.newick

        """
        # imported lazily to avoid an import cycle at module load time
        import skbio.io

        return skbio.io.read(file, into=BPTree, format=format, **kwargs)

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
        cdef:
            Py_ssize_t n_nodes, ptr
            cnp.ndarray[BOOL_t, ndim=1] topo
            cnp.ndarray[object, ndim=1] names
            cnp.ndarray[DOUBLE_t, ndim=1] lengths
            cnp.ndarray[INT32_t, ndim=1] edges
            set seen

        n_nodes = len(list(tree.traverse(include_self=True)))

        topo = np.zeros(n_nodes * 2, dtype=np.uint8)
        names = np.full(n_nodes * 2, None, dtype=object)
        lengths = np.zeros(n_nodes * 2, dtype=np.double)
        edges = np.zeros(n_nodes * 2, dtype=np.int32)

        ptr = 0
        seen = set()
        for n in tree.pre_and_postorder(include_self=True):
            if n not in seen:
                topo[ptr] = 1
                names[ptr] = n.name
                lengths[ptr] = n.length or 0.0
                edges[ptr] = getattr(n, 'edge_num', None) or 0

                if n.is_tip():
                    ptr += 1

                seen.add(n)

            ptr += 1
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
        cdef:
            Py_ssize_t i, n
            Py_ssize_t chi_ptr, cur_index
            SIZE_t node_idx, first_child, last_child, sib_idx
            cnp.ndarray[DOUBLE_t, ndim=1] length
            cnp.ndarray[UINT32_t, ndim=1] node_ids
            cnp.ndarray[object, ndim=1] name
            dict id_index

        class mock_node:
            def __init__(self, id, is_tip):
                self.is_tip_ = is_tip
                self.id = id

            def is_tip(self):
                return self.is_tip_

        n = <Py_ssize_t>self.data.sum()

        child_index = np.zeros((n - self.count(tips=True), 3), dtype=np.int64)
        length = np.zeros(n, dtype=np.double)
        node_ids = np.zeros(self.data.size, dtype=np.uint32)
        name = np.full(n, None, dtype=object)

        # TreeNode.assign_ids, decompose target
        chi_ptr = 0
        cur_index = 0  # the index into node_ids, equivalent to TreeNode.assign_ids
        id_index = dict.fromkeys(set(range(n)))  # map a node's "id" to an object which indicates if it is a leaf or not
        for i in range(n):
            node_idx = self.postorder_select(i + 1)  # the index within the BP of the node

            if not self.is_tip(node_idx):
                first_child = self.first_child(node_idx)
                last_child = self.last_child(node_idx)

                sib_idx = first_child  # the sibling index wtihin the BP of the node
                while sib_idx != 0 and sib_idx <= last_child:
                    node_ids[sib_idx] = cur_index
                    id_index[cur_index] = mock_node(cur_index, self.is_tip(sib_idx))
                    length[cur_index] = self.length(sib_idx)
                    name[cur_index] = self.name(sib_idx)

                    cur_index += 1
                    sib_idx = self.next_sibling(sib_idx)

                child_index[chi_ptr] = [node_idx, node_ids[first_child], node_ids[last_child]]
                chi_ptr += 1

        # make sure to capture root
        id_index[n - 1] = mock_node(cur_index, False)

        node_ids[0] = cur_index
        child_index[:, 0] = node_ids[child_index[:, 0]]
        child_index = child_index[np.argsort(child_index[:, 0])]

        return {'child_index': child_index, 'length': length, 'id_index': id_index,
                'name': name}

    def set_names(self, cnp.ndarray[object, ndim=1] names):
        self._names = names

    def set_lengths(self, cnp.ndarray[DOUBLE_t, ndim=1] lengths):
        self._lengths = lengths

    cdef void _set_edges(self, cnp.ndarray[INT32_t, ndim=1] edges):
        cdef:
            int i, n
            INT32_t edge
            cnp.ndarray[SIZE_t, ndim=1] _edge_lookup
            cnp.ndarray[BOOL_t, ndim=1] b

        b = self.data
        n = b.size

        _edge_lookup = np.full(n, 0, dtype=SIZE)
        for i in range(n):
            if b[i] == 1:
                edge = edges[i]
                _edge_lookup[edge] = i

        self._edge_lookup = _edge_lookup
        self._edges = edges

    def set_edges(self, cnp.ndarray[INT32_t, ndim=1] edges):
        self._set_edges(edges)

    cpdef inline unicode name(self, SIZE_t i):
        return self._names[i]

    cpdef inline DOUBLE_t length(self, SIZE_t i):
        return self._lengths[i]

    cpdef inline INT32_t edge(self, SIZE_t i):
        return self._edges[i]

    cpdef SIZE_t edge_from_number(self, INT32_t n):
        return self._edge_lookup[n]

    cdef inline SIZE_t rank(self, SIZE_t t, SIZE_t i) nogil:
        """Determine the rank order of the ith bit t

        Rank is the order of the ith bit observed, from left to right. For
        t=1, this is a preorder traversal of the tree.

        Parameters
        ----------
        t : SIZE_t
            The bit value, either 0 or 1 where 0 is a closing parenthesis and
            1 is an opening.
        i : SIZE_T
            The position to evaluate

        Returns
        -------
        SIZE_t
            The rank order of the position.
        """
        cdef int k
        cdef int r = 0
        cdef int lower_bound
        cdef int upper_bound
        cdef int j
        cdef int node

        k = i // self._rmm.b

        lower_bound = k * self._rmm.b

        # upper_bound is block boundary or end of tree
        upper_bound = min((k + 1) * self._rmm.b, self.size)
        upper_bound = min(upper_bound, i + 1)

        # collect rank from within the block
        for j in range(lower_bound, upper_bound):
            r += self._b_ptr[j]

        # collect the rank at the left end of the block
        node = bt_node_from_left(k, self._rmm.height)
        r += self._rmm.r[node]

        if t:
            return r
        else:
            return (i - r) + 1

    cdef inline SIZE_t select(self, SIZE_t t, SIZE_t k) nogil:
        """The position in B of the kth occurrence of the bit t."""
        if t:
            return self._k_index_1[k]
        else:
            return self._k_index_0[k]

    cdef SIZE_t _excess(self, SIZE_t i) nogil:
        """Actually compute excess"""
        if i < 0:
            return 0  # wasn't stated as needed but appears so given testing
        return (2 * self.rank(1, i) - i) - 1

    cdef SIZE_t excess(self, SIZE_t i) nogil:
        """the number of opening minus closing parentheses in B[1, i]"""
        # same as: self.rank(1, i) - self.rank(0, i)
        return self._e_index[i]

    cpdef inline SIZE_t close(self, SIZE_t i) nogil:
        """The position of the closing parenthesis that matches B[i]"""
        if not self._b_ptr[i]:
            # identity: the close of a closed parenthesis is itself
            return i

        return self.fwdsearch(i, -1)

    cdef inline SIZE_t open(self, SIZE_t i) nogil:
        """The position of the opening parenthesis that matches B[i]"""
        if self._b_ptr[i] or i <= 0:
            # identity: the open of an open parenthesis is itself
            # the open of 0 is open. A negative index cannot be open, so just return
            return i

        return self.bwdsearch(i, 0) + 1

    cdef inline SIZE_t enclose(self, SIZE_t i) nogil:
        """The opening parenthesis of the smallest matching pair that contains position i"""
        if self._b_ptr[i]:
            return self.bwdsearch(i, -2) + 1
        else:
            return self.bwdsearch(i - 1, -2) + 1

    cpdef SIZE_t rmq(self, SIZE_t i, SIZE_t j) nogil:
        """The leftmost minimum excess in i -> j"""
        cdef:
            SIZE_t k, min_k
            SIZE_t min_v, obs_v

        min_k = i
        min_v = self.excess(i)  # a value larger than what will be tested
        for k in range(i, j + 1):
            obs_v = self.excess(k)
            if obs_v < min_v:
                min_k = k
                min_v = obs_v
        return min_k

    cpdef SIZE_t rMq(self, SIZE_t i, SIZE_t j) nogil:
        """The leftmost maximmum excess in i -> j."""
        cdef:
            SIZE_t k, max_k
            SIZE_t max_v, obs_v

        max_k = i
        max_v = self.excess(i)  # a value larger than what will be tested
        for k in range(i, j + 1):
            obs_v = self.excess(k)
            if obs_v > max_v:
                max_k = k
                max_v = obs_v

        return max_k

    def __len__(self):
        """The number of nodes in the tree."""
        return self.size / 2

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
        cdef total_nodes = len(self)
        cdef tip_count = self.count(tips=True)

        return "<BPTree, name: %s, internal node count: %d, tips count: %d>" % \
                (self.name(0), total_nodes - tip_count, tip_count)

    def __reduce__(self):
        return (BPTree, (self.data, self._lengths, self._names))

    cpdef SIZE_t depth(self, SIZE_t i) nogil:
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
        return self._e_index[i]

    cpdef SIZE_t root(self) nogil:
        """The index of the root node of the tree."""
        return 0

    cpdef SIZE_t parent(self, SIZE_t i) nogil:
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
        if i == self.root() or i == (self.size - 1):
            return -1
        else:
            return self.enclose(i)

    cpdef BOOL_t is_tip(self, SIZE_t i) nogil:
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
        return self._b_ptr[i] and (not self._b_ptr[i + 1])

    cpdef SIZE_t first_child(self, SIZE_t i) nogil:
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
        if self._b_ptr[i]:
            if self.is_tip(i):
                return 0
            else:
                return i + 1
        else:
            return self.first_child(self.open(i))

    cpdef SIZE_t last_child(self, SIZE_t i) nogil:
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
        if self._b_ptr[i]:
            if self.is_tip(i):
                return 0
            else:
                return self.open(self.close(i) - 1)
        else:
            return self.last_child(self.open(i))

    def mincount(self, SIZE_t i, SIZE_t j):
        """number of occurrences of the minimum in excess(i), excess(i + 1), . . . , excess(j)."""
        excess, counts = np.unique([self.excess(k) for k in range(i, j + 1)], return_counts=True)
        return counts[excess.argmin()]

    def minselect(self, SIZE_t i, SIZE_t j, SIZE_t q):
        """position of the qth minimum in excess(i), excess(i + 1), . . . , excess(j)."""
        counts = np.array([self.excess(k) for k in range(i, j + 1)])
        index = counts == counts.min()

        if index.sum() < q:
            return None
        else:
            return i + index.nonzero()[0][q - 1]

    cpdef SIZE_t next_sibling(self, SIZE_t i) nogil:
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
        cdef SIZE_t pos

        if self._b_ptr[i]:
            pos = self.close(i) + 1
        else:
            pos = self.next_sibling(self.open(i))

        if pos >= self.size:
            return 0
        elif self._b_ptr[pos]:
            return pos
        else:
            return 0

    cpdef SIZE_t previous_sibling(self, SIZE_t i) nogil:
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
        cdef SIZE_t pos

        if self._b_ptr[i]:
            if self._b_ptr[max(0, i - 1)]:
                return 0

            pos = self.open(i - 1)
        else:
            pos = self.previous_sibling(self.open(i))

        if pos < 0:
            return 0
        elif self._b_ptr[pos]:
            return pos
        else:
            return 0

    cpdef SIZE_t preorder_rank(self, SIZE_t i) nogil:
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
        if self._b_ptr[i]:
            return self.rank(1, i)
        else:
            return self.preorder_rank(self.open(i))

    cpdef SIZE_t preorder_select(self, SIZE_t k) nogil:
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
        return self.select(1, k)

    cpdef SIZE_t postorder_rank(self, SIZE_t i) nogil:
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
        if self._b_ptr[i]:
            return self.rank(0, self.close(i))
        else:
            return self.rank(0, i)

    cpdef SIZE_t postorder_select(self, SIZE_t k) nogil:
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
        return self.open(self.select(0, k))

    cpdef BOOL_t is_ancestor(self, SIZE_t i, SIZE_t j) nogil:
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
        if i == j:
            return False

        if not self._b_ptr[i]:
            i = self.open(i)

        return i <= j < self.close(i)

    cpdef SIZE_t count(self, SIZE_t i=0, bint tips=False) nogil:
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
        cdef:
            SIZE_t last, j, c

        if not self._b_ptr[i]:
            i = self.open(i)

        if not tips:
            return (self.close(i) - i + 1) / 2

        # a tip is an open parenthesis immediately followed by a close; count
        # them within the span [i, close(i)] of the subtree
        last = self.close(i)
        j = i
        c = 0
        while j < last:
            if self._b_ptr[j] and not self._b_ptr[j + 1]:
                c += 1
                j += 1
            j += 1

        return c

    cpdef SIZE_t level_ancestor(self, SIZE_t i, SIZE_t d) nogil:
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
        if d <= 0:
            return -1

        if not self._b_ptr[i]:
            i = self.open(i)

        return self.bwdsearch(i, -d - 1) + 1

    cpdef SIZE_t level_next(self, SIZE_t i) nogil:
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
        return self.fwdsearch(self.close(i), 1)

    cpdef SIZE_t lca(self, SIZE_t i, SIZE_t j) nogil:
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
           The index of the lowest common ancestor
        """
        if self.is_ancestor(i, j):
            return i
        elif self.is_ancestor(j, i):
            return j
        else:
            return self.parent(self.rmq(i, j) + 1)

    cpdef SIZE_t deepest_node(self, SIZE_t i) nogil:
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
        return self.rMq(self.open(i), self.close(i))

    cpdef SIZE_t height(self, SIZE_t i) nogil:
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
        return self.excess(self.deepest_node(i)) - self.excess(self.open(i))

    cpdef BPTree shear(self, set tips):
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
        cdef:
            SIZE_t i, n = len(tips)
            SIZE_t p, t, count = 0
            BIT_ARRAY* mask
            BPTree new_bp

        mask = bit_array_create(self.data.size)
        bit_array_set_bit(mask, self.root())
        bit_array_set_bit(mask, self.close(self.root()))

        for i in range(self.data.size):
            # is_tip is only defined on the open parenthesis
            if self.is_tip(i):
                if self.name(i) in tips:  # gil is required for set operation
                    with nogil:
                        count += 1
                        bit_array_set_bit(mask, i)
                        bit_array_set_bit(mask, i + 1)

                        p = self.parent(i)
                        while p != 0 and bit_array_get_bit(mask, p) == 0:
                            bit_array_set_bit(mask, p)
                            bit_array_set_bit(mask, self.close(p))

                            p = self.parent(p)

        if count == 0:
            bit_array_free(mask)
            raise ValueError("No requested tips found")

        new_bp = self._mask_from_self(mask, self._lengths)
        bit_array_free(mask)
        return new_bp

    cdef BPTree _mask_from_self(self, BIT_ARRAY* mask,
                            cnp.ndarray[DOUBLE_t, ndim=1] lengths):
        cdef:
            SIZE_t i, k, n, mask_sum
            cnp.ndarray[BOOL_t, ndim=1] new_b
            cnp.ndarray[object, ndim=1] new_names
            cnp.ndarray[object, ndim=1] names = self._names
            cnp.ndarray[DOUBLE_t, ndim=1] new_lengths
            BOOL_t* new_b_ptr
            DOUBLE_t* lengths_ptr
            DOUBLE_t* new_lengths_ptr

        n = bit_array_length(mask)
        mask_sum = bit_array_num_bits_set(mask)

        k = 0

        lengths_ptr = &lengths[0]

        new_b = np.empty(mask_sum, dtype=BOOL)
        new_names = np.empty(mask_sum, dtype=object)
        new_lengths = np.empty(mask_sum, dtype=DOUBLE)

        new_b_ptr = &new_b[0]
        new_lengths_ptr = &new_lengths[0]

        for i in range(n):
            if bit_array_get_bit(mask, i):
                new_b_ptr[k] = self._b_ptr[i]

                # since names is dtype=object, gil is required
                new_names[k] = names[i]
                new_lengths_ptr[k] = lengths_ptr[i]
                k += 1

        return BPTree(np.asarray(new_b), names=new_names, lengths=new_lengths)

    cpdef BPTree collapse(self):
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
        cdef:
            SIZE_t i, n = self.data.sum()
            SIZE_t current, first, last
            cnp.ndarray[DOUBLE_t, ndim=1] new_lengths
            BIT_ARRAY* mask
            DOUBLE_t* new_lengths_ptr
            BPTree new_bp

        mask = bit_array_create(self.data.size)
        bit_array_set_bit(mask, self.root())
        bit_array_set_bit(mask, self.close(self.root()))

        new_lengths = self._lengths.copy()
        new_lengths_ptr = <DOUBLE_t*>new_lengths.data

        with nogil:
            for i in range(n):
                current = self.preorder_select(i)

                if self.is_tip(current):
                    bit_array_set_bit(mask, current)
                    bit_array_set_bit(mask, self.close(current))
                else:
                    first = self.first_child(current)
                    last = self.last_child(current)

                    if first == last:
                        new_lengths_ptr[first] = new_lengths_ptr[first] + \
                                new_lengths_ptr[current]
                    else:
                        bit_array_set_bit(mask, current)
                        bit_array_set_bit(mask, self.close(current))

        new_bp = self._mask_from_self(mask, new_lengths)
        bit_array_free(mask)
        return new_bp

    cdef int scan_block_forward(self, int i, int k, int b, int d) nogil:
        """Scan a block forward from i.

        Parameters
        ----------
        i : int
            The index position to start from in the tree
        k : int
            The block to explore
        b : int
            The block size
        d : int
            The depth to search for

        Returns
        -------
        int
            The index position of the result. -1 is returned if a result is not
            found.
        """
        cdef int lower_bound
        cdef int upper_bound
        cdef int j

        # lower_bound is block boundary or right of i
        lower_bound = max(k, 0) * b
        lower_bound = max(i + 1, lower_bound)

        # upper_bound is block boundary or end of tree
        upper_bound = min((k + 1) * b, self.size)

        for j in range(lower_bound, upper_bound):
            if self._e_index[j] == d:
                return j

        return -1

    cdef int scan_block_backward(self, int i, int k, int b, int d) nogil:
        """Scan a block backward from i.

        Parameters
        ----------
        i : int
            The index position to start from in the tree
        k : int
            The block to explore
        b : int
            The block size
        d : int
            The depth to search for

        Returns
        -------
        int
            The index position of the result. -1 is returned if a result is not
            found.
        """
        cdef int lower_bound
        cdef int upper_bound
        cdef int j

        # range stop is exclusive, so need to set "stop" at -1 of boundary
        lower_bound = max(k, 0) * b - 1

        # include the right most position of the k-1 block so we can identify
        # closures spanning blocks.
        if lower_bound >= 0:
            lower_bound -= 1

        # upper bound is block boundary or left of i, whichever is less
        # note that this is an inclusive boundary since this is a backward search
        upper_bound = min((k + 1) * b, self.size) - 1
        upper_bound = min(i - 1, upper_bound)

        if upper_bound <= 0:
            return -1

        for j in range(upper_bound, lower_bound, -1):
            if self.excess(j) == d:
                return j

        return -1

    cdef SIZE_t fwdsearch(self, SIZE_t i, int d) nogil:
        """Search forward from i for desired excess.

        Parameters
        ----------
        i : int
            The index to search forward from
        d : int
            The excess difference to search for (relative to E[i])

        Returns
        -------
        int
            The index of the result, or -1 if no result was found
        """
        cdef int k  # the block being interrogated
        cdef int result = -1 # the result of a scan within a block
        cdef int node  # the node within the binary tree being examined

        # get the block of parentheses to check
        k = i // self._rmm.b

        # desired excess
        d += self._e_index[i]

        # determine which node our block corresponds too
        node = bt_node_from_left(k, self._rmm.height)

        # see if our result is in our current block
        if self._rmm.mM[node, self._rmm.m_idx] <= d <= self._rmm.mM[node, self._rmm.M_idx]:
            result = self.scan_block_forward(i, k, self._rmm.b, d)

        # if we do not have a result, we need to begin traversal of the tree
        if result == -1:
            # walk up the tree
            while not bt_is_root(node):
                if bt_is_left_child(node):
                    node = bt_right_sibling(node)
                    if self._rmm.mM[node, self._rmm.m_idx] <= d  <= self._rmm.mM[node, self._rmm.M_idx]:
                        break
                node = bt_parent(node)

            if bt_is_root(node):
                return -1

            # descend until we hit a leaf node
            while not bt_is_leaf(node, self._rmm.height):
                node = bt_left_child(node)

                # evaluate right, if not found, pick left
                if not (self._rmm.mM[node, self._rmm.m_idx] <= d <= self._rmm.mM[node, self._rmm.M_idx]):
                    node = bt_right_sibling(node)

            # we have found a block with contains our solution. convert from the
            # node index back into the block index
            k = node - <int>(pow(2, self._rmm.height) - 1)

            # scan for a result using the original d
            result = self.scan_block_forward(i, k, self._rmm.b, d)

        return result

    cdef SIZE_t bwdsearch(self, SIZE_t i, int d) nogil:
        """Search backward from i for desired excess

        Parameters
        ----------
        i : int
            The index to search forward from
        d : int
            The excess difference to search for (relative to E[i])

        Returns
        -------
        int
            The index of the result, or -1 if no result was found
        """
        cdef int k  # the block being interrogated
        cdef int result = -1 # the result of a scan within a block
        cdef int node  # the node within the binary tree being examined

        # get the block of parentheses to check
        k = i // self._rmm.b

        # desired excess
        d += self.excess(i)

        # see if our result is in our current block
        result = self.scan_block_backward(i, k, self._rmm.b, d)

        # determine which node our block corresponds too
        node = bt_node_from_left(k, self._rmm.height)

        # special case: check sibling
        if result == -1 and bt_is_right_child(node):
            node = bt_left_sibling(node)
            k = node - <int>(pow(2, self._rmm.height) - 1)
            result = self.scan_block_backward(i, k, self._rmm.b, d)

           # reset node and k in the event that result == -1
            k = i // self._rmm.b
            node = bt_right_sibling(node)

        # if we do not have a result, we need to begin traversal of the tree
        if result == -1:
            while not bt_is_root(node):
                # right nodes cannot contain the solution as we are searching left
                # As such, if we are the right node already, evaluate its sibling.
                if bt_is_right_child(node):
                    node = bt_left_sibling(node)
                    if self._rmm.mM[node, self._rmm.m_idx] <= d <= self._rmm.mM[node, self._rmm.M_idx]:
                        break

                # if we did not find a valid node, adjust for the relative
                # excess of the current node, and ascend to the parent
                node = bt_parent(node)

            if bt_is_root(node):
                return -1

            # descend until we hit a leaf node
            while not bt_is_leaf(node, self._rmm.height):
                node = bt_right_child(node)

                # evaluate right, if not found, pick left
                if not (self._rmm.mM[node, self._rmm.m_idx] <= d <= self._rmm.mM[node, self._rmm.M_idx]):
                    node = bt_left_sibling(node)

            # we have found a block with contains our solution. convert from the
            # node index back into the block index
            k = node - <int>(pow(2, self._rmm.height) - 1)

            # scan for a result
            result = self.scan_block_backward(i, k, self._rmm.b, d)

        return result
