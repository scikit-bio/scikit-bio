# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int64
ctypedef np.int64_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def _tip_distances(np.ndarray[np.double_t, ndim=1] a, object t,
                   np.ndarray[DTYPE_t, ndim=1] tip_indices):
    """Sets each tip to its distance from the root

    Parameters
    ----------
    a : np.ndarray of double
        A matrix in which each row corresponds to a node in ``t``.
    t : skbio.tree.TreeNode
        The tree that corresponds to the rows in ``a``.
    tip_indices : np.ndarray of int
        The index positions in ``a`` of the tips in ``t``.

    Returns
    -------
    np.ndarray of double
        A matrix in which each row corresponds to a node in ``t``, Only the
        rows that correspond to tips are nonzero, and the values in these rows
        are the distance from that tip to the root of the tree.
    """
    cdef:
        object n
        Py_ssize_t i, p_i, n_rows
        np.ndarray[np.double_t, ndim=1] mask
        np.ndarray[np.double_t, ndim=1] tip_ds = a.copy()

    # preorder reduction over the tree to gather distances at the tips
    n_rows = tip_ds.shape[0]
    for n in t.preorder(include_self=False):
        i = n.id
        p_i = n.parent.id

        tip_ds[i] += tip_ds[p_i]

    # construct a mask that represents the locations of the tips
    mask = np.zeros(n_rows, dtype=np.double)
    for i in range(tip_indices.shape[0]):
        mask[tip_indices[i]] = 1.0

    # apply the mask such that tip_ds only includes values which correspond to
    # the tips of the tree.
    for i in range(n_rows):
        tip_ds[i] *= mask[i]

    return tip_ds


@cython.boundscheck(False)
@cython.wraparound(False)
cdef _traverse_reduce(np.ndarray[DTYPE_t, ndim=2] child_index,
                      np.ndarray[DTYPE_t, ndim=2] a):
    """Apply a[k] = sum[i:j]

    Parameters
    ----------
    child_index: np.array of int
        A matrix in which the first column corresponds to an index position in
        ``a``, which represents a node in a tree. The second column is the
        starting index in ``a`` for the node's children, and the third column
        is the ending index in ``a`` for the node's children.
    a : np.ndarray of int
        A matrix of the environment data. Each row corresponds to a node in a
        tree, and each column corresponds to an environment. On input, it is
        assumed that only tips have counts.

    Notes
    -----
    This is effectively a postorder reduction over the tree. For example,
    given the following tree:

                            /-A
                  /E-------|
                 |          \-B
        -root----|
                 |          /-C
                  \F-------|
                            \-D

    And assuming counts for [A, B, C, D] in environment FOO of [1, 1, 1, 0] and
    counts for environment BAR of [0, 1, 1, 1], the input counts matrix ``a``
    would be:

        [1 0  -> A
         1 1  -> B
         1 1  -> C
         0 1  -> D
         0 0  -> E
         0 0  -> F
         0 0] -> root

    The method will perform the following reduction:

        [1 0     [1 0     [1 0     [1 0
         1 1      1 1      1 1      1 1
         1 1      1 1      1 1      1 1
         0 1  ->  0 1  ->  0 1  ->  0 1
         0 0      2 1      2 1      2 1
         0 0      0 0      1 2      1 2
         0 0]     0 0]     0 0]     3 3]

    The index positions of the above are encoded in ``child_index`` which
    describes the node to aggregate into, and the start and stop index
    positions of the nodes immediate descendents.

    This method operates inplace on ``a``
    """
    cdef:
        Py_ssize_t i, j, k
        DTYPE_t node, start, end
        DTYPE_t n_envs = a.shape[1]

    # possible GPGPU target
    for i in range(child_index.shape[0]):
        node = child_index[i, 0]
        start = child_index[i, 1]
        end = child_index[i, 2]

        for j in range(start, end + 1):
            for k in range(n_envs):
                a[node, k] += a[j, k]


@cython.boundscheck(False)
@cython.wraparound(False)
def _nodes_by_counts(np.ndarray counts,
                     np.ndarray tip_ids,
                     dict indexed):
    """Construct the count array, and the counts up the tree

    Parameters
    ----------
    counts : np.array of int
        A 1D or 2D vector in which each row corresponds to the observed counts
        in an environment. The rows are expected to be in order with respect to
        `tip_ids`.
    tip_ids : np.array of str
        A vector of tip names that correspond to the columns in the `counts`
        matrix.
    indexed : dict
        The result of `index_tree`.

    Returns
    -------
    np.array of int
        The observed counts of every node and the counts if its descendents.

    """
    cdef:
        np.ndarray nodes, observed_ids
        np.ndarray[DTYPE_t, ndim=2] count_array, counts_t
        np.ndarray[DTYPE_t, ndim=1] observed_indices, otus_in_nodes
        Py_ssize_t i, j
        set observed_ids_set
        object n
        dict node_lookup
        DTYPE_t n_count_vectors, n_count_otus

    nodes = indexed['name']

    # allow counts to be a vector
    counts = np.atleast_2d(counts)
    counts = counts.astype(DTYPE, copy=False)

    # determine observed IDs. It may be possible to unroll these calls to
    # squeeze a little more performance
    observed_indices = counts.sum(0).nonzero()[0]
    observed_ids = tip_ids[observed_indices]
    observed_ids_set = set(observed_ids)

    # construct mappings of the observed to their positions in the node array
    node_lookup = {}
    for i in range(nodes.shape[0]):
        n = nodes[i]
        if n in observed_ids_set:
            node_lookup[n] = i

    # determine the positions of the observed IDs in nodes
    otus_in_nodes = np.zeros(observed_ids.shape[0], dtype=DTYPE)
    for i in range(observed_ids.shape[0]):
        n = observed_ids[i]
        otus_in_nodes[i] = node_lookup[n]

    # count_array has a row per node (not tip) and a column per env.
    n_count_vectors = counts.shape[0]
    count_array = np.zeros((nodes.shape[0], n_count_vectors), dtype=DTYPE)

    # populate the counts array with the counts of each observation in each
    # env
    counts_t = counts.transpose()
    n_count_otus = otus_in_nodes.shape[0]
    for i in range(n_count_otus):
        for j in range(n_count_vectors):
            count_array[otus_in_nodes[i], j] = counts_t[observed_indices[i], j]

    child_index = indexed['child_index'].astype(DTYPE, copy=False)
    _traverse_reduce(child_index, count_array)

    return count_array
