import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int64
ctypedef np.int64_t DTYPE_t


@cython.boundscheck(False) 
@cython.wraparound(False) 
def tip_distances(np.ndarray[np.double_t, ndim=2] a, object t, 
                  np.ndarray[DTYPE_t, ndim=1] tip_indices):
    """Sets each tip to its distance from the root.
    Note: This will need its own unittest"""
    cdef:
        object n
        Py_ssize_t i, p_i, j, n_rows, n_cols
        np.ndarray[np.double_t, ndim=1] mask

    ### preorder_reduce
    ##### take bind_to_parent_array result, dump down to cythonz
    n_rows = a.shape[0]
    n_cols = a.shape[1]
    for n in t.preorder(include_self=False):
        i = n.id
        p_i = n.parent.id

        for j in range(n_cols):
            a[i, j] += a[p_i, j]

    #for i, s in bound_indices:
    #    i += s
    mask = np.zeros(n_rows, dtype=np.double)
    for i in range(tip_indices.shape[0]):
        mask[tip_indices[i]] = 1.0
    #np.put(mask, tip_indices, 1)

    for i in range(n_rows):
        for j in range(n_cols):
            a[i, j] *= mask[i]
    #a *= mask[:, np.newaxis]

@cython.boundscheck(False) 
@cython.wraparound(False) 
cdef traverse_reduce(np.ndarray[DTYPE_t, ndim=2] child_index, 
                     np.ndarray[DTYPE_t, ndim=2] a):
    """Apply a[k] = sum[i:j]
    ### postorder_reduce
    TODO: describe in detail, and include a hand example of how the reduction
    works. It's awesome.
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
def nodes_by_counts(np.ndarray counts, 
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
    counts = counts.astype(DTYPE)

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
    
    traverse_reduce(indexed['child_index'], count_array)

    return count_array
