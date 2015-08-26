import numpy as np
cimport numpy as cnp
cimport cython


@cython.boundscheck(False) 
@cython.wraparound(False) 
cdef traverse_reduce(tuple child_index, cnp.ndarray[cnp.int64_t, ndim=2] a):
    """Apply a[k] = sum[i:j]

    TODO: describe in detail, and include a hand example of how the reduction
    works. It's awesome.
    """
    cdef:
        cnp.int64_t node, start, end
        Py_ssize_t i, j
        tuple child_data

    for i in range(len(child_index)):
        child_data = child_index[i]
        node = child_data[0]
        start = child_data[1]
        end = child_data[2]

        for j in range(start, end + 1):
            a[node] += a[j]


@cython.boundscheck(False) 
@cython.wraparound(False) 
def nodes_by_counts(cnp.ndarray[cnp.int64_t, ndim=2] counts, 
                    cnp.ndarray tip_ids, 
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
        cnp.ndarray nodes, observed_ids
        cnp.ndarray[cnp.int64_t, ndim=2] count_array, counts_t
        cnp.ndarray[cnp.int64_t, ndim=1] observed_indices, otus_in_nodes
        Py_ssize_t i, j
        set observed_ids_set
        object n
        dict node_lookup
        cnp.int64_t n_count_vectors, n_count_otus

    nodes = indexed['name']

    # allow counts to be a vector
    counts = np.atleast_2d(counts)

    # determine observed IDs. It may be worth breaking apart the reduce call
    # and optimizing it futher.
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
    otus_in_nodes = np.zeros(observed_ids.shape[0], dtype=np.int64)
    for i in range(observed_ids.shape[0]):
        n = observed_ids[i]
        otus_in_nodes[i] = node_lookup[n]

    # count_array has a row per node (not tip) and a column per env.
    n_count_vectors = counts.shape[0]
    count_array = np.zeros((nodes.shape[0], n_count_vectors), dtype=np.int64)
   
    # populate the counts array with the counts of each observation in each
    # env
    counts_t = counts.transpose()
    n_count_otus = otus_in_nodes.shape[0]
    for i in range(n_count_otus):
        for j in range(n_count_vectors):
            count_array[otus_in_nodes[i], j] = counts_t[observed_indices[i], j]
    
    traverse_reduce(indexed['child_index'], count_array)

    return count_array
