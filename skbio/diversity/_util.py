# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import scipy as sp

from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity._phylogenetic import _nodes_by_counts


def _validate_counts(counts, cast_int=False):
    """Validate and convert input to an acceptable counts vector/matrix type.

    Parameters
    ----------
    counts : array_like of shape (n_taxa,) or (n_samples, n_taxa)
        Vector or matrix of counts.
    cast_int : bool, optional
        Whether cast values into integers, if not already. Default is False.

    Returns
    -------
    ndarray of shape (n_taxa,) or (n_samples, n_taxa)
        Valid counts vector or matrix.

    Raises
    ------
    ValueError
        If input array has an invalid data type.
    ValueError
        If input array contains negative values.

    Notes
    -----
    This function will return the original ``counts`` if it is already a valid counts
    vector. Otherwise it will return an edited copy that is valid.

    The data type of counts must be any subtype of ``np.integer`` (integers) or
    ``np.floating`` (floating-point numbers; excluding complex numbers) [1]_.

    References
    ----------
    .. [1] https://numpy.org/doc/stable/reference/arrays.scalars.html

    """
    counts = np.asarray(counts)

    dtype = counts.dtype
    if np.issubdtype(dtype, np.floating):
        if cast_int:
            counts = counts.astype(int)
    elif not np.issubdtype(dtype, np.integer) and dtype is not np.dtype("bool"):
        raise ValueError("Counts must be integers or floating-point numbers.")

    # This is more efficient that `(counts < 0).any()`, but a more efficient way is to
    # iterate by element and exit early on the first negative value.
    if counts.size > 0 and counts.min() < 0:
        raise ValueError("Counts cannot contain negative values.")

    return counts


def _validate_counts_vector(counts, cast_int=False):
    """Validate and convert input to an acceptable counts vector type.

    Parameters
    ----------
    counts : array_like of int or float of shape (n_taxa,)
        Vector of counts.
    cast_int : bool, optional
        Whether cast values into integers, if not already. Default is False.

    Returns
    -------
    ndarray of int or float of shape (n_taxa,)
        Valid counts vector.

    Raises
    ------
    ValueError
        If counts has more than 1 dimension.

    """
    counts = _validate_counts(counts, cast_int=cast_int)
    if counts.ndim != 1:
        raise ValueError("`counts` must be a 1-D array (vector).")
    return counts


def _validate_counts_matrix(counts, cast_int=False):
    """Validate and convert input to an acceptable counts matrix type.

    Parameters
    ----------
    counts : array_like of shape (n_samples, n_taxa)
        Matrix of counts.
    cast_int : bool, optional
        Whether cast values into integers, if not already. Default is False.

    Returns
    -------
    ndarray of shape (n_samples, n_taxa)
        Valid counts matrix.

    Raises
    ------
    ValueError
        If counts has more than 2 dimensions.

    """
    msg = "`counts` has {} dimensions whereas up to 2 dimensions are allowed."
    counts = _validate_counts(counts, cast_int=cast_int)
    counts = np.atleast_2d(counts)
    if counts.ndim > 2:
        raise ValueError(msg.format(counts.ndim))
    return counts


def vectorize_counts_and_tree(counts, taxa, tree):
    """Index tree and convert counts to np.array in corresponding order.

    Parameters
    ----------
    counts : array_like of shape (n_samples, n_taxa) or (n_taxa,)
        Counts/abundances of taxa in one or multiple samples.
    taxa : array_like of shape (n_taxa,)
        Taxon IDs corresponding to tip names in `tree`.
    tree : skbio.TreeNode
        Tree relating taxa. The set of tip names in the tree can be a superset
        of `taxa`, but not a subset.

    Returns
    -------
    ndarray of shape (n_samples, n_nodes)
        Total counts/abundances of taxa descending from individual nodes of the tree.
    dict of array
        Indexed tree. See `to_array`.
    ndarray of shape (n_nodes,)
        Branch lengths of corresponding nodes of the tree.

    See Also
    --------
    skbio.tree.TreeNode.to_array

    Notes
    -----
    Leveraging internal node counts in the tree (in addition to tip abundances) can
    double the accuracy in downstream machine learning pipelines [1]_.

    References
    ----------
    .. [1] Martino C., McDonald D., Cantrell K., Dilmore AH., Vázquez-Baeza Y.,
       Shenhav L., Shaffer J.P., Rahman G., Armstrong G., Allaband C., Song S.J.,
       Knight R. Compositionally aware phylogenetic beta-diversity measures better
       resolve microbiomes associated with phenotype. mSystems. 7(3) (2022).

    """
    tree_index = tree.to_array(nan_length_value=0.0)
    taxa = np.asarray(taxa)
    counts = np.atleast_2d(counts)
    counts_by_node = _nodes_by_counts(counts, taxa, tree_index)
    branch_lengths = tree_index["length"]

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths


def vectorize_counts_and_tree_sparse(count, taxa, tree):
    # convert tree into a dict
    tree_index = tree.to_array(nan_length_value=0.0)

    # convert the taxa into an array
    taxa = np.asarray(taxa)

    # get the branch lengths from the tree dict. these are an array
    branch_lengths = tree_index["length"]

    # count has to be a sparse array, so convert it to a sparse array
    # count can be an array and this would convert it to a sparse array... again
    # if count already is a sparse array, then this has no effect (redundant basically)
    count = sp.csr_array(count)

    # best to use a helper method to do the rest of the calculations
    # think of this as a 3 step process
    # 1.) get everything as the apropriate types
    # 2.) create empty sparse array
    # 3.) get counts from tips to root and do inplace manipulations
    # so we are on step 1 going to step 2
    counts_by_node = _nodes_by_counts_sparse(count, taxa, tree_index)

    # the returned counts_by_node is a placeholder right now
    # what exactly do we want to return? data, indices, and/or indptr?
    # maybe just return the whole thing since those 3 are attributes that
    # can be accessed with .data and such
    # but it's called VECTORIZE so we should return a vector...
    # need clarification here
    return counts_by_node.T, tree_index, branch_lengths


def _nodes_by_counts_sparse(counts, tip_ids, indexed):
    # get the names of the nodes of the tree. includes root
    # taxa, features, nodes, etc. go by many names
    # structure: tips, parent, grandparent, root
    nodes = indexed["name"]

    # might not need to reverse it if we expect results to be
    # descendants, parents, grandparents, root
    # vibe code seems like it wants root, grandparents, parents, descendants
    nodes = nodes.reverse()

    # tuple (rows, cols)
    # gives the shape of the matrix
    # may not need
    counts_dimensions = counts.shape
    num_of_samples = counts_dimensions[0]
    num_of_features = counts_dimensions[1]

    # scipy has a similar sum() and nonzero() functions
    # used to get the indices of where the tips are
    observed_indices = counts.sum(axis=0).nonzero()[0]

    # get a set of the tip names. these are the descendants
    tree_tips = tip_ids[observed_indices]
    tree_tips_set = set(tree_tips)

    # generate list of sample ids. unfortunately, they are just numerical since
    # we aren't given the sample ids (only taxa)
    sample_ids = []
    for i in range(num_of_samples):
        sample_ids.append(i)

    # construct mappings of the observed to their positions in the node array
    # so make the arrays for each: data, row, col
    node_lookup = {}
    rows, cols, data = [], [], []

    # for vibe code i am assuming the outer for loop is going through all nodes of tree
    # inner for loop is for going through the samples in the specified node
    # we likely wouldnt need the loops after this loop (as seen in og code)
    # since here we are making the arrays needed for the sparse matrix we will make
    # once this loop completes
    # so then we can do the matrices caluclations and return the resulting
    # sparse array

    # MODIFY THIS LOOP!!!! IT IS THE MEAT/HEART OF THE WHOLE FUNCTION!!!
    # Not completed and does not work properly
    # we would want something like root, root->c, root->c->a, root->c->b
    # we would use the path() method from TreeNode to get the path and properly
    # label rows and cols on the sparse matrix
    for i, feature in enumerate(nodes):
        for j in tree_tips_set:
            idx = node_lookup.setdefault(j, len(node_lookup))
            rows.append(i)
            cols.append(idx)
            data.append(1.0)

    # apparently matrix is a legacy function and array is recommended to use
    # Source: https://numpy.org/doc/2.3/reference/generated/numpy.matrix.html
    new_csr_array = sp.csr_array(
        data, (rows, cols), shape=(len(nodes), len(node_lookup))
    )

    final_csr_array = counts @ new_csr_array

    # Note: the rows of original and the resulting arrays are always the same for the
    # original function. Columns change since input does not include
    # parents and root and the output does include them
    return final_csr_array


def _get_phylogenetic_kwargs(kwargs, taxa):
    if taxa is None:
        raise ValueError("`taxa` is required for phylogenetic diversity metrics.")
    try:
        tree = kwargs.pop("tree")
    except KeyError:
        raise ValueError("`tree` is required for phylogenetic diversity metrics.")

    return taxa, tree, kwargs


def _qualify_counts(counts):
    return counts > 0.0
