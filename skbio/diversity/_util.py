# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity._phylogenetic import _nodes_by_counts


def _validate_counts(counts, cast_int=False):
    """Validate and convert input to an acceptable counts vector type.

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


def _validate_taxa_and_tree(counts, taxa, tree, rooted=True):
    """Validate taxa and tree prior to calculating phylogenetic diversity metrics."""
    len_taxa = len(taxa)
    set_taxa = set(taxa)
    if len_taxa != len(set_taxa):
        raise ValueError("``taxa`` cannot contain duplicated ids.")

    if len(counts) != len_taxa:
        raise ValueError("``taxa`` must be the same length as ``counts`` vector(s).")

    if len(tree.root().children) == 0:
        raise ValueError("``tree`` must contain more than just a root node.")

    if rooted is True and len(tree.root().children) > 2:
        # this is an imperfect check for whether the tree is rooted or not.
        # can this be improved?
        raise ValueError("``tree`` must be rooted.")

    # all nodes (except the root node) have corresponding branch lengths
    # all tip names in tree are unique
    # all taxa correspond to tip names in tree
    branch_lengths = []
    tip_names = []
    for e in tree.traverse():
        if not e.is_root():
            branch_lengths.append(e.length)
        if e.is_tip():
            tip_names.append(e.name)
    set_tip_names = set(tip_names)
    if len(tip_names) != len(set_tip_names):
        raise DuplicateNodeError("All tip names in tree must be unique.")

    if np.array([branch is None for branch in branch_lengths]).any():
        raise ValueError("All non-root nodes in tree must have a branch length.")

    if n_missing := len(set_taxa - set_tip_names):
        raise MissingNodeError(
            f"{n_missing} taxa are not present as tip names in tree."
        )


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


def _validate_table(counts, ids, kwargs):
    """Disallow overriding of sample and feature IDs.

    WARNING: this implicitly adds an entry to kwargs IF `tree` is present.

    """
    if ids is not None:
        raise ValueError("Cannot provide a `Table` as `counts` and `ids`.")

    if "taxa" in kwargs:
        raise ValueError("Cannot provide a `Table` as `counts` and `taxa`.")

    from skbio.table._base import _table_to_numpy

    dense_counts, sample_ids, feature_ids = _table_to_numpy(counts)

    if "tree" in kwargs:
        kwargs["taxa"] = feature_ids

    return dense_counts, sample_ids
