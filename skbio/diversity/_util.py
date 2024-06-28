# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections.abc

import numpy as np
import pandas as pd

from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity._phylogenetic import _nodes_by_counts


def _validate_counts_vector(counts, cast_int=False):
    """Validate and convert input to an acceptable counts vector type.

    Parameters
    ----------
    counts : array_like of int or float of shape (n_taxa,)
        Vector of counts.
    cast_int : bool, optional
        Cast values into integers, if not already. ``False`` by default.

    Returns
    -------
    ndarray of int or float of shape (n_taxa,)
        Valid counts vector.

    Raises
    ------
    ValueError
        If input array has an invalid data type.
    ValueError
        If input array is not 1-D.
    ValueError
        If there are negative values.

    Notes
    -----
    This function will return the original ``counts`` if it is already a valid counts
    vector. Otherwise it will return an edited copy that is valid.

    The data type of counts must be any subtype of ``np.integer`` (integers) or
    ``np.floating`` (floating-point numbers; excluding complex numbers) [1]_.

    See Also
    --------
    _validate_counts_matrix

    References
    ----------
    .. [1] https://numpy.org/doc/stable/reference/arrays.scalars.html

    """
    counts = np.asarray(counts)

    # counts must be int or float
    if np.issubdtype(dtype := counts.dtype, np.floating):
        # cast values into integers
        if cast_int:
            counts = counts.astype(int)

    elif not np.issubdtype(dtype, np.integer) and dtype is not np.dtype("bool"):
        raise ValueError("Counts must be integers or floating-point numbers.")

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    if (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts


def _validate_counts_matrix(counts, ids=None, cast_int=False):
    """Validate and convert input to an acceptable counts matrix type.

    Parameters
    ----------
    counts : array_like of shape (n_samples, n_taxa)
        Matrix of counts.
    ids : array_like of shape (n_samples,), optional
        Sample IDs to check against counts dimensions.
    cast_int : bool, optional
        Cast values into integers, if not already. ``False`` by default.

    Returns
    -------
    ndarray of shape (n_samples, n_taxa)
        Valid counts matrix.

    See Also
    --------
    _validate_counts_vector

    """
    lenerr = "Number of rows in `counts` must be equal to number of provided `ids`."

    # handle pandas data frame
    if isinstance(counts, pd.DataFrame):
        if ids is not None and counts.shape[0] != len(ids):
            raise ValueError(lenerr)
        counts = counts.to_numpy()

    else:
        # convert counts into a 2-D array
        # will raise ValueError if row lengths are unequal
        counts = np.atleast_2d(counts)

        if counts.ndim > 2:
            raise ValueError(
                "Only 1-D and 2-D array-like objects can be provided as input. "
                f"Provided object has {counts.ndim} dimensions."
            )

        if ids is not None and counts.shape[0] != len(ids):
            raise ValueError(lenerr)

    # counts must be int or float
    if np.issubdtype(dtype := counts.dtype, np.floating):
        # cast values into integers
        if cast_int:
            counts = counts.astype(int)

    elif not np.issubdtype(dtype, np.integer) and dtype is not np.dtype("bool"):
        raise ValueError("Counts must be integers or floating-point numbers.")

    # negative values are not allowed
    # TODO: `counts < 0` creates a Boolean array of the same shape, which could be
    # memory-inefficient if the input array is very large. Should optimize.
    # See: https://stackoverflow.com/questions/75553212/
    if (counts < 0).any():
        raise ValueError("Counts cannot contain negative values.")

    return counts


def _validate_taxa_and_tree(counts, taxa, tree, rooted=True):
    """Validate taxa and tree prior to calculating phylogenetic diversity metrics."""
    len_taxa = len(taxa)
    set_taxa = set(taxa)
    if len_taxa != len(set_taxa):
        raise ValueError("``taxa`` cannot contain duplicated ids.")

    if len(counts) != len_taxa:
        raise ValueError("``taxa`` must be the same length as ``counts`` " "vector(s).")

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
        raise DuplicateNodeError("All tip names must be unique.")

    if np.array([branch is None for branch in branch_lengths]).any():
        raise ValueError("All non-root nodes in ``tree`` must have a branch " "length.")
    missing_tip_names = set_taxa - set_tip_names
    if missing_tip_names != set():
        n_missing_tip_names = len(missing_tip_names)
        raise MissingNodeError(
            "All ``taxa`` must be present as tip names "
            "in ``tree``. ``taxa`` not corresponding to "
            "tip names (n=%d): %s" % (n_missing_tip_names, " ".join(missing_tip_names))
        )


def _vectorize_counts_and_tree(counts, taxa, tree):
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

    """
    tree_index = tree.to_array(nan_length_value=0.0)
    taxa = np.asarray(taxa)
    counts = np.atleast_2d(counts)
    counts_by_node = _nodes_by_counts(counts, taxa, tree_index)
    branch_lengths = tree_index["length"]

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths


def _get_phylogenetic_kwargs(counts, **kwargs):
    try:
        taxa = kwargs.pop("taxa")
    except KeyError:
        raise ValueError("``taxa`` is required for phylogenetic diversity " "metrics.")
    try:
        tree = kwargs.pop("tree")
    except KeyError:
        raise ValueError("``tree`` is required for phylogenetic diversity " "metrics.")

    return taxa, tree, kwargs


def _quantitative_to_qualitative_counts(counts):
    return counts > 0.0


def _check_taxa_alias(taxa, tree, otu_ids):
    # make `taxa` an alias of `taxa`; for backward compatibility
    if taxa is None:
        if otu_ids is None:
            raise ValueError("A list of taxon IDs must be provided.")
        taxa = otu_ids
    if tree is None:
        raise ValueError("A phylogenetic tree must be provided.")
    return taxa


def _table_to_numpy(table):
    """Convert a skbio.table.Table to a dense representation.

    This is a stop-gap solution to allow current Table objects to interoperate
    with existing driver methods, until they transition to be "sparse" aware.
    """
    sample_ids = list(table.ids())
    obs_ids = list(table.ids(axis="observation"))

    if table.is_empty():
        counts = np.array([[]] * len(sample_ids))
    else:
        counts = table.matrix_data.T.toarray()

    return counts, sample_ids, obs_ids


def _validate_table(counts, ids, kwargs):
    """Disallow overriding of sample and feature IDs.

    WARNING: this implicitly adds an entry to kwargs IF `tree` is present.
    """
    if ids is not None:
        raise ValueError("Cannot provide a `Table` as `counts` and `ids`")

    if "taxa" in kwargs:
        raise ValueError("Cannot provide a `Table` as `counts` and `taxa`")

    dense_counts, sample_ids, feature_ids = _table_to_numpy(counts)
    if "tree" in kwargs:
        kwargs["taxa"] = feature_ids

    return dense_counts, sample_ids
