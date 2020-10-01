# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections

import numpy as np
import pandas as pd

from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity._phylogenetic import _nodes_by_counts


def _validate_counts_vector(counts, suppress_cast=False):
    """Validate and convert input to an acceptable counts vector type.

    Note: may not always return a copy of `counts`!

    """
    counts = np.asarray(counts)
    if not np.all(np.isreal(counts)):
        raise ValueError("Counts vector must contain real-valued entries.")
    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts


def _validate_counts_matrix(counts, ids=None, suppress_cast=False):
    results = []

    # handle case of where counts is a single vector by making it a matrix.
    # this has to be done before forcing counts into an ndarray because we
    # don't yet know that all of the entries are of equal length
    if isinstance(counts, pd.core.frame.DataFrame):
        if ids is not None and len(counts.index) != len(ids):
            raise ValueError(
                "Number of rows in ``counts``"
                " must be equal to number of provided ``ids``."
            )
        return np.asarray(counts)
    else:

        if len(counts) == 0 or not isinstance(counts[0], collections.Iterable):
            counts = [counts]
        counts = np.asarray(counts)
        if counts.ndim > 2:
            raise ValueError(
                "Only 1-D and 2-D array-like objects can be provided "
                "as input. Provided object has %d dimensions." %
                counts.ndim)

        if ids is not None and len(counts) != len(ids):
            raise ValueError(
                "Number of rows in ``counts`` must be equal "
                "to number of provided ``ids``."
            )

        lens = []
        for v in counts:
            results.append(_validate_counts_vector(v, suppress_cast))
            lens.append(len(v))
        if len(set(lens)) > 1:
            raise ValueError(
                "All rows in ``counts`` must be of equal length."
            )
        return np.asarray(results)


def _validate_otu_ids_and_tree(counts, otu_ids, tree):
    len_otu_ids = len(otu_ids)
    set_otu_ids = set(otu_ids)
    if len_otu_ids != len(set_otu_ids):
        raise ValueError("``otu_ids`` cannot contain duplicated ids.")

    if len(counts) != len_otu_ids:
        raise ValueError("``otu_ids`` must be the same length as ``counts`` "
                         "vector(s).")

    if len(tree.root().children) == 0:
        raise ValueError("``tree`` must contain more than just a root node.")

    if len(tree.root().children) > 2:
        # this is an imperfect check for whether the tree is rooted or not.
        # can this be improved?
        raise ValueError("``tree`` must be rooted.")

    # all nodes (except the root node) have corresponding branch lengths
    # all tip names in tree are unique
    # all otu_ids correspond to tip names in tree
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
        raise ValueError("All non-root nodes in ``tree`` must have a branch "
                         "length.")
    missing_tip_names = set_otu_ids - set_tip_names
    if missing_tip_names != set():
        n_missing_tip_names = len(missing_tip_names)
        raise MissingNodeError("All ``otu_ids`` must be present as tip names "
                               "in ``tree``. ``otu_ids`` not corresponding to "
                               "tip names (n=%d): %s" %
                               (n_missing_tip_names,
                                " ".join(missing_tip_names)))


def _vectorize_counts_and_tree(counts, otu_ids, tree):
    """ Index tree and convert counts to np.array in corresponding order
    """
    tree_index = tree.to_array(nan_length_value=0.0)
    otu_ids = np.asarray(otu_ids)
    counts = np.atleast_2d(counts)
    counts_by_node = _nodes_by_counts(counts, otu_ids, tree_index)
    branch_lengths = tree_index['length']

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths


def _get_phylogenetic_kwargs(counts, **kwargs):
    try:
        otu_ids = kwargs.pop('otu_ids')
    except KeyError:
        raise ValueError("``otu_ids`` is required for phylogenetic diversity "
                         "metrics.")
    try:
        tree = kwargs.pop('tree')
    except KeyError:
        raise ValueError("``tree`` is required for phylogenetic diversity "
                         "metrics.")

    return otu_ids, tree, kwargs
