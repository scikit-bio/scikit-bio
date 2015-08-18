# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from skbio.tree import DuplicateNodeError, MissingNodeError


def _validate_counts_vector(counts, suppress_cast=False):
    """Validate and convert input to an acceptable counts vector type.

    Note: may not always return a copy of `counts`!

    """
    counts = np.asarray(counts)

    if not suppress_cast:
        counts = counts.astype(int, casting='safe', copy=False)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts


def _validate_counts_vectors(*args, **kwargs):
    results = []
    lens = []
    # py2-compatible mechanism for specifying a keyword argument when also
    # passing *args derived from SO answer:
    # http://stackoverflow.com/a/15302038/3424666
    suppress_cast = kwargs.pop('suppress_cast', False)
    for counts in args:
        results.append(_validate_counts_vector(counts, suppress_cast))
        lens.append(len(counts))
    if len(set(lens)) > 1:
        raise ValueError("Input vectors u_counts and v_counts must be of "
                         "equal length.")

    return results


def _validate_otu_ids_and_tree(counts, otu_ids, tree):
    # all otu_ids are unique
    # len(otu_ids) == len(counts)
    len_otu_ids = len(otu_ids)
    set_otu_ids = set(otu_ids)
    if len_otu_ids != len(set_otu_ids):
        raise ValueError("OTU IDs vector cannot contain duplicated ids.")
    if len(counts) != len_otu_ids:
        raise ValueError("OTU IDs vector must be the same length as counts "
                         "vector(s).")

    # the tree is rooted
    if len(tree.root().children) > 2:
        # this is an imperfect check for whether the tree is rooted or not.
        # can this be improved?
        raise ValueError("Tree must be rooted.")

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
    if np.array([l is None for l in branch_lengths]).any():
        raise ValueError("All non-root nodes in tree must have a branch "
                         "length.")
    missing_tip_names = set_otu_ids - set_tip_names
    if missing_tip_names != set():
        raise MissingNodeError("All otu_ids must be present as tip names in "
                               "tree. Tree is missing tips with names: %s"
                               % " ".join(missing_tip_names))
