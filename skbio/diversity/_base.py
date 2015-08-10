# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np


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


def _validate_counts_vectors(u_counts, v_counts, suppress_cast=False):
    u_counts = _validate_counts_vector(u_counts, suppress_cast)
    v_counts = _validate_counts_vector(v_counts, suppress_cast)

    len_u_counts = len(u_counts)
    len_v_counts = len(v_counts)
    if len_u_counts != len_v_counts:
        raise ValueError("Input vectors u_counts and v_counts must be of "
                         "equal length.")

    return u_counts, v_counts


def _validate_otu_ids_and_tree(counts, otu_ids, tree):
    len_otu_ids = len(otu_ids)
    set_otu_ids = set(otu_ids)
    if len_otu_ids != len(set_otu_ids):
        raise ValueError("OTU IDs vector cannot contain duplicated ids.")
    if len(counts) != len_otu_ids:
        raise ValueError("OTU IDs vector must be the same length as counts "
                         "vector(s).")

    branch_lengths = [e.length for e in tree.traverse() if not e.is_root()]
    if np.array([l is None for l in branch_lengths]).any():
        raise ValueError("All non-root nodes in tree must have a branch "
                         "length.")
