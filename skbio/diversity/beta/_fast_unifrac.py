# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from skbio.util._decorator import experimental
from ._unifrac import _validate
from skbio.diversity._fast_base import _counts_and_length


def unifrac(branch_lengths, i, j):
    """Calculates unifrac(i,j) from branch lengths and cols i and j of m.

    This is the original, unweighted UniFrac metric.

    Parameters
    ----------
    branch_lengths: np.array
         branch_lengths should be row vector, same length as # nodes in tree.
    i : np.array
        a slice of states from m, same length as # nodes in tree.
        Slicing m (e.g. m[:,i]) returns a vector in the right format; note
        that it should be a row vector (the default), not a column vector.
    j : np.array
        a slice of states from m, same length as # nodes in tree.
        Slicing m (e.g. m[:,j]) returns a vector in the right format; note
        that it should be a row vector (the default), not a column vector.

    Returns
    -------
    float
        Unweighted unifrac metric

    Notes
    -----
    This is cogent.maths.unifrac.fast_tree.unifrac, but there are
    other metrics that can (should?) be ported like:
     - unnormalized_unifrac
     - G
     - unnormalized_G
    """

    _or, _and = np.logical_or(i, j), np.logical_and(i, j)
    return 1 - ((branch_lengths * _and).sum() / (branch_lengths * _or).sum())


def w_unifrac(branch_lengths, tip_indices, i, j):
    """Calculates weighted unifrac(i,j) from branch lengths and cols i,j of m.
    """
    i_sum = (np.take(i, tip_indices)).sum()
    j_sum = (np.take(j, tip_indices)).sum()

    if i_sum:
        i_ = i / i_sum
    else:
        i_ = 0.0

    if j_sum:
        j_ = j / j_sum
    else:
        j_ = 0.0

    return (branch_lengths * abs(i_ - j_)).sum()


@experimental(as_of="0.4.0")
def unweighted_unifrac_fast(u_counts, v_counts, otu_ids, tree,
                            validate=True, indexed=None, **kwargs):
    """TODO: add or just pull from educational's __doc__?"""
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    u_counts = np.asarray(u_counts)
    v_counts = np.asarray(v_counts)
    otu_ids = np.asarray(otu_ids)

    boundary = _boundary_case(u_counts, v_counts)
    if boundary is not None:
        return boundary

    counts = np.vstack([u_counts, v_counts])
    count_array, length = _counts_and_length(counts, otu_ids, tree, indexed)

    return unifrac(length, count_array[:, 0], count_array[:, 1])


@experimental(as_of="0.4.0")
def weighted_unifrac_fast(u_counts, v_counts, otu_ids, tree, normalized=False,
                          validate=True, indexed=None, **kwargs):
    """TODO: add or just pull from educational's __doc__?"""
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    u_counts = np.asarray(u_counts)
    v_counts = np.asarray(v_counts)
    otu_ids = np.asarray(otu_ids)

    boundary = _boundary_case(u_counts, v_counts, normalized, unweighted=False)
    if boundary is not None:
        return boundary

    counts = np.vstack([u_counts, v_counts])
    count_array, length = _counts_and_length(counts, otu_ids, tree, indexed)

    tip_indices = [n.id for n in tree.tips()]
    tip_ds = length.copy()[:, np.newaxis]
    bindings = bind_to_parent_array(tree, tip_ds)
    tip_distances(tip_ds, bindings, tip_indices)

    u_counts = count_array[:, 0]
    if count_array.shape[1] == 1:
        v_counts = np.zeros(count_array.shape[0])
    else:
        v_counts = count_array[:, 1]

    u_sum = (np.take(u_counts, tip_indices)).sum()
    v_sum = (np.take(v_counts, tip_indices)).sum()

    u = w_unifrac(length, tip_indices, u_counts, v_counts)
    if normalized:
        u /= _branch_correct(tip_ds, u_counts, v_counts, u_sum, v_sum)

    return u


def _to_numpy(u_counts, v_counts, otu_ids):
    """Cast objects of importance to numpy arrays"""
    return np.asarray(u_counts), np.asarray(v_counts), np.asarray(otu_ids)


def _boundary_case(u_counts, v_counts, normalized=False, unweighted=True):
    """Test for boundary conditions

    Parameters
    ----------
    u_counts, v_counts: np.array
        Vectors of counts of OTUs for two samples.
    normalized: bool
        Indicates if the method is normalized.
    unweighted: bool
        Indicates if the method is weighted.

    Returns
    -------
    float or None
        Specifically, one of `[0.0, 1.0, None]`. `None` indicates that a
        boundary condition was not observed.

    Notes
    -----
    The following boundary conditions are tested:

        * if u_counts or v_counts is all zeros
        * if both u_counts and v_counts are all zeros
        * both of the above conditions with respect to normalized and weighted
    """
    u_sum = u_counts.sum()
    v_sum = v_counts.sum()

    if not u_sum or not v_sum:
        if u_sum + v_sum:
            # u or v counts are all zeros
            if normalized or unweighted:
                # u or v counts are all zeros
                return 1.0
            elif unweighted:
                # u and v are zero
                return 0.0

            # NOTE: we cannot handle the unnormalized case here yet as it
            # requires operations on the tree vector
        else:
            # u and v are zero
            return 0.0

    return None


def _branch_correct(tip_distances, i, j, i_sum, j_sum):
    """Calculates weighted unifrac branch length correction.

    tip_distances  must be 0 except for tips.
    """
    result = tip_distances.ravel()*((i / i_sum)+(j / j_sum))
    return result.sum()


def bind_to_parent_array(t, a):
    """Binds tree to array a, returning result in list.

    Takes as input tree t with id set.

    Returns list of (node_row, parent_row such that node_row points to the
    row of a that corresponds to the current row, and parent_row points to
    the row of the parent.

    Order will be preorder traversal, i.e. for propagating attributes from
    the root to the tip.

    Typical usage of this function is to set up an array structure for many
    preorder traversals on the same tree, especially where you plan to change
    the data between traversals.
    """
    result = []
    for n in t.traverse(self_before=True, self_after=False):
        if n is not t:
            result.append([a[n.id], a[n.parent.id]])
    return result


def tip_distances(a, bound_indices, tip_indices):
    """Sets each tip to its distance from the root.
    Note: This will need its own unittest"""
    for i, s in bound_indices:
        i += s
    mask = np.zeros(len(a))
    np.put(mask, tip_indices, 1)
    a *= mask[:, np.newaxis]
