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
from skbio.diversity._fast_base import (_fast_unifrac_setup, bind_to_array,
                                        bool_descendants,
                                        _skbio_counts_to_envs, traverse_reduce)


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
                            validate=False, **kwargs):
    """fit to unifrac API"""
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    u_sum = sum(u_counts)
    v_sum = sum(v_counts)

    if not u_sum or not v_sum:
        if u_sum + v_sum:
            # u or v counts are all zeros
            return 1.0
        else:
            # u and v are zero
            return 0.0

    envs = _skbio_counts_to_envs(otu_ids, u_counts, v_counts)

    return fast_unifrac(tree, envs)


@experimental(as_of="0.4.0")
def weighted_unifrac_fast(u_counts, v_counts, otu_ids, tree, normalized=False,
                          validate=False, **kwargs):
    u_sum = sum(u_counts)
    v_sum = sum(v_counts)

    if not u_sum or not v_sum:
        if u_sum + v_sum:
            # u or v counts are all zeros
            if normalized:
                return 1.0
            # NOTE: we cannot handle the unnormalized case here yet as it
            # requires operations on the tree vector
        else:
            # u and v are zero
            return 0.0
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)
    envs = _skbio_counts_to_envs(otu_ids, u_counts, v_counts)

    return fast_unifrac(tree, envs, weighted=True, normalized=normalized)


def fast_unifrac(t, envs, weighted=False, metric=unifrac, normalized=False):
    # weighted_unifrac_f=_weighted_unifrac,make_subtree=True):
    """
    Run fast unifrac.

    Parameters
    ----------
    t: skbio.TreeNode
        phylogenetic tree relating the sequences.
    envs: dict
        dict of {sequence:{env:count}} showing environmental abundance.
    weighted: bool
        if True, performs the weighted UniFrac procedure.
    metric: function
        distance metric to use.  currently you must use unifrac only
        if weighted=True.
        see fast_tree.py for metrics (e.g.: G, unnormalized_G, unifrac, etc.)

    Returns
    -------
    u : np.array
        Unifrac distance matrix
    """
    (envs, count_array,
     unique_envs, env_to_index,
     node_to_index, env_names,
     branch_lengths, nodes, t) = _fast_unifrac_setup(t, envs)

    bound_indices = bind_to_array(nodes, count_array)
    # weighted unifrac
    if weighted:
        tip_indices = [n._leaf_index for n in t.tips()]
        sum_descendants(bound_indices)
        tip_ds = branch_lengths.copy()[:, np.newaxis]
        bindings = bind_to_parent_array(t, tip_ds)
        tip_distances(tip_ds, bindings, tip_indices)

        # NOTE: in the unnormalized case, it is possible for count_array to
        # have a single column of v_counts come in as all zeros.
        u_counts = count_array[:, 0]
        if count_array.shape[1] == 1:
            v_counts = np.zeros(count_array.shape[0])
        else:
            v_counts = count_array[:, 1]
        u_sum = (np.take(u_counts, tip_indices)).sum()
        v_sum = (np.take(v_counts, tip_indices)).sum()

        u = w_unifrac(branch_lengths, tip_indices, u_counts, v_counts)
        if normalized:
            u /= _branch_correct(tip_ds, u_counts, v_counts, u_sum, v_sum)
    # unweighted unifrac
    else:
        bool_descendants(bound_indices)
        u = unifrac(branch_lengths, count_array[:, 0], count_array[:, 1])

    return u


def _branch_correct(tip_distances, i, j, i_sum, j_sum):
    """Calculates weighted unifrac branch length correction.

    tip_distances  must be 0 except for tips.
    """
    result = tip_distances.ravel()*((i / i_sum)+(j / j_sum))
    return result.sum()


# specific to weighted
def sum_descendants(bound_indices):
    """For each internal node, sets col to sum of values in descendants."""
    traverse_reduce(bound_indices, sum)


def bind_to_parent_array(t, a):
    """Binds tree to array a, returning result in list.

    Takes as input tree t with _leaf_index set.

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
            result.append([a[n._leaf_index], a[n.parent._leaf_index]])
    return result


def unifrac_tasks_from_matrix(u, env_names):
    """Returns the UniFrac matrix, PCoA, and/or cluster from the matrix."""
    result = {}
    result['distance_matrix'] = (u, env_names)
    return result


def tip_distances(a, bound_indices, tip_indices):
    """Sets each tip to its distance from the root.
    Note: This will need its own unittest"""
    for i, s in bound_indices:
        i += s
    mask = np.zeros(len(a))
    np.put(mask, tip_indices, 1)
    a *= mask[:, np.newaxis]
