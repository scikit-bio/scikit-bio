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
from skbio.diversity._fast_base_cy import tip_distances


def unifrac(m, i, j):
    """Calculates unifrac(i,j) from m.

    Parameters
    ----------
    m : np.array
        A 1D vector that represents the lengths in the tree in postorder.
    i,j : np.array
        A slice of states from m in postorder.

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
    _or = np.logical_or(i, j),
    _and = np.logical_and(i, j)
    return 1 - ((m * _and).sum() / (m * _or).sum())


def w_unifrac(m, i, j, i_sum, j_sum):
    """Calculates weighted unifrac(i, j) from m

    Parameters
    ----------
    m : np.array
        A 1D vector that represents the lengths in the tree in postorder.
    i, j : np.array
        A slice of states from m in postorder.
    i_sum, j_sum: float
        Counts of the observations in each environment

    Returns
    -------
    float
        The weighted unifrac score
    """
    if i_sum:
        i_ = i / i_sum
    else:
        i_ = 0.0

    if j_sum:
        j_ = j / j_sum
    else:
        j_ = 0.0

    return (m * abs(i_ - j_)).sum()


def make_pdist(counts, obs_ids, tree, indexed=None, metric=unifrac,
               normalized=False, **kwargs):
    """Construct a metric for use with skbio.diversity.beta.pw_distances

    Parameters
    ----------
    counts : np.ndarray
        A matrix of environment counts where each row is an environment, and
        each column is an observation. The observations are expected to be in
        index order w.r.t. obs_ids, but do not need to be in tree order.
    obs_ids : np.ndarray
        A vector of observation IDs. These IDs must map to tips in the tree.
    tree : skbio.tree.TreeNode
        A tree that represents the relationships between the observations.
    indexed : dict, optiona;
        The result of skbio.diversity.beta.index_tree
    metric : function, {unweighted_unifrac_fast, weighted_unifrac_fast}
        The specific metric to use.
    normalized : bool
        Whether the weighted unifrac calculation should be normalized.

    Notes
    -------
    example usage

    metric, counts, length = make_unweighted_pdist(input_counts, tip_ids, tree)
    mat = pw_distances(metric, counts, ids=['%d' % i for i in range(10)])
    """
    count_array, indexed = _counts_and_length(counts, obs_ids, tree, indexed)
    length = indexed['length']

    if metric is unifrac:
        def f(u, v):
            boundary = _boundary_case(u.sum(), v.sum())
            if boundary is not None:
                return boundary
            return unifrac(length, u, v)

    elif metric is w_unifrac:
        # This block is duplicated in weighted_unifrac_fast -- possibly should
        # be decomposed.
        # There is a lot in common with both of these methods, but pulling the
        # normalized check out reduces branching logic.
        if normalized:
            def f(u, v):
                u_sum = u.sum()
                v_sum = v.sum()

                boundary = _boundary_case(u_sum, v_sum)
                if boundary is not None:
                    return boundary

                tip_idx = np.array([n.id for n in indexed['id_index'].values()
                                    if n.is_tip()])
                tip_ds = tip_distances(length, tree, tip_idx)

                u = w_unifrac(length, u, v, u_sum, v_sum)
                u /= _branch_correct(tip_ds, u, v, u_sum, v_sum)
                return u
        else:
            def f(u, v):
                u_sum = u.sum()
                v_sum = v.sum()

                boundary = _boundary_case(u_sum, v_sum)
                if boundary is not None:
                    return boundary

                u = w_unifrac(length, u, v, u_sum, v_sum)
                return u
    else:
        raise AttributeError("Unknown metric: %s" % metric)

    return f, count_array.T, length


@experimental(as_of="0.4.0")
def unweighted_unifrac_fast(u_counts, v_counts, otu_ids, tree,
                            validate=True, indexed=None, **kwargs):
    """TODO: add or just pull from educational's __doc__?"""
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    # cast to numpy types
    u_counts = np.asarray(u_counts)
    v_counts = np.asarray(v_counts)
    otu_ids = np.asarray(otu_ids)

    boundary = _boundary_case(u_counts.sum(), v_counts.sum())
    if boundary is not None:
        return boundary

    # aggregate state information up the tree (stored in counts_array), and
    # retrieve the aggregated state information for each input count vector
    counts = np.vstack([u_counts, v_counts])
    count_array, indexed = _counts_and_length(counts, otu_ids, tree, indexed)
    u_counts = count_array[:, 0]
    v_counts = count_array[:, 1]

    length = indexed['length']

    return unifrac(length, u_counts, v_counts)


@experimental(as_of="0.4.0")
def weighted_unifrac_fast(u_counts, v_counts, otu_ids, tree, normalized=False,
                          validate=True, indexed=None, **kwargs):
    """TODO: add or just pull from educational's __doc__?"""
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    # convert to numpy types
    u_counts = np.asarray(u_counts)
    v_counts = np.asarray(v_counts)
    otu_ids = np.asarray(otu_ids)

    u_sum = u_counts.sum()
    v_sum = v_counts.sum()

    # check boundary conditions and shortcut if necessary
    boundary = _boundary_case(u_sum, v_sum, normalized, unweighted=False)
    if boundary is not None:
        return boundary

    # aggregate state information up the tree (stored in counts_array), and
    # retrieve the aggregated state information for each input count vector
    counts = np.vstack([u_counts, v_counts])
    count_array, indexed = _counts_and_length(counts, otu_ids, tree, indexed)
    u_counts = count_array[:, 0]
    v_counts = count_array[:, 1]

    # fetch the lengths
    length = indexed['length']

    u = w_unifrac(length, u_counts, v_counts, u_sum, v_sum)
    if normalized:
        # get the index positions for tips in counts_array, and determine the
        # tip distances to the root
        tip_indices = np.array([n.id for n in indexed['id_index'].values()
                                if n.is_tip()])

        tip_ds = tip_distances(length, tree, tip_indices)
        u /= _branch_correct(tip_ds, u_counts, v_counts, u_sum, v_sum)

    return u


def _boundary_case(u_sum, v_sum, normalized=False, unweighted=True):
    """Test for boundary conditions

    Parameters
    ----------
    u_sum, v_sum: float
        The sum of the observations in both environments.
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

        * if u_sum or v_sum are zero
        * if both u_sum and v_sum are zero
        * both of the above conditions with respect to normalized and weighted
    """
    if u_sum and v_sum:
        return None

    if u_sum + v_sum:
        # u or v counts are all zeros
        if unweighted or normalized:
            # u or v counts are all zeros
            return 1.0
        # NOTE: we cannot handle the unnormalized case here yet as it
        # requires operations on the tree vector
    else:
        # u and v are zero
        return 0.0

    return None


def _branch_correct(tip_dists, i, j, i_sum, j_sum):
    """Calculates weighted unifrac branch length correction.

    Parameters
    ----------
    tip_dists : np.ndarray
        1D column vector of branch lengths in post order form. Only tips
        should be non-zero as this represents the distance from each tip to
        root.
    i, j : np.ndarray
        Aggregated environment counts. This vector is expected to be in index
        order with tip_dists
    i_sum, j_sum: float
        Counts of the observations in each environment

    Returns
    -------
    np.ndarray
        The corrected branch lengths
    """
    return (tip_dists.ravel() * ((i / i_sum) + (j / j_sum))).sum()


# TODO: define a block pdist method such that "blocks" of the resulting
# distance matrix can be computed independently. This would allow for an
# embarassingly parallel strategy for computing large distance matrices. This
# method needs to be paired with a DistanceMatrix assembly method, however,
# it may be beneficial to "assemble" directly into a memory mapped data
# structure.


def block_pdist(observation_data):
    """Compute a distance matrix in parallel

    Parameters
    ----------
    observation_data : scipy.sparse.spmatrix
        The observation data where the columns are the samples and the rows
        are the observations

    Notes
    -----
    If the resulting distance matrix looked like the following, where the
    IDs of the matrix are [1, 2, 3, 4]:

      1 2 3 4
    1 0 A B C
    2 A 0 D E
    3 B D 0 F
    4 C E F 0

    The goal of block_pdist is to compute allow "blocks" of the distance matrix
    to be computed independently. For instance, the above distnace matrix
    could be computed by just valuating the following "blocks":

    #### DO LOWER TRI INSTEAD


    B C  0 A  0 F
    D E, A 0, F 0

    This translates into computing the distances between the following IDs:


    #### DO LOWER TRI INSTEAD
    B C
    D E [(1, 3), (1, 4), (2, 3), (2, 4)]

    0 A
    A 0 [(1, 1), (1, 2), (2, 1), (2, 2)]

    0 F
    F 0 [(3, 3), (3, 4), (4, 3), (4, 4)]

    There is an easily (but minor) optimization which is that we do not need
    to compute d(u, u), and if we have computed d(u, v) in a block, then it is
    unnecessary to compute d(v, ).
    """
    pass


def _block_coordinates(n, k):
    """Yield block coordinates for block_pdist

    Parameters
    ----------
    n : int
        The size of the resulting distance matrix. n > 0
    k : int
        The size of a block. 0 < k <= n

    Notes
    -----
    Does not cross the diagonal for points yielded

    Does not yield the diagonal

    k does not need to divide evenly into n. Any remainder is returned. For
    instance, if the resulting distance matrix were 4 x 4:

    0 A B C D E
    A 0 E F G H
    B E 0 I J K
    C F I 0 L M
    D G J L 0 N
    E H K M N 0

    and k were 4, the following blocks would be yielded (values are index
    positions i,j):


    First set:

        [(5, 0), (5, 1), (5, 2), (5, 3),
         (4, 0), (4, 1), (4, 2), (4, 3),
         (3, 0), (3, 1), (3, 2), (3, 3),
         (2, 0), (2, 1)]

    Second set:

        [(1, 0)]

    """
    pass
