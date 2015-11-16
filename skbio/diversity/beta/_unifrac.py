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
from skbio.diversity._driver import (_validate_counts_vectors,
                                     _validate_otu_ids_and_tree,
                                     _vectorize_counts_and_tree)
from skbio.diversity._base_cy import _tip_distances


@experimental(as_of="0.4.0-dev")
def unweighted_unifrac(u_counts, v_counts, otu_ids, tree, validate=True):
    """ Compute unweighted UniFrac

    Parameters
    ----------
    u_counts, v_counts: list, np.array
        Vectors of counts of OTUs for two samples. Must be equal length.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u_counts`` and ``v_counts``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results, so this step
        should not be bypassed all together.

    Returns
    -------
    float
        The unweighted UniFrac distance between the two samples.

    Raises
    ------
    ValueError
        If ``u_counts``, ``v_counts``, and ``otu_ids`` are not all of equal
        length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    See Also
    --------
    weighted_unifrac

    Notes
    -----
    Unweighted UniFrac was originally described in [1]_. A discussion of
    unweighted (qualitative) versus weighted (quantitiative) diversity metrics
    is presented in [2]_. Deeper mathemtical discussions of this metric is
    presented in [3]_.

    This implementation differs from that in PyCogent (and therefore QIIME
    versions less than 2.0.0) by imposing a few additional restrictions on the
    inputs. First, the input tree must be rooted. In PyCogent, if an unrooted
    tree was provided that had a single trifurcating node (a newick convention
    for unrooted trees) that node was considered the root of the tree. Next,
    all OTU IDs must be tips in the tree. PyCogent would silently ignore OTU
    IDs that were not present the tree. To reproduce UniFrac results from
    PyCogent with scikit-bio, ensure that your PyCogent UniFrac calculations
    are performed on a rooted tree and that all OTU IDs are present in the
    tree.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-8235
       (2005).

    .. [2] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).

    .. [3] Lozupone, C., Lladser, M. E., Knights, D., Stombaugh, J. & Knight,
       R. UniFrac: an effective distance metric for microbial community
       comparison. ISME J. 5, 169-172 (2011).

    """
    normalized = False
    unweighted = True
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    boundary = _boundary_case(sum(u_counts), sum(v_counts),
                              normalized=normalized, unweighted=unweighted)
    if boundary is not None:
        return boundary

    u_node_counts, v_node_counts, u_total_count, v_total_count, tree_index =\
        _setup_single_unifrac(u_counts, v_counts, otu_ids, tree, validate,
                              normalized=normalized, unweighted=unweighted)
    return _unweighted_unifrac(u_node_counts, v_node_counts,
                               tree_index['length'])


@experimental(as_of="0.4.0-dev")
def weighted_unifrac(u_counts, v_counts, otu_ids, tree, normalized=False,
                     validate=True):
    """ Compute weighted UniFrac with or without branch length normalization

    Parameters
    ----------
    u_counts, v_counts: list, np.array
        Vectors of counts of OTUs for two samples. Must be equal length.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u_counts`` and ``v_counts``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    normalized: boolean, optional
        If ``True``, apply branch length normalization, which is described in
        [1]_. Resulting distances will then be in the range ``[0, 1]``.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results, so this step
        should not be bypassed all together.

    Returns
    -------
    float
        The weighted UniFrac distance between the two samples.

    Raises
    ------
    ValueError
        If ``u_counts``, ``v_counts``, and ``otu_ids`` are not all equal in
        length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    See Also
    --------
    unweighted_unifrac

    Notes
    -----
    Weighted UniFrac was originally described in [1]_, which includes a
    discussion of unweighted (qualitative) versus weighted (quantitiative)
    diversity metrics. Deeper mathemtical discussions of this metric is
    presented in [2]_.

    This implementation differs from that in PyCogent (and therefore QIIME
    versions less than 2.0.0) by imposing a few additional restrictions on the
    inputs. First, the input tree must be rooted. In PyCogent, if an unrooted
    tree was provided that had a single trifurcating node (a newick convention
    for unrooted trees) that node was considered the root of the tree. Next,
    all OTU IDs must be tips in the tree. PyCogent would silently ignore OTU
    IDs that were not present the tree. To reproduce UniFrac results from
    PyCogent with scikit-bio, ensure that your PyCogent UniFrac calculations
    are performed on a rooted tree and that all OTU IDs are present in the
    tree.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576-1585 (2007).

    .. [2] Lozupone, C., Lladser, M. E., Knights, D., Stombaugh, J. & Knight,
       R. UniFrac: an effective distance metric for microbial community
       comparison. ISME J. 5, 169-172 (2011).

    """
    unweighted = False
    if validate:
        _validate(u_counts, v_counts, otu_ids, tree)

    boundary = _boundary_case(sum(u_counts), sum(v_counts),
                              normalized=normalized, unweighted=unweighted)
    if boundary is not None:
        return boundary

    u_node_counts, v_node_counts, u_total_count, v_total_count, tree_index =\
        _setup_single_unifrac(u_counts, v_counts, otu_ids, tree, validate,
                              normalized=normalized, unweighted=unweighted)

    if normalized:
        return _weighted_unifrac_normalized(u_node_counts, v_node_counts,
                                            u_total_count, v_total_count,
                                            tree, tree_index)
    else:
        return _weighted_unifrac(u_node_counts, v_node_counts,
                                 u_total_count, v_total_count,
                                 tree_index['length'])


def _validate(u_counts, v_counts, otu_ids, tree):
    _validate_counts_vectors(u_counts, v_counts, suppress_cast=True)
    _validate_otu_ids_and_tree(counts=u_counts, otu_ids=otu_ids, tree=tree)


def _setup_single_unifrac(u_counts, v_counts, otu_ids, tree, validate,
                          normalized, unweighted):

    # temporarily store u_counts and v_counts in a 2D array as that's what
    # _vectorize_counts_and_tree takes
    u_counts = np.asarray(u_counts)
    v_counts = np.asarray(v_counts)
    counts = np.vstack([u_counts, v_counts])
    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)
    # unpack counts vectors for single pairwise UniFrac calculation
    u_node_counts = counts_by_node[0]
    v_node_counts = counts_by_node[1]

    u_total_count = u_counts.sum()
    v_total_count = v_counts.sum()

    return (u_node_counts, v_node_counts, u_total_count, v_total_count,
            tree_index)


def _unweighted_unifrac(u_node_counts, v_node_counts, branch_lengths):
    """
    Parameters
    ----------
    u_node_counts, v_node_counts : np.array
        Vectors indicating presense (value greater than zero) and absense
        (value equal to zero) of nodes in two samples, `u` and `v`. Order is
        assumed to be the same as in `branch_lengths`.
    branch_lengths : np.array
        Vector of branch lengths of all nodes (tips and internal nodes) in
        postorder representation of their tree.

    Returns
    -------
    float
        Unweighted UniFrac distance between samples.

    Notes
    -----
    The count vectors passed here correspond to all nodes in the tree, not
    just the tips.

    """
    unique = np.logical_xor(u_node_counts, v_node_counts)
    observed = np.logical_or(u_node_counts, v_node_counts)
    return (branch_lengths * unique).sum() / (branch_lengths * observed).sum()


def _weighted_unifrac(u_node_counts, v_node_counts, u_total_count,
                      v_total_count, branch_lengths):
    """
    Parameters
    ----------
    u_node_counts, v_node_counts : np.array
        Vectors indicating presense (value greater than zero) and absense
        (value equal to zero) of nodes in two samples, `u` and `v`. Order is
        assumed to be the same as in `branch_lengths`.
    u_total_count, v_total_counts : int
        The sum of ``u_node_counts`` and ``v_node_counts`` vectors,
        respectively. This could be computed internally, but since this is a
        private method and the calling function has already generated these
        values, this saves an iteration over each of these vectors.
    branch_lengths : np.array
        Vector of branch lengths of all nodes (tips and internal nodes) in
        postorder representation of their tree.

    Returns
    -------
    float
        Weighted UniFrac distance between samples.

    Notes
    -----
    The count vectors passed here correspond to all nodes in the tree, not
    just the tips.

    """
    if u_total_count > 0:
        # convert to relative abundances if there are any counts
        u_ = u_node_counts / u_total_count
    else:
        # otherwise, we'll just do the computation with u_node_counts, which
        # is necessarily all zeros
        u_ = u_node_counts

    if v_total_count > 0:
        v_ = v_node_counts / v_total_count
    else:
        v_ = v_node_counts

    return (branch_lengths * abs(u_ - v_)).sum()


def _weighted_unifrac_normalized(u_node_counts, v_node_counts, u_total_count,
                                 v_total_count, tree, tree_index):
    """
    Parameters
    ----------
    u_node_counts, v_node_counts : np.array
         Vectors indicating presense (value greater than zero) and absense
         (value equal to zero) of nodes in two samples, `u` and `v`. Order is
         assumed to be the same as in `branch_lengths`.
    u_total_count, v_total_counts : int
         The sum of ``u_node_counts`` and ``v_node_counts`` vectors,
         respectively. This could be computed internally, but since this is a
         private method and the calling function has already generated these
         values, this saves an iteration over each of these vectors.
    tree: skbio.TreeNode
         Tree relating the OTUs.
    tree_index:
         Result of calling ``tree.to_array``. This is expensive to compute, so
         it's provided as a parameter here as the calling function will have
         already generated it.

    Returns
    -------
    float
        Normalized weighted UniFrac distance between samples.

    Notes
    -----
    The count vectors passed here correspond to all nodes in the tree, not
    just the tips.

    """
    branch_lengths = tree_index['length']
    u = _weighted_unifrac(u_node_counts, v_node_counts, u_total_count,
                          v_total_count, branch_lengths)
    # get the index positions for tips in counts_array, and determine the
    # tip distances to the root
    tip_indices = _get_tip_indices(tree_index)
    node_to_root_distances = _tip_distances(branch_lengths, tree, tip_indices)
    u /= _weighted_unifrac_branch_correction(
        node_to_root_distances, u_node_counts, v_node_counts, u_total_count,
        v_total_count)

    return u


def _unweighted_unifrac_pdist_f(counts, otu_ids, tree):
    """ Create optimized pdist-compatible unweighted UniFrac function

    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u_counts`` and ``v_counts``. These IDs do not need to
        be in tip order with respect to the tree.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.

    Returns
    -------
    function
        Optimized pairwise unweighted UniFrac calculator that can be passed
        to ``scipy.spatial.distance.pdist``.
    2D np.array of ints, floats
        Counts of all nodes in ``tree``.
    2D np.array of floats
        Branch length leading to all nodes in ``tree``.

    """
    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)

    def f(u_node_counts, v_node_counts):
        boundary = _boundary_case(u_node_counts.sum(), v_node_counts.sum())
        if boundary is not None:
            return boundary
        return _unweighted_unifrac(u_node_counts, v_node_counts,
                                   branch_lengths)

    return f, counts_by_node, branch_lengths


def _weighted_unifrac_pdist_f(counts, otu_ids, tree, normalized):
    """ Create optimized pdist-compatible weighted UniFrac function

    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u_counts`` and ``v_counts``. These IDs do not need to
        be in tip order with respect to the tree.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.

    Returns
    -------
    function
        Optimized pairwise unweighted UniFrac calculator that can be passed
        to ``scipy.spatial.distance.pdist``.
    2D np.array of ints, floats
        Counts of all nodes in ``tree``.
    2D np.array of floats
        Branch length leading to all nodes in ``tree``.

    """
    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)
    tip_indices = _get_tip_indices(tree_index)

    if normalized:
        node_to_root_distances = _tip_distances(branch_lengths, tree,
                                                tip_indices)

    def f(u_node_counts, v_node_counts):
        u_total_count = np.take(u_node_counts, tip_indices).sum()
        v_total_count = np.take(v_node_counts, tip_indices).sum()

        boundary = _boundary_case(u_total_count, v_total_count,
                                  unweighted=False)
        if boundary is not None:
            return boundary

        u = _weighted_unifrac(u_node_counts, v_node_counts, u_total_count,
                              v_total_count, branch_lengths)
        if normalized:
            u /= _weighted_unifrac_branch_correction(
                    node_to_root_distances, u_node_counts, v_node_counts,
                    u_total_count, v_total_count)
        return u

    return f, counts_by_node, branch_lengths


def _boundary_case(u_total_count, v_total_count, normalized=False,
                   unweighted=True):
    """Test for boundary conditions

    Parameters
    ----------
    u_total_count, v_total_count : float
        The sum of the observations in samples ``u`` and ``v``, respectively.
    normalized : bool, optional
        Indicates if the method is normalized.
    unweighted : bool, optional
        Indicates if the method is unweighted.

    Returns
    -------
    float or None
        Specifically, one of `[0.0, 1.0, None]`. `None` indicates that a
        boundary condition was not observed, and the float values are distances
        between ``u`` and ``v``.

    Notes
    -----
    The following boundary conditions are tested:

        * if u_total_count or v_total_count are zero
        * if both u_total_count and v_total_count are zero
    """
    if u_total_count and v_total_count:
        return None

    if u_total_count + v_total_count:
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


def _get_tip_indices(tree_index):
    tip_indices = np.array([n.id for n in tree_index['id_index'].values()
                            if n.is_tip()])
    return tip_indices


def _weighted_unifrac_branch_correction(node_to_root_distances, u_node_counts,
                                        v_node_counts, u_total_count,
                                        v_total_count):
    """Calculates weighted unifrac branch length correction.

    Parameters
    ----------
    node_to_root_distances : np.ndarray
        1D column vector of branch lengths in post order form. There should be
        positions in this vector for all nodes in the tree, but only tips
        should be non-zero.
    u_node_counts, v_node_counts : np.ndarray
        Counts of observations of all nodes in the tree in samples ``u`` and
        ``v``, respectively.
    u_total_count, v_total_count : float
        The sum of the observations in samples ``u`` and ``v``, respectively.

    Returns
    -------
    np.ndarray
        The corrected branch lengths
    """
    return (node_to_root_distances.ravel() *
            ((u_node_counts / u_total_count) +
             (v_node_counts / v_total_count))).sum()
