# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from functools import partial


def unweighted_unifrac(u, v, otu_ids, tree):
    """ Compute unweighted unifrac

    Parameters
    ----------
    u, v: list, np.array
        Vectors of counts of OTUs for two samples. Must be equal length.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u`` and ``v``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.

    Returns
    -------
    float
        The unweighted unifrac distance between the two samples.

    Raises
    ------
    ValueError
        If ``u``, ``v``, and ``otu_ids`` are not all equal in length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    See Also
    --------
    weighted_unifrac

    Notes
    -----
    Unweighted unifrac was originally described in [1]_. A discussion of
    unweighted (qualitative) versus weighted (quantitiative) diversity metrics
    is presented in [2]_.

    References
    ----------
    .. [1] Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for
       comparing microbial communities. Appl. Environ. Microbiol. 71, 8228–8235
       (2005).

    .. [2] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576–1585 (2007).

    """
    if len(u) != len(v) != len(otu_ids):
        raise ValueError("Input vectors u, v, and otu_ids must be equal"
                         " length.")
    u_observed_otus = {o: u for o, u in zip(otu_ids, u) if u >= 1}
    v_observed_otus = {o: v for o, v in zip(otu_ids, v) if v >= 1}
    observed_nodes1 = set(tree.observed_node_counts(u_observed_otus))
    observed_nodes2 = set(tree.observed_node_counts(v_observed_otus))
    observed_branch_length = sum(
        o.length for o in observed_nodes1 | observed_nodes2
        if o.length is not None)
    if observed_branch_length == 0:
        # boundary case where both communities have no members
        return 0.0
    shared_branch_length = sum(
        o.length for o in observed_nodes1 & observed_nodes2
        if o.length is not None)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_branch_length
    return unweighted_unifrac


def _sample_branch_weight(observed_nodes, total_count):
    if total_count == 0:
        return 0
    else:
        return observed_nodes / total_count


def weighted_unifrac(u, v, otu_ids, tree, normalized=False):
    """ Compute weighted unifrac with or without branch length normalization

    Parameters
    ----------
    u, v: list, np.array
        Vectors of counts of OTUs for two samples. Must be equal length.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u`` and ``v``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    normalized: boolean, optional
        If ``True``, apply branch length normalization. Resulting distances
        will then be in the range ``[0, 1]``.

    Returns
    -------
    float
        The weighted unifrac distance between the two samples.

    Raises
    ------
    ValueError
        If ``u``, ``v``, and ``otu_ids`` are not all equal in length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    See Also
    --------
    unweighted_unifrac

    Notes
    -----
    Unweighted unifrac was originally described in [1]_, which includes a
    discussion of unweighted (qualitative) versus weighted (quantitiative)
    diversity metrics.

    References
    ----------
    .. [1] Lozupone, C. A., Hamady, M., Kelley, S. T. & Knight, R. Quantitative
       and qualitative beta diversity measures lead to different insights into
       factors that structure microbial communities. Appl. Environ. Microbiol.
       73, 1576–1585 (2007).

    """
    u_observed_otus = {o: u for o, u in zip(otu_ids, u) if u >= 1}
    u_total_count = sum(u)
    v_observed_otus = {o: v for o, v in zip(otu_ids, v) if v >= 1}
    v_total_count = sum(v)
    if u_total_count == 0 and v_total_count == 0:
        # boundary case where both communities have no members
        return 0.0
    u_observed_nodes = tree.observed_node_counts(u_observed_otus)
    v_observed_nodes = tree.observed_node_counts(v_observed_otus)
    weighted_unifrac = 0
    D = 0
    for o in set(u_observed_nodes) | set(v_observed_nodes):
        # handle the case of o.length is None
        b = o.length or 0
        u_branch_weight = _sample_branch_weight(
            u_observed_nodes[o], u_total_count)
        v_branch_weight = _sample_branch_weight(
            v_observed_nodes[o], v_total_count)
        branch_weight = abs(u_branch_weight - v_branch_weight)
        weighted_unifrac += b * branch_weight
        if o.is_tip() and normalized:
            d = o.accumulate_to_ancestor(tree.root())
            normed_weight = u_branch_weight + v_branch_weight
            D += (d * normed_weight)
    if normalized:
        return weighted_unifrac / D
    else:
        return weighted_unifrac


def _to_pdist_metric(metric, **kwargs):
    result = partial(metric, **kwargs)
    return result
