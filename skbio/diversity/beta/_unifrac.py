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


def _observed_otu_counts(counts, otu_ids):
    return {o: c for o, c in zip(otu_ids, counts) if c >= 1}


def _validate(u_counts, v_counts, otu_ids, tree):
    if (np.asarray(u_counts) < 0).any() or (np.asarray(v_counts) < 0).any():
        raise ValueError('Counts vectors can only contain integers >= 0.')

    len_otu_ids = len(otu_ids)
    set_otu_ids = set(otu_ids)
    if len_otu_ids != len(set_otu_ids):
        raise ValueError('Duplicated OTU ids are not allowed in otu_ids.')

    len_u_counts = len(u_counts)
    len_v_counts = len(v_counts)
    if len_u_counts != len_v_counts:
        raise ValueError("Input vectors u_counts and v_counts must be of "
                         "equal length.")
    if len_u_counts != len_otu_ids:
        raise ValueError("Input vector otu_ids must be the same length as "
                         "u_counts and v_counts.")

    branch_lengths = [e.length for e in tree.traverse() if not e.is_root()]
    if np.array([l is None for l in branch_lengths]).any():
        raise ValueError('All non-root nodes in tree must have a branch '
                         'length.')


@experimental(as_of="0.4.0-dev")
def unweighted_unifrac(u_counts, v_counts, otu_ids, tree,
                       suppress_validation=False):
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
    suppress_validation: bool, optional
        If `True`, validation of the input won't be performed. This step can
        be slow, so if validation is being before elsewhere it can be disabled
        here. However, invalid input data can lead to invalid results, so this
        step should not be bypassed all together.

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
    if not suppress_validation:
        _validate(u_counts, v_counts, otu_ids, tree)
    u_obs_otu_counts = _observed_otu_counts(u_counts, otu_ids)
    v_obs_otu_counts = _observed_otu_counts(v_counts, otu_ids)
    u_obs_nodes = set(tree.observed_node_counts(u_obs_otu_counts))
    v_obs_nodes = set(tree.observed_node_counts(v_obs_otu_counts))
    obs_branch_length = sum(o.length or 0.0 for o in u_obs_nodes | v_obs_nodes)
    if obs_branch_length == 0:
        # boundary case where both communities have no members
        return 0.0
    shared_branch_length = sum(
        o.length or 0.0 for o in u_obs_nodes & v_obs_nodes)
    unique_branch_length = obs_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / obs_branch_length
    return unweighted_unifrac


@experimental(as_of="0.4.0-dev")
def weighted_unifrac(u_counts, v_counts, otu_ids, tree, normalized=False,
                     suppress_validation=False):
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
        If ``True``, apply branch length normalization. Resulting distances
        will then be in the range ``[0, 1]``.
    suppress_validation: bool, optional
        If `True`, validation of the input won't be performed. This step can
        be slow, so if validation is being before elsewhere it can be disabled
        here. However, invalid input data can lead to invalid results, so this
        step should not be bypassed all together.

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
    Unweighted UniFrac was originally described in [1]_, which includes a
    discussion of unweighted (qualitative) versus weighted (quantitiative)
    diversity metrics. Deeper mathemtical discussions of this metric is
    presented in [2]_.

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
    if not suppress_validation:
        _validate(u_counts, v_counts, otu_ids, tree)
    u_obs_otu_counts = _observed_otu_counts(u_counts, otu_ids)
    u_total_count = sum(u_counts)
    v_obs_otu_counts = _observed_otu_counts(v_counts, otu_ids)
    v_total_count = sum(v_counts)
    u_obs_nodes = tree.observed_node_counts(u_obs_otu_counts)
    v_obs_nodes = tree.observed_node_counts(v_obs_otu_counts)
    uv_obs_nodes = set(u_obs_nodes) | set(v_obs_nodes)
    if len(uv_obs_nodes) == 0:
        # boundary case where both communities have no members
        return 0.0
    weighted_unifrac = 0
    D = 0
    for o in uv_obs_nodes:
        # handle the case of o.length is None
        b = o.length or 0
        u_branch_weight = _sample_branch_weight(u_obs_nodes[o], u_total_count)
        v_branch_weight = _sample_branch_weight(v_obs_nodes[o], v_total_count)
        branch_weight = abs(u_branch_weight - v_branch_weight)
        weighted_unifrac += b * branch_weight
        if normalized and o.is_tip():
            d = o.accumulate_to_ancestor(tree.root())
            normed_weight = u_branch_weight + v_branch_weight
            D += (d * normed_weight)
    if normalized:
        return weighted_unifrac / D
    else:
        return weighted_unifrac


def _sample_branch_weight(observed_nodes, total_count):
    if total_count == 0:
        return 0
    else:
        return observed_nodes / total_count
