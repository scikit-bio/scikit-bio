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
