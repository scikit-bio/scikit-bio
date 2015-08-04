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
    observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2 if o.length is not None)
    shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2 if o.length is not None)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_branch_length
    return unweighted_unifrac

def weighted_unifrac(u, v, otu_ids, tree):
    u_observed_otus = {o: u for o, u in zip(otu_ids, u) if u >= 1}
    v_observed_otus = {o: v for o, v in zip(otu_ids, v) if v >= 1}
    observed_nodes1 = tree.observed_node_counts(u_observed_otus)
    observed_nodes2 = tree.observed_node_counts(v_observed_otus)
    differential_branch_length = 0
    observed_branch_length = 0
    for o in set(observed_nodes1) | set(observed_nodes2):
        observed_count = observed_nodes1[o] + observed_nodes2[o]
        observed_count_diff = abs(observed_nodes1[o] - observed_nodes2[o])
        # handle the case of o.length is None
        o_length = o.length or 0
        differential_branch_length += (o_length * observed_count_diff)
        observed_branch_length += (o_length * observed_count)
    weighted_unifrac = differential_branch_length / observed_branch_length
    return weighted_unifrac

def _to_pdist_metric(metric, **kwargs):
    result = partial(metric, **kwargs)
    return result
