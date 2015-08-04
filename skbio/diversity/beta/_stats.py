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
    u_observed_otus = dict(zip(otu_ids, u))
    v_observed_otus = dict(zip(otu_ids, v))
    observed_nodes1 = set(tree.observed_nodes(u_observed_otus))
    observed_nodes2 = set(tree.observed_nodes(v_observed_otus))
    observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
    shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_branch_length
    return unweighted_unifrac

def weighted_unifrac(u, v, otu_ids, tree):
    u_observed_otus = dict(zip(otu_ids, u))
    v_observed_otus = dict(zip(otu_ids, v))
    observed_nodes1 = tree.observed_nodes(u_observed_otus)
    observed_nodes2 = tree.observed_nodes(v_observed_otus)
    differential_branch_length = 0
    observed_branch_length = 0
    for o in set(observed_nodes1) | set(observed_nodes2):
        observed_count = observed_nodes1[o] + observed_nodes2[o]
        observed_count_diff = abs(observed_nodes1[o] - observed_nodes2[o])
        differential_branch_length += (o.length * observed_count_diff)
        observed_branch_length += (o.length * observed_count)
    weighted_unifrac = differential_branch_length / observed_branch_length
    return weighted_unifrac

def _to_pdist_metric(metric, **kwargs):
    result = partial(metric, **kwargs)
    return result
