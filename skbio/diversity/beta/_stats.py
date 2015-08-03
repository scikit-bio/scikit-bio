# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from functools import partial

def unweighted_unifrac(u, v, otu_ids, tree, min_count=1):
    u_observed_otus = []
    v_observed_otus = []
    for u_count, v_count, id_ in zip(u, v, otu_ids):
        if u_count >= min_count:
            u_observed_otus.append(id_)
        if v_count >= min_count:
            v_observed_otus.append(id_)
    observed_nodes1 = tree.observed_nodes(u_observed_otus)
    observed_nodes2 = tree.observed_nodes(v_observed_otus)
    observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
    shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_branch_length
    return unweighted_unifrac

def _to_pdist_metric(metric, **kwargs):
    result = partial(metric, **kwargs)
    return result
