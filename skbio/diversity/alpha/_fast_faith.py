# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util._decorator import experimental
from skbio.diversity._fast_base import (_fast_unifrac_setup, bind_to_array,
                                        bool_descendants,
                                        _skbio_counts_to_envs)
from skbio.diversity._base import (_validate_counts_vector,
                                   _validate_otu_ids_and_tree)


def PD_whole_tree(t, envs):
    """Run PD on t and envs for each env.

    Note: this is specific for PD per se, use PD_generic_whole_tree if you
    want to calculate a related metric.
    """
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)
    count_array = count_array.astype(bool)
    bound_indices = bind_to_array(nodes, count_array)
    # initialize result
    bool_descendants(bound_indices)
    result = (branch_lengths * count_array.T).sum(1)
    return unique_envs, result


@experimental(as_of="0.4.0")
def faith_pd_fast(counts, otu_ids, tree, validate=True):
    """skbio api"""
    if validate:
        counts = _validate_counts_vector(counts)
        _validate_otu_ids_and_tree(counts, otu_ids, tree)

    if sum(counts) == 0:
        return 0.0

    envs = _skbio_counts_to_envs(otu_ids, counts)

    return PD_whole_tree(tree, envs)[1][0]  # just return the result
