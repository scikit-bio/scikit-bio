# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


def phylogenetic_diversity(u, otu_ids, tree):
    observed_otus = {o: u for o, u in zip(otu_ids, u) if u >= 1}
    observed_nodes = tree.observed_node_counts(observed_otus)
    result = sum(o.length for o in observed_nodes if o.length is not None)
    return result
