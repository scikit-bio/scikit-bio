# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


def phylogenetic_diversity(u, otu_ids, tree, min_count=1):
    observed_otus = []
    for count, id_ in zip(u, otu_ids):
        if count >= min_count:
            observed_otus.append(id_)
    observed_nodes = tree.observed_nodes(observed_otus)
    result = sum(o.length for o in observed_nodes)
    return result
