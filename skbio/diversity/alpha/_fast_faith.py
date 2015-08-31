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
from skbio.diversity._fast_base import counts_and_index
from skbio.diversity._base import (_validate_counts_vector,
                                   _validate_otu_ids_and_tree)


@experimental(as_of="0.4.0")
def faith_pd_fast(counts, otu_ids, tree, validate=True, indexed=None):
    """TODO: add or just pull from educational's __doc__?"""
    if validate:
        counts = _validate_counts_vector(counts)
        _validate_otu_ids_and_tree(counts, otu_ids, tree)

    counts = np.asarray(counts)
    otu_ids = np.asarray(otu_ids)

    if counts.sum() == 0:
        return 0.0

    # If Faith's PD is being calculated for a large number of samples, the
    # count_array could be produced a single time for the samples. This would
    # be much faster than producing it for each sample.
    count_array, indexed = counts_and_index(counts, otu_ids, tree, indexed)
    length = indexed['length']

    count_array = np.where(count_array > 0, 1, 0)
    result = (length * count_array.T).sum()

    return result
