# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""BIOM table specifications and utilities."""

import numpy as np
from biom import Table, example_table
from skbio.io.descriptors import Read, Write

Table.default_write_format = "biom"

# Define read and write methods for the Table class.
Table.read = Read()
Table.write = Write()


def _table_to_numpy(table):
    """Convert a Table to a dense representation.

    This is a stop-gap solution to allow current Table objects to interoperate
    with existing driver methods, until they transition to be "sparse" aware.

    """
    sample_ids = list(table.ids())
    obs_ids = list(table.ids(axis="observation"))

    if table.is_empty():
        counts = np.array([[]] * len(sample_ids))
    else:
        counts = table.matrix_data.T.toarray()

    return counts, sample_ids, obs_ids
