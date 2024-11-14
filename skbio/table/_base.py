# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from biom import Table, example_table
from skbio.io.registry import Read, Write

Table.default_write_format = "biom"

# Define read and write methods for the Table class.
Table.read = Read()
Table.write = Write()
