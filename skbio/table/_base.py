# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from biom import Table, example_table
from skbio.io.util import DocDescriptor

Table.default_write_format = "biom"

read = DocDescriptor("read")
write = DocDescriptor("write")
