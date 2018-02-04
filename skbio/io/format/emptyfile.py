r"""
Empty Files (:mod:`skbio.io.format.emptyfile`)
==============================================

.. currentmodule:: skbio.io.format.emptyfile

This format exists to make debugging simpler, often an empty file is a mistake
which can take an embarrasing amount of time to notice. This format has only
a sniffer and no readers or writers, so error messages will indicate as such
if an empty file is accidentally used as input.

Format Support
--------------
**Has Sniffer: Yes**

Format Specification
--------------------
An empty file consists of only whitespace characters.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.io import create_format

emptyfile = create_format('<emptyfile>')


@emptyfile.sniffer()
def _empty_file_sniffer(fh):
    for line in fh:
        if line.strip():
            return False, {}
    return True, {}
