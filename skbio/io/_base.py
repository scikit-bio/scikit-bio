# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


def _chunk_str(s, n, char):
    """Insert `char` character every `n` characters in string `s`."""
    # Modified from http://stackoverflow.com/a/312464/3776794
    return char.join((s[i:i+n] for i in range(0, len(s), n)))
