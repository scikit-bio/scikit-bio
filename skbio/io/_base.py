# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range


def _chunk_str(s, n, char):
    """Insert `char` character every `n` characters in string `s`.

    Canonically pronounced "chunkster".

    """
    # Modified from http://stackoverflow.com/a/312464/3776794
    if n < 1:
        raise ValueError(
            "Cannot split string into chunks with n=%d. n must be >= 1." % n)
    return char.join((s[i:i+n] for i in range(0, len(s), n)))


def _int_to_ordinal_str(n):
    """Return ordinal string version of int `n`.

    For example:

    1 -> '1st'
    2 -> '2nd'
    3 -> '3rd'
    ...

    """
    # Taken from http://stackoverflow.com/a/20007730/3776794
    # Originally from http://codegolf.stackexchange.com/a/4712 by Gareth
    if n < 0:
        raise ValueError("Cannot convert negative integer %d to ordinal "
                         "string." % n)
    return "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
