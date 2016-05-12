# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def is_id_pair(index):
    return (
        isinstance(index, tuple) and
        len(index) == 2 and
        isinstance(index[0], str) and
        isinstance(index[1], str)
    )
