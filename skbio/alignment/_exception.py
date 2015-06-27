# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


class SequenceCollectionError(Exception):
    """General error for sequence collection validation failures."""
    pass


class AlignmentError(SequenceCollectionError):
    """General error for alignment validation failures."""
    pass
