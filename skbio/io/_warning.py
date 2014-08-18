from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class UnprovenFormatWarning(Warning):
    """Warn when the identifer of a format cannot confirm expected value."""
    pass


class ArgumentOverrideWarning(Warning):
    """Warn when a user provided kwarg differs from a guessed kwarg."""
    pass
