# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class FormatIdentificationWarning(Warning):
    """Warn when the sniffer of a format cannot confirm the format."""
    pass


class ArgumentOverrideWarning(Warning):
    """Warn when a user provided kwarg differs from a guessed kwarg."""
    pass
