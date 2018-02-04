# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class TestingUtilError(Exception):
    """Raised when an exception is needed to test exception handling."""
    pass


class OverrideError(AssertionError):
    """Raised when a property does not exist in the parent class."""
    pass
