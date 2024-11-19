# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, simplefilter


class SkbioWarning(Warning):
    """Filter our warnings from warnings given by 3rd parties."""

    pass


class EfficiencyWarning(SkbioWarning):
    """Warn about potentially accidental use of inefficient code.

    For example, if a user doesn't have an optimized version of a
    function/algorithm available in their scikit-bio installation, a slower,
    pure-Python implementation may be used instead. This warning can be used to
    let the user know they are using a version of the function that could be
    potentially orders of magnitude slower.

    """

    pass


class DeprecationWarning(DeprecationWarning, SkbioWarning):
    """Used to indicate deprecated functionality in scikit-bio."""

    pass


def _warn_deprecated(func, ver, msg=None):
    """Warn of deprecated status."""
    if not hasattr(func, "warned"):
        simplefilter("once", DeprecationWarning)
        message = f"{func.__name__} is deprecated as of {ver}."
        if msg:
            message += f" {msg}"
        warn(message, DeprecationWarning)
        func.warned = True
