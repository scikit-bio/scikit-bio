# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, simplefilter


def _warn_deprecated(func, since, msg=None):
    """Warn of deprecated status."""
    if not hasattr(func, "warned"):
        simplefilter("once", DeprecationWarning)
        message = f"{func.__name__} has been deprecated since {since}."
        if msg:
            message += f" {msg}"
        warn(message, DeprecationWarning)
        func.warned = True


def _warn_renamed(func, oldname, since=None, until=None):
    """Warn of renamed status."""
    if not hasattr(func, "warned"):
        simplefilter("once", DeprecationWarning)
        msg = (
            f"`{oldname}` was renamed to `{func.__name__}` in {since}. The old name "
            "is kept as an alias but is deprecated"
        )
        if until:
            msg += f" and will be removed in {until}."
        msg += "."
        warn(msg, DeprecationWarning)
        func.warned = True
