# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, simplefilter


def _warn_deprecated(func, ver, msg=None):
    """Warn of deprecated status."""
    if not hasattr(func, "warned"):
        simplefilter("once", DeprecationWarning)
        message = f"{func.__name__} is deprecated as of {ver}."
        if msg:
            message += f" {msg}"
        warn(message, DeprecationWarning)
        func.warned = True


def _warn_renamed(func, oldname, ver=None):
    """Warn of renaming status."""
    if not hasattr(func, "warned"):
        simplefilter("once", DeprecationWarning)
        if not ver:
            ver = "a future release"
        msg = (
            f"`{oldname}` has been renamed as `{func.__name__}`. The old name is "
            f"deprecated and will be removed in {ver}."
        )
        warn(msg, DeprecationWarning)
        func.warned = True
