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


def _warn_renamed(func, oldname, since=None, until=None, param=None):
    """Warn of renamed status."""
    # warning message
    msg_until = f" It will be removed in {until}." if until else ""
    msg_param = f"`{func.__name__}`'s parameter " if param else ""
    msg = (
        msg_param
        + (
            f"`{oldname}` was renamed to `{param or func.__name__}` in {since}."
            " The old name is kept as an alias but is deprecated."
        )
        + msg_until
    )

    # alias of the function
    if not param:
        if not hasattr(func, "warned"):
            simplefilter("once", DeprecationWarning)
            warn(msg, DeprecationWarning)
            func.warned = True

    # alias of a parameter
    else:
        if not hasattr(func, "params_warned"):
            func.params_warned = set()
        if param not in func.params_warned:
            simplefilter("once", DeprecationWarning)
            warn(msg, DeprecationWarning)
            func.params_warned.add(param)
