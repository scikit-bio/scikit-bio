# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, simplefilter


def _warn_once(func, warning, message, param=None):
    r"""Warn once

    Parameters
    ----------
    func : callable
        The function that triggers the warning.
    warning : class
        Warning type.
    message : str
        Warning message.
    param : str, optional
        A parameter of the function that triggers the warning.

    Notes
    -----
    The warning is raised only at the first call of the function. Meanwhile, a "warned"
    attribute is added to the function, which prevents subsequent warnings.

    If a parameter is given, the warning is raised at the first use of the particular
    parameter, while not affect the host function or other parameters. Meanwhile, a
    "warned_params" attribute of the function will contain the parameter to prevent
    subsequent warnings.

    """
    if not param:
        if not hasattr(func, "warned"):
            simplefilter("once", warning)
            warn(message, warning)
            func.warned = True
    else:
        if not hasattr(func, "warned_params"):
            func.warned_params = set()
        if param not in func.warned_params:
            simplefilter("once", warning)
            warn(message, warning)
            func.warned_params.add(param)


def _deprecation_message(name, ver=None, msg=None, append=True):
    """Create a message indicating deprecation."""
    if msg and not append:
        return msg
    if ver:
        message = f"{name} has been deprecated since {ver}."
    else:
        message = f"{name} is deprecated."
    if msg:
        message += " " + msg
    return message


def _warn_deprecated(func, ver=None, msg=None, append=True):
    r"""Warn of deprecated status of a function.

    Parameters
    ----------
    func : callable
        The function that triggers the warning.
    ver : str, optional
        The version when deprecation became effective.
    msg : str, optional
        A custom warning message.
    append : bool, optional
        Append the custom message to the end of the default message (True, default),
        or replace the entire default message with the custom message (False).

    See Also
    --------
    _warn_param_deprecated

    """
    msg = _deprecation_message(f"`{func.__name__}`", ver, msg, append)
    _warn_once(func, DeprecationWarning, msg)


def _warn_param_deprecated(func, param, ver=None, msg=None, append=True):
    r"""Warn of deprecated status of a parameter of a function.

    Parameters
    ----------
    func : callable
        The function that hosts the parameter.
    param : str
        The parameter that triggers the warning.
    ver : str, optional
    msg : str, optional
    append : bool, optional
        Refer to :func:`_warn_deprecated`.

    See Also
    --------
    _warn_deprecated

    """
    msg = _deprecation_message(
        f"`{func.__name__}`'s parameter `{param}`", ver, msg, append
    )
    _warn_once(func, DeprecationWarning, msg, param)


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
        _warn_once(func, DeprecationWarning, msg)

    # alias of a parameter
    else:
        _warn_once(func, DeprecationWarning, msg, param)
