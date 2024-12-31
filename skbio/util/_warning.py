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
        Function that triggers the warning.
    warning : class
        Warning type.
    message : str
        Warning message.
    param : str, optional
        Parameter of the function that triggers the warning.

    Notes
    -----
    The warning is raised only at the first call of the function. Meanwhile, a "_warned"
    attribute is added to the function, which prevents subsequent warnings.

    If a parameter is given, the warning is raised at the first use of the particular
    parameter, while not affect the host function or other parameters. Meanwhile, a
    "_warned_params" attribute of the function will contain the parameter to prevent
    subsequent warnings.

    """
    if not param:
        if not hasattr(func, "_warned"):
            simplefilter("once", warning)
            warn(message, warning)
            func._warned = True
    else:
        if not hasattr(func, "_warned_params"):
            func._warned_params = set()
        if param not in func._warned_params:
            simplefilter("once", warning)
            warn(message, warning)
            func._warned_params.add(param)


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
        Function that triggers the warning.
    ver : str, optional
        Version when deprecation became effective.
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
        Function that hosts the parameter.
    param : str
        Parameter that triggers the warning.
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


def _renaming_message(oldname, newname, ver=None):
    """Create a message indicating renaming."""
    if ver:
        message = f"{oldname} was renamed to {newname} in {ver}."
    else:
        message = f"{oldname} has been renamed to {newname}."
    message += " The old name is kept as an alias but is deprecated."
    return message


def _warn_renamed(func, oldname, ver=None):
    """Warn of renaming status of a function."""
    msg = _renaming_message(f"`{oldname}`", f"`{func.__name__}`", ver)
    _warn_once(func, DeprecationWarning, msg)


def _warn_param_renamed(func, param, oldname, ver=None):
    """Warn of renaming status of a parameter of a function."""
    msg = _renaming_message(
        f"`{func.__name__}`'s parameter `{oldname}`", f"`{param}`", ver
    )
    _warn_once(func, DeprecationWarning, msg, param)
