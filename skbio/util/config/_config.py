# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

_BACKEND_OPTIONS = {"output": "pandas"}


def set_config(option: str, value: str):
    """Set a scikit-bio config option.

    This function enables users to set the configuration of scikit-bio functions
    globally.

    Parameters
    ----------
    option : str
        The configuration option to be modified. Currently there is only one
        configurable option, ``"output"``.
    value : str
        The value to update the configuration dictionary with. For the
        ``"output"`` option, ``value`` may be set to ``"pandas"``,
        ``"polars"``, or ``"numpy"``. Defaults to ``"pandas"``.

    Raises
    ------
    ValueError
        If an unkown option is used or if an unsupported value for an option is used.

    Examples
    --------
    >>> from skbio import set_config
    >>> set_config("output", "numpy")

    """
    if option not in _BACKEND_OPTIONS:
        raise ValueError(f"Unknown option: '{option}'")
    # possible options for now
    pos_opts = ["pandas", "polars", "numpy"]  # , "biom"]
    if value not in pos_opts:
        raise ValueError(f"Unsupported value '{value}' for '{option}'")
    _BACKEND_OPTIONS[option] = value


def get_config(option: str) -> str:
    """Get the current value of an skbio config option.

    Parameters
    ----------
    option : str
        The configuration option to be found.

    Returns
    -------
    str
        The current value of the configuration option supplied.

    """
    if option not in _BACKEND_OPTIONS:
        raise ValueError(f"Unknown option: '{option}'")
    return _BACKEND_OPTIONS[option]
