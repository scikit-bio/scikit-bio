r"""Configuration Options
=====================

.. currentmodule:: skbio

This module provides a configuration system that controls the global behavior of all
scikit-bio functionalities.

Functions
---------

.. autosummary::
   :toctree: generated/

   get_config
   set_config


Configuration options
---------------------
At present, only one configuration option is available:

**table_output** : *{"pandas", "numpy", "polars"}, default="pandas"*
    The preferred output format of tables. See :ref:`details <table_output>`.

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Any


_SKBIO_OPTIONS = {
    "table_output": "pandas",
}


def set_config(option: str, value: Any):
    """Set a scikit-bio configuration option.

    This function enables users to set the configuration of scikit-bio functions
    globally.

    Parameters
    ----------
    option : str
        The configuration option to be modified.
    value : str
        The value to update the configuration dictionary with.

    Raises
    ------
    ValueError
        If an unknown option is used or if an unsupported value for an option is used.

    Examples
    --------
    >>> from skbio import set_config
    >>> set_config("table_output", "numpy")  # doctest: +SKIP

    """
    if option not in _SKBIO_OPTIONS:
        raise KeyError(f"Unknown option: '{option}'.")

    # Validate option-specific values.
    match option:
        case "table_output":
            pos_opts = ["pandas", "polars", "numpy"]  # , "biom"]
            if value not in pos_opts:
                raise ValueError(f"Unsupported value '{value}' for '{option}'.")

    _SKBIO_OPTIONS[option] = value


def get_config(option: str) -> Any:
    """Get the current value of a scikit-bio configuration option.

    Parameters
    ----------
    option : str
        The configuration option to be found.

    Returns
    -------
    str
        The current value of the configuration option supplied.

    """
    try:
        return _SKBIO_OPTIONS[option]
    except KeyError:
        raise KeyError(f"Unknown option: '{option}'.")
