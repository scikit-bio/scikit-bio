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

**table_output** : *{"pandas", "numpy", "polars"}, default="pandas"*
    The preferred output format of tables. See :ref:`details <table_output>`.

**engine** : *{"cython", "numba"}, default="cython"*
    The default compute engine for functions that support multiple engines.
    Functions accept an ``engine`` argument to override this per call. The
    ``"numba"`` engine requires the optional Numba dependency.

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
    "engine": "cython",
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
        case "engine":
            pos_opts = ["cython", "numba"]
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


def _resolve_engine(engine, supported):
    """Resolve the compute engine for a function call.

    Parameters
    ----------
    engine : str or None
        The engine requested by the caller. If None, the global default
        (``get_config("engine")``) is used.
    supported : tuple of str
        The engines this function supports (e.g. ``("cython", "numba")``).

    Returns
    -------
    str
        The resolved engine name.

    Raises
    ------
    ValueError
        If the resolved engine is not in ``supported``.
    ImportError
        If ``"numba"`` is requested but Numba is not installed.

    """
    if engine is None:
        engine = get_config("engine")
    if engine not in supported:
        raise ValueError(
            f"engine='{engine}' is not supported here; choose from {supported}."
        )
    if engine == "numba":
        try:
            import numba  # noqa: F401
        except ImportError:
            raise ImportError(
                "engine='numba' requires the optional numba dependency."
            )
    return engine
