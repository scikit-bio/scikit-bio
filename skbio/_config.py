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
   reset_config


.. _configuration:

Configuration options
---------------------

Settings apply to the current Python process and are not saved between sessions.
Explicit function arguments override the corresponding global setting; passing
None uses that setting. Only functions supporting an option are affected.
See :ref:`compute_engines` for compute engine selection and requirements.

Inspect all current values with ``get_config()``, or one with
``get_config('compute_engine')``. The returned dictionary is a copy.

>>> from skbio import get_config, set_config, reset_config
>>> previous = get_config('compute_engine')
>>> set_config('compute_engine', 'fast')
>>> get_config('compute_engine')
'fast'
>>> reset_config('compute_engine')
>>> get_config('compute_engine')
'cython'
>>> set_config('compute_engine', previous)

Use ``reset_config()`` to restore all defaults. The available options, their
accepted values, and their defaults are listed below.

{option_catalog}

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Any


# The conservative engine. Also what "fast" degrades to when a function offers
# nothing faster, so it cannot depend on the current value of the option.
_DEFAULT_ENGINE = "cython"

# Default, accepted values, and description for each option.
_OPTION_DEFINITIONS = {
    "table_output": (
        "pandas",
        ("pandas", "numpy", "polars"),
        "Preferred table output format. See :ref:`table_output`.",
    ),
    "compute_engine": (
        _DEFAULT_ENGINE,
        ("cython", "numba", "fast"),
        "Default compute engine. See :ref:`compute_engines`.",
    ),
}
_SKBIO_OPTIONS = {key: spec[0] for key, spec in _OPTION_DEFINITIONS.items()}

# Keep the documented catalog in sync with validation and reset defaults.
__doc__ = __doc__.replace(
    "{option_catalog}",
    "\n\n".join(
        f"**{key}** : {{{', '.join(repr(value) for value in values)}}}, "
        f"default={default!r}\n    {description}"
        for key, (default, values, description) in _OPTION_DEFINITIONS.items()
    ),
)


def set_config(option: str, value: Any):
    """Set a scikit-bio configuration option.

    Parameters
    ----------
    option : str
        Option to modify. See :ref:`configuration` for available options.
    value : str
        New value. Explicit function arguments override this global setting.

    Raises
    ------
    KeyError
        If the option is unknown.
    ValueError
        If the value is unsupported for this option.

    See Also
    --------
    get_config
    reset_config

    Notes
    -----
    Settings affect the current Python process only. Optional compute engine
    dependencies are checked when a function uses the engine, not when setting
    the option.

    .. versionchanged:: 0.7.4
        Added ``compute_engine``, accepting 'cython', 'numba', and 'fast'.

    Examples
    --------
    >>> from skbio import get_config, set_config
    >>> previous = get_config('table_output')
    >>> set_config('table_output', 'numpy')
    >>> get_config('table_output')
    'numpy'
    >>> set_config('table_output', previous)

    """
    if option not in _OPTION_DEFINITIONS:
        raise KeyError(f"Unknown option: '{option}'.")
    if value not in _OPTION_DEFINITIONS[option][1]:
        raise ValueError(f"Unsupported value '{value}' for '{option}'.")
    _SKBIO_OPTIONS[option] = value


def get_config(option: str | None = None) -> Any:
    """Get one or all current scikit-bio configuration values.

    Parameters
    ----------
    option : str or None, optional
        Option to inspect. If None (default), return all options.
        See :ref:`configuration` for available options.

        .. versionchanged:: 0.7.4
            Can be omitted to return all options.

    Returns
    -------
    str or dict of str to str
        Current value, or a copy of all current values keyed by option name.
        Changing the returned dictionary does not change configuration.

    Raises
    ------
    KeyError
        If the option is unknown.

    See Also
    --------
    set_config
    reset_config

    Examples
    --------
    >>> from skbio import get_config
    >>> sorted(get_config())
    ['compute_engine', 'table_output']

    """
    if option is None:
        return _SKBIO_OPTIONS.copy()
    try:
        return _SKBIO_OPTIONS[option]
    except KeyError:
        raise KeyError(f"Unknown option: '{option}'.")


def reset_config(option: str | None = None):
    """Restore one or all scikit-bio configuration options to their defaults.

    .. versionadded:: 0.7.4

    Parameters
    ----------
    option : str or None, optional
        Option to reset. If None (default), reset all options.
        See :ref:`configuration` for options and their defaults.

    Raises
    ------
    KeyError
        If the option is unknown. No settings are changed.

    See Also
    --------
    get_config
    set_config

    Examples
    --------
    >>> from skbio import get_config, set_config, reset_config
    >>> previous = get_config('table_output')
    >>> reset_config('table_output')
    >>> get_config('table_output')
    'pandas'
    >>> set_config('table_output', previous)

    """
    if option is None:
        for key, spec in _OPTION_DEFINITIONS.items():
            _SKBIO_OPTIONS[key] = spec[0]
    elif option in _OPTION_DEFINITIONS:
        _SKBIO_OPTIONS[option] = _OPTION_DEFINITIONS[option][0]
    else:
        raise KeyError(f"Unknown option: '{option}'.")


def _resolve_engine(engine, supported, fast=None):
    """Resolve the compute engine for a function call.

    Parameters
    ----------
    engine : str or None
        The engine requested by the caller. If None, the global default
        (``get_config('compute_engine')``) is used.
    supported : tuple of str
        The engines this function supports (e.g. ``('cython', 'numba')``).
    fast : str, optional
        What ``engine='fast'`` resolves to for this function. The caller
        decides, since which engine is fastest depends on the function and on
        what is installed. If not given, ``'fast'`` resolves to the
        conservative default engine, which makes it a no-op for functions that
        have nothing faster to offer.

    Returns
    -------
    str
        The resolved engine name.

    Raises
    ------
    ValueError
        If the resolved engine is not in ``supported``.
    ImportError
        If ``'numba'`` is requested but Numba is not installed.

    """
    if engine is None:
        engine = get_config("compute_engine")
    # Resolved after the global default is read, so a single branch handles
    # "fast" wherever it came from. The fallback is the conservative engine
    # rather than a re-read of the option, so a function that offers nothing
    # faster degrades instead of raising.
    if engine == "fast":
        engine = fast if fast is not None else _DEFAULT_ENGINE
    if engine not in supported:
        raise ValueError(
            f"engine='{engine}' is not supported here; choose from {supported}."
        )
    if engine == "numba":
        try:
            import numba  # noqa: F401
        except ImportError:
            raise ImportError("engine='numba' requires the optional numba dependency.")
    return engine
