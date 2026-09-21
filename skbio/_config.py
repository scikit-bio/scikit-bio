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


.. _configuration:

Available options
-----------------

The available options, their accepted values, and their defaults are listed below.

{option_catalog}

.. versionchanged:: 0.7.4
    Added option ``compute_engine``.


How to configure
----------------

Settings apply to the current Python session and are not saved between sessions.
Explicit function parameters override the corresponding global setting. Passing
None uses that setting. Only functions supporting an option are affected. Inspect
all current settings with ``get_config()``, or one with ``get_config(option)``.
The returned dictionary is a copy.

The following example gets and sets the :ref:`compute engine <compute_engines>`.

>>> from skbio import get_config, set_config

>>> option = 'compute_engine'
>>> previous = get_config(option)
>>> set_config(option, 'fast')
>>> get_config(option)
'fast'

>>> set_config(option, previous)
>>> get_config(option)
'cython'

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

# Keep the documented catalog in sync with validation and defaults.
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
        Option to modify. See :ref:`available options <configuration>`.
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
        Option to inspect. See :ref:`available options <configuration>`. If None
        (default), return all options.

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
