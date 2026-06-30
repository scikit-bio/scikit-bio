# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Optional dependencies."""

import importlib


def get_package(name, raise_error=True):
    r"""Attempt to import a package for functions that require the said package.

    Parameters
    ----------
    name : str
        Name of the package.
    raise_error : bool, optional
        If the package cannot be imported, raise an ``ImportError`` if True (default),
        otherwise return None.

    Returns
    -------
    module
        Imported package (None if not available).

    Raises
    ------
    ImportError
        If the package cannot be imported.

    """
    msg = (
        f'Optional dependency "{name}" is not found. '
        "Install it to use relevant functionalities."
    )
    try:
        return importlib.import_module(name)
    except (ImportError, ModuleNotFoundError):
        if raise_error:
            raise ImportError(msg)


def _is_usable(mod, probe):
    """Return True if `mod` is importable AND `probe()` runs without error."""
    if mod is None:
        return False
    try:
        probe()
        return True
    except Exception:
        return False
