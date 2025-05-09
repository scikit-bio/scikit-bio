# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import importlib


def _get_package(name, raise_error=True, no_bool=True):
    """Attempt to import a package for functions which require said package."""
    msg = f"Using the {name} backend requires the {name} package to be installed."
    try:
        pkg = importlib.import_module(f"{name}")
        if no_bool:
            return pkg
        else:
            return pkg, True
    except (ImportError, ModuleNotFoundError):
        if raise_error:
            raise ImportError(msg)
        else:
            return None, False
