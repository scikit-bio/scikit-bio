# for handling optional dependencies, such as polars

import importlib


def _get_polars():
    """Import polars."""
    msg = "Using the polars backend requires the polars package to be installed."
    try:
        polars = importlib.import_module("polars")
    except (ImportError, ModuleNotFoundError):
        raise ImportError(msg)

    return polars
