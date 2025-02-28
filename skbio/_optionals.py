# for handling optional dependencies, such as polars

import importlib


def _get_package(name):
    """Attempt to import a package for functions which require said package."""
    msg = f"Using the {name} backend requires the {name} package to be installed."
    try:
        pkg = importlib.import_module(f"{name}")
    except (ImportError, ModuleNotFoundError):
        raise ImportError(msg)

    return pkg


# def _test_package(name):
#     """Check if a package is available for testing purposes."""
#     try:
#         import name
#         import polars.testing as plt
#     except (ImportError, ModuleNotFoundError):
#         has_polars = False
#     else:
#         has_polars = True
