# for handling optional dependencies, such as polars

import importlib


def _get_package(name):
    msg = f"Using the {name} backend requires the {name} package to be installed."
    try:
        pkg = importlib.import_module(f"{name}")
    except (ImportError, ModuleNotFoundError):
        raise ImportError(msg)

    return pkg
