_BACKEND_OPTIONS = {"output": "pandas"}


def set_config(option: str, value: str):
    """Set an skbio config option.

    Parameters
    ----------
    option : str
        The configuration option to be modified.
    value : str
        The value to update the configuration dictionary with.

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
