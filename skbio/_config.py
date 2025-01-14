_BACKEND_OPTIONS = {"table_backend": "pandas"}


def set_option(option_name: str, value: str):
    """Set an skbio config option."""
    if option_name not in _BACKEND_OPTIONS:
        raise ValueError(f"Unkown option: {option_name}")
    _BACKEND_OPTIONS[option_name] = value


def get_option(option_name: str) -> str:
    """Get the current value of an skbio config option."""
    if option_name not in _BACKEND_OPTIONS:
        raise ValueError(f"Unknown option: {option_name}")
    return _BACKEND_OPTIONS[option_name]
