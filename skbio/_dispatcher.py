import numpy as np
import pandas as pd

from ._config import get_option
from ._optionals import _get_polars


def create_table(data, columns=None, index=None, backend=None):
    """Create a table object using the specified backend.

    Parameters
    ----------
    data : ndarray
    columns : array-like
        Column labels to use if data does not have them.
    index : array-like
        Index labels to use if data does not have them.
    backend : str
        The desired data structure to be used within scikit-bio functions.

    Returns
    -------
    pd.DataFrame or np.array
        Representation of the data in the appropriate format depending on the
        underlying configuration option.

    """
    if backend is None:
        backend = get_option("tabular_backend")

    if backend == "pandas":
        return pd.DataFrame(data, index=index, columns=columns)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        polars = _get_polars()
        return polars.LazyFrame(data)
    else:
        raise ValueError(f"Unsupported backend '{backend}'")


def create_table_1d(data, index=None, backend=None):
    """Create a 1d array using the specified backend.

    Parameters
    ----------
    data : ndarray
    columns : array-like
        Column labels to use if data does not have them.
    index : array-like
        Index labels to use if data does not have them.
    backend : str
        The desired data structure to be used within scikit-bio functions.

    Returns
    -------
    pd.Series or numpy array
        Representation of the data in the appropriate format depending on the
        underlying configuration option.

    """
    if backend is None:
        backend = get_option("tabular_backend")

    if backend == "pandas":
        return pd.Series(data, index=index)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        polars = _get_polars()
        return polars.Series(values=data)
    else:
        raise ValueError(f"Unsupported backend '{backend}'")
