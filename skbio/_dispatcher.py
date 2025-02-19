import numpy as np
import pandas as pd

from skbio.table import Table
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
        return polars.DataFrame(data, schema=columns)
    # elif backend == "biom":
    #     return Table(data, observation_ids=index, sample_ids=columns)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


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

    if backend in ("pandas"):  # , "biom"):
        return pd.Series(data, index=index)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        polars = _get_polars()
        return polars.Series(values=data)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


def ingest_array(input_data, row_ids=None, col_ids=None, dtype=None):
    # pd.DataFrame
    if isinstance(input_data, pd.DataFrame):
        data_ = input_data.values
        row_ids = list(input_data.index) if row_ids is None else row_ids
        col_ids = list(input_data.columns) if col_ids is None else col_ids
    # numpy array
    elif isinstance(input_data, np.ndarray):
        data_ = input_data
    # BIOM (skbio) Table
    elif isinstance(input_data, Table):
        # BIOM puts samples as columns and features as rows, so need to handle
        # accordingly
        # Or, maybe it's better just to enforce that rows must be samples and
        # columns must be features (observations), then we don't worry about it
        print(
            "BIOM format uses samples as columns and features as rows. Most "
            "scikit-bio functions expect samples as rows and features as columns."
            "Please ensure your input is in the correct orientation.\n"
        )
        data_ = input_data.to_dataframe().values
        row_ids = (
            list(input_data.ids(axis="observation")) if row_ids is None else row_ids
        )
        col_ids = list(input_data.ids()) if col_ids is None else col_ids
    # pl.DataFrame
    elif hasattr(input_data, "schema"):
        pl = _get_polars()
        if isinstance(input_data, pl.DataFrame):
            data_ = input_data.to_numpy()
            col_ids = list(input_data.schema) if col_ids is None else col_ids
    # ndarray

    else:
        raise TypeError(
            "Input data must be pandas DataFrame, polars DataFrame, or numpy ndarray"
        )

    return data_, row_ids, col_ids
