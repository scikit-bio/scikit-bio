import numpy as np
import pandas as pd

from skbio.table import Table
from ._config import get_option
from ._optionals import _get_package


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
        pl = _get_package(backend)
        return pl.DataFrame(data, schema=columns)
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
        pl = _get_package(backend)
        return pl.Series(values=data)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


def ingest_array(input_data, row_ids=None, col_ids=None, dtype=None):
    """Process an input dataframe, table, or array into individual components.

    Parameters
    ----------
    input_data : tabular
        The original source of data. May be pandas or polars DataFrame, numpy array,
        or BIOM table.
    row_ids : list of str
        IDs corresponding to the rows of the input data. If ``None``, extraction from
        input data will be attempted. In the case that IDs may not be extracted, they
        will be assigned integer values starting at 0.
    col_ids : list of str
        IDs corresponding to the columns of the input data. If ``None``, extraction
        from input data will be atttempted. In the case that IDs may not be extracted,
        they will be assigned integer values starting at 0.

    Returns
    -------
    data_ : ndarray
        The raw numeric values from the input data.
    row_ids : list of str
        The extracted or provided row_ids.
    col_ids : list of str
        The extracted or provided col_ids.

    """
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
        # Can't do an explicit check until polars is imported,
        # so check for schema first
        pl = _get_package("polars")
        if isinstance(input_data, pl.DataFrame):
            data_ = input_data.to_numpy()
            col_ids = list(input_data.schema) if col_ids is None else col_ids
    # anndata

    else:
        raise TypeError(
            "Input data must be pandas DataFrame, polars DataFrame, or numpy ndarray"
        )

    return data_, row_ids, col_ids


def extract_row_ids(input_data):
    """Extract row ids from a dataframe or table."""
    if isinstance(input_data, pd.DataFrame):
        return list(input_data.index)
    # for right now, just going to worry about pandas/polars/numpy,
    # which is to say that if it's not pandas, then it doesn't have ids
    else:
        return list(range(input_data.shape[0]))
