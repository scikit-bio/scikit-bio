# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn

import numpy as np
import pandas as pd

from ._base import _table_to_numpy
from skbio._config import get_config
from skbio.util import get_package


def _create_table(data, columns=None, index=None, backend=None):
    """Create a table object using the specified backend.

    Parameters
    ----------
    data : table_like
        Input data.
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
        backend = get_config("output")

    if backend == "pandas":
        return pd.DataFrame(data, index=index, columns=columns)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        pl = get_package(backend)
        return pl.DataFrame(data, schema=columns)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


def _create_table_1d(data, index=None, backend=None):
    """Create a 1d array using the specified backend.

    Parameters
    ----------
    data : table_like
        Input data.
    columns : array-like
        Column labels to use if data does not have them.
    index : array-like
        Index labels to use if data does not have them.
    backend : str
        The desired data structure to be used within scikit-bio functions.

    Returns
    -------
    pd.Series or 1-D ndarray
        Representation of the data in the appropriate format depending on the
        underlying configuration option.

    """
    if backend is None:
        backend = get_config("output")

    if backend in ("pandas"):  # , "biom"):
        return pd.Series(data, index=index)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        pl = get_package(backend)
        return pl.Series(values=data)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


def _extract_row_ids(input_data, warn_ids=False):
    """Extract row ids from a dataframe or table."""
    if isinstance(input_data, pd.DataFrame):
        return list(input_data.index)
    # for right now, just going to worry about pandas/polars/numpy,
    # which is to say that if it's not pandas, then it doesn't have ids
    else:
        # Raise warning if sample_ids and feature_ids are both None, as this means
        # that both will have arbitrary integer IDs starting at 0.
        if warn_ids:
            warn(
                (
                    "sample_ids and feature_ids were both None. As a "
                    "result, both have been set to integer IDs "
                    "starting at 0. Namespaces for sample_ids and "
                    "feature_ids are no longer mutually exclusive."
                )
            )
        return list(range(input_data.shape[0]))


def _ingest_table(table, sample_ids=None, feature_ids=None):
    """Process an input data table into individual components.

    Parameters
    ----------
    table : table_like
        The input data table. May be any of the supported formats.
    sample_ids : sequence of str, optional
        IDs corresponding to samples (rows). If ``None``, extraction from input data
        will be attempted. In the case that IDs may not be extracted, they will be
        assigned integer values starting at 0.
    feature_ids : sequence of str, optional
        IDs corresponding to features (columns). If ``None``, extraction from input
        data will be atttempted. In the case that IDs may not be extracted, they will
        be assigned integer values starting at 0.

    Returns
    -------
    data : ndarray of shape (n_samples, n_features)
        The raw numeric values from the input data.
    sample_ids : list of str
        The extracted or provided sample IDs.
    feature_ids : list of str
        The extracted or provided feature IDs.

    """
    # pandas DataFrame
    if isinstance(table, pd.DataFrame):
        data = table.values
        sample_ids = list(table.index) if sample_ids is None else sample_ids
        feature_ids = list(table.columns) if feature_ids is None else feature_ids

    # NumPy array
    elif isinstance(table, np.ndarray):
        data = table

    # BIOM (skbio) Table
    # Check the BIOM-specific attribute "generated_by" before lazy-loading BIOM.
    elif hasattr(table, "generated_by"):
        from skbio.table import Table

        if isinstance(table, Table):
            data, sample_ids, feature_ids = _table_to_numpy(table)

    # Polars DataFrame
    # Can't do an explicit check until polars is imported, so check for schema first.
    elif hasattr(table, "schema"):
        pl = get_package("polars")
        if isinstance(table, pl.DataFrame):
            data = table.to_numpy()
            feature_ids = list(table.schema) if feature_ids is None else feature_ids

    # AnnData object
    elif hasattr(table, "X"):
        adt = get_package("anndata")
        if isinstance(table, adt.AnnData):
            data = table.X
            sample_ids = list(table.obs.index) if sample_ids is None else sample_ids
            feature_ids = list(table.var.index) if feature_ids is None else feature_ids

    else:
        raise TypeError(
            "Input data must be pandas DataFrame, polars DataFrame, numpy ndarray, "
            "skbio.Table, or anndata.AnnData."
        )

    return data, sample_ids, feature_ids
