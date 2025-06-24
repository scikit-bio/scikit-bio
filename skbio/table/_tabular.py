# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
from collections.abc import Sequence

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


def _ingest_table(table, sample_ids=None, feature_ids=None, expand=True):
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
    expand : bool, optional
        If table has only one dimension and this parameter is True, expand the data to
        two dimensions (i.e., a single-row matrix). Otherwise, raise an error.

    Returns
    -------
    data : ndarray of shape (n_samples, n_features)
        The raw numeric values from the input data.
    sample_ids : list of str
        The extracted or provided sample IDs.
    feature_ids : list of str
        The extracted or provided feature IDs.

    Raises
    ------
    TypeError
        If input table format is not supported.
    ValueError
        If number of provided sample/feature IDs doesn't match table dimensions.

    """
    data, samples, features = None, None, None

    # Python (nested) list, tuple, etc.
    if isinstance(table, Sequence) and not isinstance(table, (str, bytes)):
        data = np.asarray(table)

    # NumPy array
    # to be replaced with `aac.is_array_api_obj(table)`
    elif isinstance(table, np.ndarray):
        data = table

    # pandas DataFrame
    elif isinstance(table, pd.DataFrame):
        data = table.to_numpy()
        samples = table.index
        features = table.columns

    # BIOM (skbio) Table
    # Check the BIOM-specific attribute "generated_by" before lazy-loading BIOM.
    elif hasattr(table, "generated_by"):
        from skbio.table import Table

        if isinstance(table, Table):
            data, samples, features = _table_to_numpy(table)

    # Polars DataFrame
    # Can't do an explicit check until polars is imported, so check for schema first.
    elif hasattr(table, "schema"):
        pl = get_package("polars")
        if isinstance(table, pl.DataFrame):
            data = table.to_numpy()
            features = table.schema

    # AnnData object
    elif hasattr(table, "X"):
        adt = get_package("anndata")
        if isinstance(table, adt.AnnData):
            data = table.X
            samples = table.obs.index
            features = table.var.index

    if data is None:
        raise TypeError("Input table format is not supported.")
    if data.ndim == 1 and expand:
        data = data.reshape(1, -1)
    if data.ndim < 2:
        raise ValueError("Input table has less than 2 dimensions.")

    lenerr = "Input table has {0} {1}s whereas {2} {1} IDs were provided."
    if sample_ids is None:
        sample_ids = samples
    elif len(sample_ids) != data.shape[0]:
        raise ValueError(lenerr.format(data.shape[0], "sample", len(sample_ids)))
    if feature_ids is None:
        feature_ids = features
    elif len(feature_ids) != data.shape[1]:
        raise ValueError(lenerr.format(data.shape[1], "feature", len(feature_ids)))

    # cast to lists
    if sample_ids is not None and not isinstance(sample_ids, list):
        sample_ids = list(sample_ids)
    if feature_ids is not None and not isinstance(feature_ids, list):
        feature_ids = list(feature_ids)

    return data, sample_ids, feature_ids


def _ingest_vector(vector, ids=None):
    """Process an input vector into individual components.

    Parameters
    ----------
    vector : vector_like
        The input vector. May be any of the supported formats.
    ids : sequence of str, optional
        Positional IDs of the vector. If ``None``, extraction from input vector will
        will be attempted.

    Returns
    -------
    data : ndarray of shape (n,)
        The raw numeric values from the input data.
    ids : list of str
        The extracted or provided IDs.

    Raises
    ------
    TypeError
        If input table format is not supported.
    ValueError
        If number of provided IDs doesn't match vector length.

    See Also
    --------
    _ingest_table

    """
    lenerr = "Vector has a length of {} but {} IDs are provided."
    data, index = None, None

    # Python (nested) list, tuple, etc.
    if isinstance(vector, Sequence) and not isinstance(vector, (str, bytes)):
        data = np.asarray(vector)

    # NumPy array
    # to be replaced with `aac.is_array_api_obj(table)`
    elif isinstance(vector, np.ndarray):
        data = vector

    # pandas Series
    elif isinstance(vector, pd.Series):
        data = vector.to_numpy()
        if ids is None:
            index = vector.index
    else:
        raise TypeError("Input vector format is not supported.")
    if data.ndim < 1:
        raise ValueError("Input vector has less than 1 dimension.")

    if ids is None:
        if index is not None:
            ids = index.tolist()
    else:
        if not isinstance(ids, list):
            ids = list(ids)
        if len(ids) != data.shape[0]:
            raise ValueError(lenerr.format(data.shape[0], len(ids)))

    return data, ids
