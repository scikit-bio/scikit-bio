# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from ._base import _table_to_numpy
from skbio.util._array import ingest_array
from skbio._config import get_config
from skbio.util import get_package


def _create_table(data, columns=None, index=None, backend=None):
    """Create a table object using the specified backend.

    Parameters
    ----------
    data : table_like
        Input data.
    columns : array_like
        Column labels to use if data does not have them.
    index : array_like
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
        backend = get_config("table_output")

    if backend == "pandas":
        return pd.DataFrame(data, index=index, columns=columns)
    elif backend == "numpy":
        return np.array(data)
    elif backend == "polars":
        pl = get_package(backend)
        # Polars doesn't support Pandas Index or Series as schema. Convertion needed.
        if columns is not None and not isinstance(columns, Sequence):
            columns = list(columns)
        return pl.DataFrame(data, schema=columns)
    else:
        raise ValueError(f"Unsupported backend: '{backend}'")


def _create_table_1d(data, index=None, backend=None):
    """Create a 1d array using the specified backend.

    Parameters
    ----------
    data : table_like
        Input data.
    columns : array_like
        Column labels to use if data does not have them.
    index : array_like
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
        backend = get_config("table_output")

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
        data will be attempted. In the case that IDs may not be extracted, they will
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
            data = np.asarray(table.X)
            samples = table.obs.index
            features = table.var.index
    # array-like object
    if data is None:
        xp, data = ingest_array(table)
    else:
        xp, data = ingest_array(data)

    # zero-dimensional arrays are considered invalid
    if data.ndim == 0:
        raise TypeError(
            f"'{table.__class__.__name__}' is not a supported table format."
        )

    # convert a 1-D vector into a 2-D array
    if data.ndim == 1 and expand:
        data = xp.reshape(data, (1, -1))

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

    return data, sample_ids, feature_ids


def _aggregate_features(data, aggregator, features=None):
    """Aggregate features in a data table.

    Features assigned the same aggregate ID are summed.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        The raw numeric values from the input data.
    aggregator : callable, mapping, or 1-D array_like
        Rule for aggregating features. Two options:

        1. A callable or mapping (including pandas Series) that maps each feature ID to
           an aggregate ID. Requires `features` supplied.
        2. A 1D array-like that provides one aggregate ID per feature in data order.
           Does not require `features`.

    features : sequence of shape (n_features,), optional
        Feature IDs. If None, no feature IDs are available.

    Returns
    -------
    agg_data : ndarray of shape (n_samples, n_agg_features)
        Aggregated data table.
    agg_features : sequence of shape (n_agg_features,)
        Aggregated feature IDs, in order of first appearance.

    Raises
    ------
    ValueError
        If `aggregator` cannot assign an aggregate ID to every feature.
    TypeError
        If `aggregator` is in an invalid format.

    """
    if features is not None and len(features) != data.shape[1]:
        raise ValueError(
            f"Input table has {data.shape[1]} columns whereas {len(features)} "
            "features were provided."
        )

    # Callable (a function that converts the original ID into an aggregated ID)
    if callable(aggregator):
        if features is None:
            raise ValueError("A callable aggregator requires named features.")
        agg_ids = [aggregator(feature) for feature in features]
        agg_ids = np.asarray(agg_ids, dtype=object)

    # Mapping (e.g., a dictionary of original ID -> aggregated ID)
    elif isinstance(aggregator, (Mapping, pd.Series)):
        if features is None:
            raise ValueError("A mapping aggregator requires named features.")
        try:
            agg_ids = [aggregator[feature] for feature in features]
        except KeyError as e:
            raise ValueError(
                f"Aggregator does not define feature {e.args[0]!r}."
            ) from e
        else:
            agg_ids = np.asarray(agg_ids, dtype=object)

    # One-dimensional array-like (does not require feature IDs)
    else:
        agg_ids = np.asarray(aggregator, dtype=object)
        if agg_ids.ndim == 0:
            raise TypeError("Invalid aggregator format.")
        if agg_ids.ndim != 1 or len(agg_ids) != data.shape[1]:
            raise ValueError(
                "An array-like aggregator must be 1D and have one entry per feature."
            )

    if pd.isna(agg_ids).any():
        raise ValueError("Aggregator must assign every feature an aggregate ID.")

    agg_codes, agg_features = pd.factorize(agg_ids, sort=False)
    agg_data = np.zeros((data.shape[0], len(agg_features)), dtype=data.dtype)
    np.add.at(agg_data, (slice(None), agg_codes), data)

    return agg_data, agg_features
