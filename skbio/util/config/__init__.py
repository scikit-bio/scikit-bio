r"""Configuration Options (:mod:`skbio.util.config`)
================================================

.. currentmodule:: skbio.util.config

This module provides a flexible configuration system which allows scikit-bio
functions to accept multiple types of input data structures and return results in the
users preferred format.

The dispatch system follows a design similar to scikit-learn's ``set_output`` API,
allowing users to work with different data formats without changing existing workflows.

Functions
---------

.. autosummary::
   :toctree: generated/

   get_config
   set_config


Supported Input Formats
-----------------------
- Numpy :class:`~numpy.ndarray` (2D)
- pandas :class:`~pandas.DataFrame`
- Polars :class:`~polars.DataFrame`
- scikit-bio :class:`~skbio.table.Table`
- :class:`~anndata.AnnData` object

Input handling of supported types is handled automatically. No need to set a
configuration variable to use whichever input type you like. When possible, row and
column identifiers will be preserved and carried through operations to the results.
When IDs are not available from the input data and not provided explicitly, integer
indices starting from 0 will be used.

Supported Output Formats
------------------------
- Numpy :class:`~numpy.ndarray` (2D and 1D)
- pandas :class:`~pandas.DataFrame` and :class:`~pandas.Series` (default)
- Polars :class:`~polars.DataFrame` and :class:`~polars.Series`

Configuring Output Format
-------------------------
There are two ways to control the output format.

The first option is to use the `set_config` function. This function will change the
global behavior of scikit-bio functions.:

.. code-block:: python

    # set_config is available as a top level import from skbio
    from skbio import set_config

    # Set output format to NumPy arrays
    set_config("output", "numpy")

    # Return to default pandas output
    set_config("output", "pandas")

The second option is to set the desired output format on a per-function basis, using
the `output_format` parameter:

.. code-block:: python

    from skbio.stats.ordination import cca

    # This specific call will return an
    # :class:`~skbio.stats.ordination.OrdinationResults` object whose attributes are
    # numpy arrays
    res = cca(Y, X, output_format="numpy")


Notes
-----
When using the default pandas backend, existing workflows remain unchanged.
Additionally, input and output formats do not have to match.

"""

from ._config import get_config, set_config
from ._dispatcher import (
    _create_table,
    _create_table_1d,
    _extract_row_ids,
    _ingest_array,
)
from ._optionals import _get_package

__all__ = [
    "get_config",
    "set_config",
    "_create_table",
    "_create_table_1d",
    "_extract_row_ids",
    "_ingest_array",
    "_get_package",
]
