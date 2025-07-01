.. meta::
   :description: Tabular data formats supported by scikit-bio

Table-like formats
==================

This module provides a flexible configuration system which allows scikit-bio
functions to accept multiple types of input data structures and return results in the
users preferred format.

The dispatch system follows a design similar to scikit-learn's ``set_output`` API,
allowing users to work with different data formats without changing existing workflows.


The TableLike Type
------------------
The :data:`~skbio.util.TableLike` type is the union of the following types of objects:

- Numpy :class:`~numpy.ndarray` (2D)
- pandas :class:`~pandas.DataFrame`
- Polars :class:`~polars.DataFrame`
- scikit-bio :class:`~skbio.table.Table`
- :class:`~anndata.AnnData` object

For all functions which accept a ``TableLike`` type as input, any objects from this
list are acceptable input. Input handling of supported types is handled automatically.
No need to set a configuration variable to use whichever input type you like. When
possible, row and column identifiers will be preserved and carried through operations
to the results. When IDs are not available from the input data and not provided
explicitly, integer indices starting from 0 will be used.


Sample and Feature Identifiers
------------------------------
In scikit-bio, data is commonly organized in a two-dimensional structure:

- **sample_ids**: Identifiers for rows, typically representing biological samples,
  observations, or experimental units (e.g., patients, sites, time points).

- **feature_ids**: Identifiers for columns, typically representing measured variables
  or characteristics (e.g., taxa, genes, OTUs, metabolites, environmental parameters).

Different data formats use different terminology for these concepts:

- ``pandas``: "index" (rows) and "columns"
- ``anndata``: "obs" (samples) and "var" (features)
- ``scikit-bio Table``: "sample" (samples) and "observation"
  (features)
- ``polars``: rows (by position) and "schema" (features)


Common TableLike Parameters
---------------------------
Many functions that accept ``TableLike`` inputs share a set of common parameters that
control how identifiers are handled and specify output format preferences:

sample_ids : list of str, optional
    List of identifiers for samples (rows). If not provided implicitly by the input
    data structure or explicitly by the user, defaults to integers starting at zero.
    This parameter is useful when the input format doesn't support row labels
    (e.g., NumPy arrays) or when you want to override existing labels.

feature_ids : list of str, optional
    List of identifiers for features (columns). If not provided implicitly by the
    input data structure or explicitly by the user, defaults to integers starting
    at zero. This parameter is useful when the input format doesn't support column
    labels (e.g., NumPy arrays) or when you want to override existing labels.

output_format : str, optional
    Specifies the desired format for the output. Valid options are:

    - ``"pandas"``: Return pandas DataFrames/Series (default)
    - ``"numpy"``: Return NumPy ndarrays
    - ``"polars"``: Return Polars DataFrames/Series

    Note that this parameter overrides the global configuration setting for
    the specific function call.


Supported Output Formats
------------------------
Currently, scikit-bio functions support outputing the following types of data. Note
that this list is not equivalent to the input types supported through the ``TableLike``
type.

- Numpy :class:`~numpy.ndarray` (2D and 1D)
- pandas :class:`~pandas.DataFrame` and :class:`~pandas.Series` (default)
- Polars :class:`~polars.DataFrame` and :class:`~polars.Series`


Configuring Output Format
-------------------------
There are two ways to control the output format.

The first option is to use the :func:`set_config` function. This function will change
the global behavior of scikit-bio functions.

.. code-block:: python

    # set_config is available as a top level import from skbio
    from skbio import set_config

    # Set output format to NumPy arrays
    set_config("output", "numpy")

    # Return to default pandas output
    set_config("output", "pandas")

The second option is to set the desired output format on a per-function basis, using
the ``output_format`` parameter.

.. code-block:: python

    from skbio.stats.ordination import cca

    # This specific call will return an
    # :class:`~skbio.stats.ordination.OrdinationResults` object whose attributes are
    # numpy arrays
    res = cca(Y, X, output_format="numpy")
