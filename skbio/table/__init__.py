r"""Data Table (:mod:`skbio.table`)
==================================

.. currentmodule:: skbio.table

Tabular data are among the most common data structures in bioinformatics and many other
fields. scikit-bio provides a powerful table class for biological data (:class:`Table`),
a flexible system that accepts multiple input formats and returns results in the user's
preferred format, and algorithms for manipulating tabular data in biological research.

See an |table_intro|_ to the specifics, nomenclature and usage of data tables in
scikit-bio.


BIOM table
----------

`Biological Observation Matrix (BIOM) <https://biom-format.org/>`_ is an efficient and
versatile sparse table format designed for biological "omic" data types. It is the
native table format in scikit-bio.

.. autosummary::
   :toctree: generated/

   Table
   example_table


Table-like formats
------------------

scikit-bio functions directly operate on various "table-like" formats, such as BIOM
table, Pandas and Polars dataframes, NumPy array and AnnData objects, without the need
for explicit format conversion.

- :ref:`table_like`
- :ref:`array_like`


Data augmentation
-----------------

Techniques for creating synthetic samples based on the current data and biological
properties. Helpful for improving the accuracy and robustness of machine learning
models.

.. autosummary::
   :toctree: generated/

   phylomix
   compos_cutmix
   aitchison_mixup
   mixup


.. |table_intro| replace:: **Introduction**
.. _table_intro:

Introduction
------------

In scikit-bio, a table is an **ID-labeled data matrix**.

Data are stored in a two-dimensional structure with **rows = samples** and
**columns = features**. Cell values represent counts, abundances, or other numeric
measurements of a given feature within a given sample. The entire table must be in
a homogeneous data type.

Internally, a table is decomposed into a 2-D data matrix, and lists of sample IDs and
feature IDs (:ref:`explained below <samples_and_features>`).

IDs are optional. When they are absent, scikit-bio uses integer indices to refer to
samples and features. When IDs are present, however, scikit-bio uses them to bridge
structured biological knowledge with numeric data, enabling convenient ID-based in
addition to index-based operations (e.g., :func:`~skbio.diversity.alpha_diversity` and
:class:`~skbio.stats.distance.DistanceMatrix`).


.. _table_like:

Supported table-like formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scikit-bio supports a variety of table-like formats used in bioinformatics, scientific
computing and data science in general. When a function parameter is annotated as
``table_like``, it accepts any of the following formats:

+--------------------------------------+---------------------+------------------------+
| Format                               | Sample IDs          | Feature IDs            |
+======================================+=====================+========================+
| skbio :class:`~skbio.table.Table`\   | "sample"            | "observation"          |
| (i.e., BIOM table)                   |                     |                        |
+--------------------------------------+---------------------+------------------------+
| anndata :class:`~anndata.AnnData`    | "obs" (observation) | "var" (variable)       |
+--------------------------------------+---------------------+------------------------+
| Pandas :class:`~pandas.DataFrame`    | "index" (rows)      | "columns"              |
+--------------------------------------+---------------------+------------------------+
| Polars |pl_df|_                      | N/A (rows)          | "schema"               |
+--------------------------------------+---------------------+------------------------+
| array-like objects \                 | N/A                 | N/A                    |
| (:ref:`detailed below <array_like>`) |                     |                        |
+--------------------------------------+---------------------+------------------------+

Input handling of supported types is automatic. No manual conversion is required. If
provided, sample and feature IDs will be preserved and propagated to downstream
results. Otherwise, integer indices starting at 0 will be used.


.. _array_like:

Supported array-like formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scikit-bio supports a wide range of array libraries in the Python ecosystem. A
parameter annotated as ``array_like`` accepts any of the following formats:

- Numpy :class:`~numpy.ndarray`, the "native" format which most of scikit-bio's
  functions are optimized for (e.g., utilizing vectorization if possible).
- Any array object compliant with the `Python array API standard
  <https://data-apis.org/array-api/latest/>`_. Examples are PyTorch ``Tensor``,
  CuPy ``ndarray``, JAX ``Array``, Dask ``Array``, and sparse ``SparseArray``.
- Any object convertible to a ≥1-dimensional NumPy array via :func:`~numpy.asarray`.
  Examples are plain Python lists and tuples.

.. note::
   A proportion of scikit-bio functions offer *native* support for alternative array
   libraries. For instance, log-ratio transformations like
   :func:`~skbio.stats.composition.clr` can consume and return GPU-resident
   tensors with arbitrary dimensions, eliminating the overhead of round-tripping
   through NumPy and significantly accelerating computation.


.. _samples_and_features:

Samples and features
^^^^^^^^^^^^^^^^^^^^

Typically, the two table dimensions represent:

- **Sample IDs** (row labels): biological samples, or any other experimental units
  (specimens, subjects, sites, time points, cells, etc.).
- **Feature IDs** (column labels): variables or characteristics measured in each
  sample (taxa, genes, molecules, environmental factors, etc.).

Terminology varies across disciplines and can sometimes be confusing. Scikit-bio
standardizes on "sample" and "feature" -- the most common terms in data science.
But you may encounter aliases:

* In the BIOM format (wrapped by :mod:`skbio.table`), features are called
  "observations" and samples stay "samples". (Beware: in anndata, "observations"
  instead refer to samples.)
* Some other scikit-bio sub-modules adopt field-specific terms, documented on their
  index pages, such as:

  - :mod:`skbio.diversity`: sample → "community", feature → "taxon"
  - :mod:`skbio.stats.composition`: sample → "composition", feature → "component"

Depending on the research task, you are free to transpose a table to swap samples and
features.


.. _table_params:

Common parameters
^^^^^^^^^^^^^^^^^

Many functions involving tabular data share a set of common parameters that control how
IDs are handled and specify output format preferences:

``sample_ids`` : *list of str, optional*
    Identifiers for samples (rows). If not provided implicitly by the input table or
    explicitly by the user, defaults to ``range(n_samples)``. This parameter is useful
    when the input format doesn't support row labels (e.g., NumPy arrays) or when you
    want to override existing labels.

``feature_ids`` : *list of str, optional*
    Identifiers for features (columns). Analogous to ``sample_ids``.  Default:
    ``range(n_features)``.

``output_format`` : *{"pandas", "numpy", "polars"}*, optional
    Preferred output format **for this call only** (default: ``"pandas"``).  This
    setting overrides the global configuration (see below).


.. _table_output:

Output formats
^^^^^^^^^^^^^^

Some functions that *produce* tables can return the result in one of three formats:

- Pandas :class:`~pandas.DataFrame` and :class:`~pandas.Series` (default)
- Polars |pl_df|_ and |pl_s|_
- NumPy :class:`~numpy.ndarray` (2-D and 1-D)

There are two ways to control the output format:

1. Set the desired output format on a per-function basis, using the ``output_format``
   parameter.

.. code-block:: python

    from skbio.stats.ordination import cca

    # This specific call will return an `OrdinationResults` object whose attributes are
    # NumPy arrays.
    res = cca(Y, X, output_format="numpy")

2. Set the ``table_output`` configuration option using :func:`skbio.set_config`. It
   will change the global behavior of all scikit-bio functions.

.. code-block:: python

    # `set_config` is available as a top level import from skbio
    from skbio import set_config

    # Set output format to NumPy arrays, or
    set_config("table_output", "numpy")

    # Return to default Pandas output
    set_config("table_output", "pandas")


..
   Polars DataFrame and Series cannot be cross-linked therefore the following
   workaround is used. See: https://github.com/pola-rs/polars/issues/7027

.. |pl_df| replace:: ``DataFrame``
.. _pl_df: https://docs.pola.rs/api/python/stable/reference/dataframe/index.html

.. |pl_s| replace:: ``Series``
.. _pl_s: https://docs.pola.rs/api/python/stable/reference/series/index.html


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.table._base import Table, example_table
from skbio.table._augment import phylomix, compos_cutmix, aitchison_mixup, mixup

__all__ = [
    "Table",
    "example_table",
    "phylomix",
    "compos_cutmix",
    "aitchison_mixup",
    "mixup",
]
