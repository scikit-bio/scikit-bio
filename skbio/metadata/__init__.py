r"""Metadata (:mod:`skbio.metadata`)
================================

.. currentmodule:: skbio.metadata

This module provides functionality for storing and working with metadata -- the data
that describes other data. While a typical data table (see :mod:`skbio.table`) stores
the measurements of features in samples, metadata provides information about the
samples or features themselves. Examples of metadata include experimental grouping,
demographic and clinical properties of subjects, functional categories of genes and
metabolites, etc.


Sample metadata
---------------

.. autosummary::
   :toctree: generated/

   SampleMetadata


Metadata columns
----------------

.. autosummary::
   :toctree: generated/

   MetadataColumn
   NumericMetadataColumn
   CategoricalMetadataColumn


Interval metadata
-----------------

.. autosummary::
   :toctree: generated/

   Interval
   IntervalMetadata

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._metadata import (
    SampleMetadata,
    MetadataColumn,
    NumericMetadataColumn,
    CategoricalMetadataColumn,
)
from ._interval import Interval, IntervalMetadata

__all__ = [
    "SampleMetadata",
    "MetadataColumn",
    "NumericMetadataColumn",
    "CategoricalMetadataColumn",
    "Interval",
    "IntervalMetadata",
]
