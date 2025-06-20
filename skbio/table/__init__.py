r"""Data Table (:mod:`skbio.table`)
==================================

.. currentmodule:: skbio.table

This module provides support for interaction with data tables.


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
for explicit format conversion. Read below on the specifics, nomenclature and usage of
supported table-like formats.

.. autosummary::
   :toctree: generated/

   table_like


Data augmentation
-----------------

Techniques for creating synthetic samples based on the current data and biological
properties. Helpful for improving the accuracy and robustness of machine learning
models.

.. autosummary::
   :toctree: generated/

   phylomix
   compositional_cutmix
   aitchison_mixup
   mixup

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.table._base import Table, example_table
from skbio.table._augment import phylomix, compositional_cutmix, aitchison_mixup, mixup
from skbio.table import _dispatcher as table_like

__all__ = [
    "Table",
    "example_table",
    "phylomix",
    "compositional_cutmix",
    "aitchison_mixup",
    "mixup",
    "table_like",  # only for documentation
]
