r"""Data Table (:mod:`skbio.table`)
==================================

.. currentmodule:: skbio.table

This module provides support for interaction with data tables.


BIOM table
----------

Biological Observation Matrix (BIOM) is an efficient and versatile table format
designed for biological "omic" data types.

.. autosummary::
   :toctree: generated/

   Table
   example_table

Please refer to the `BIOM documentation <https://biom-format.org/>`_ for the
instructions on working with BIOM tables.


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

__all__ = [
    "Table",
    "example_table",
    "phylomix",
    "compositional_cutmix",
    "aitchison_mixup",
    "mixup",
]
