r"""Metadata (:mod:`skbio.metadata`)
================================

.. currentmodule:: skbio.metadata

This module provides classes for storing and working with metadata.

Classes
-------

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

from ._interval import Interval, IntervalMetadata
from ._metadata import SampleMetadata, NumericMetadataColumn, CategoricalMetadataColumn


__all__ = ["Interval", "IntervalMetadata", "Metadata"]
