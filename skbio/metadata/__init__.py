r"""
Metadata (:mod:`skbio.metadata`)
================================

.. currentmodule:: skbio.metadata

This module provides classes for storing and working with metadata.
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util import TestRunner


from ._interval import IntervalMetadata, Interval
# from ._mixin import (IntervalMetadataMixin, PositionalMetadataMixin,
#                      MetadataMixin)

# __all__ = ['Interval', 'IntervalMetadata',
#            'MetadataMixin', 'PositionalMetadataMixin',
#            'IntervalMetadataMixin']

__all__ = ['Interval', 'IntervalMetadata',]

test = TestRunner(__file__).test
