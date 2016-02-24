# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import TestRunner

from ._feature import Feature
from ._interval import IntervalMetadata
from .intersection import Interval
from .mixin import IntervalMetadataMixin, PositionalMetadataMixin, MetadataMixin

__all__ = ['Feature', 'Interval', 'IntervalMetadata',
           'MetadataMixin', 'PositionalMetadataMixin',
           'IntervalMetadataMixin']

test = TestRunner(__file__).test
