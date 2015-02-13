"""
Statistics (:mod:`skbio.stats`)
===============================

.. currentmodule:: skbio.stats

This package contains various statistical methods, including ordination
techniques and distance matrix-based statistics.

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   distance
   ordination
   spatial
   gradient
   power

Functions
---------

.. autosummary::
   :toctree: generated/

   p_value_to_str
   subsample
   subsample_counts
   isubsample

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numpy.testing import Tester

from ._misc import p_value_to_str
from ._subsample import subsample, subsample_counts, isubsample

__all__ = ['p_value_to_str', 'subsample', 'subsample_counts', 'isubsample']

test = Tester().test
