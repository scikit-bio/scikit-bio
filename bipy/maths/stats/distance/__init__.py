#!/usr/bin/env python
"""
Distance-based statistics (:mod:`bipy.maths.stats.distance`)
============================================================

.. currentmodule:: bipy.maths.stats.distance

Distance-based statistical methods package.

Classes
-------

.. autosummary::
   :toctree: generated/

   ANOSIM
   PERMANOVA
   CategoricalStatsResults

Examples
--------
Load a 4x4 distance matrix and grouping vector denoting 2 groups of objects.
Note that these statistical methods require symmetric distances:

>>> from bipy.core.distance import SymmetricDistanceMatrix
>>> dm = SymmetricDistanceMatrix([[0, 1, 1, 4],
...                               [1, 0, 3, 2],
...                               [1, 3, 0, 3],
...                               [4, 2, 3, 0]],
...                              ['s1', 's2', 's3', 's4'])
>>> grouping = ['Group1', 'Group1', 'Group2', 'Group2']

Create an ANOSIM instance and run the method with 99 permutations:

>>> import numpy as np
>>> np.random.seed(0) # Make output deterministic; not necessary for normal use
>>> from bipy.maths.stats.distance import ANOSIM
>>> anosim = ANOSIM(dm, grouping)
>>> results = anosim(99)
>>> print results.summary() # doctest: +NORMALIZE_WHITESPACE
Method name\tSample size\tNumber of groups\tR statistic\tp-value\t\
Number of permutations
ANOSIM\t4\t2\t0.25\t0.67\t99

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .base import CategoricalStatsResults
from .anosim import ANOSIM
from .permanova import PERMANOVA

__all__ = ['ANOSIM', 'PERMANOVA', 'CategoricalStatsResults']
