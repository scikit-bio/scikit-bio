#!/usr/bin/env python
"""
Distance-based statistics (:mod:`skbio.maths.stats.distance`)
=============================================================

.. currentmodule:: skbio.maths.stats.distance

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

>>> from skbio.core.distance import SymmetricDistanceMatrix
>>> dm = SymmetricDistanceMatrix([[0, 1, 1, 4],
...                               [1, 0, 3, 2],
...                               [1, 3, 0, 3],
...                               [4, 2, 3, 0]],
...                              ['s1', 's2', 's3', 's4'])
>>> grouping = ['Group1', 'Group1', 'Group2', 'Group2']

Create an ANOSIM instance and run the method with 99 permutations:

>>> import numpy as np
>>> np.random.seed(0) # Make output deterministic; not necessary for normal use
>>> from skbio.maths.stats.distance import ANOSIM
>>> anosim = ANOSIM(dm, grouping)
>>> results = anosim(99)
>>> print results
Method name  Sample size  Number of groups  R statistic  p-value  \
Number of permutations
     ANOSIM            4                 2         0.25     0.67  \
                    99
<BLANKLINE>

It is possible to rerun a method using an existing instance. Rerun ANOSIM with
999 permutations this time. Note that we obtain the same R statistic as before:

>>> results = anosim(999)
>>> print results
Method name  Sample size  Number of groups  R statistic  p-value  \
Number of permutations
     ANOSIM            4                 2         0.25    0.667  \
                   999
<BLANKLINE>

To suppress calculation of the p-value and only obtain the R statistic, specify
zero permutations:

>>> results = anosim(0)
>>> print results
Method name  Sample size  Number of groups  R statistic  p-value  \
Number of permutations
     ANOSIM            4                 2         0.25      N/A  \
                     0
<BLANKLINE>

A statistical results object can also format its results as delimited text.
This is useful, for example, if you want to view the results in a spreadsheet
program such as Excel:

>>> print results.summary(delimiter=',')
Method name,Sample size,Number of groups,R statistic,p-value,\
Number of permutations
ANOSIM,4,2,0.25,N/A,0
<BLANKLINE>

Individual values of the results can be accessed via the attributes of the
``CategoricalStatsResults`` class:

>>> results = anosim(99)
>>> results.statistic
0.25
>>> results.p_value
0.75
>>> results.permutations
99

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .base import CategoricalStatsResults
from .anosim import ANOSIM
from .permanova import PERMANOVA

__all__ = ['ANOSIM', 'PERMANOVA', 'CategoricalStatsResults']

from numpy.testing import Tester
test = Tester().test
