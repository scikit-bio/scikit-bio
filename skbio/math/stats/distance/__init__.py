#!/usr/bin/env python
"""
Distance-based statistics (:mod:`skbio.math.stats.distance`)
============================================================

.. currentmodule:: skbio.math.stats.distance

This package contains various statistical methods that operate on distance
matrices, often relating distances (e.g., community distances) to categorical
and/or continuous variables of interest (e.g., gender or age).

Categorical Variable Stats
--------------------------

.. autosummary::
   :toctree: generated/

   ANOSIM
   PERMANOVA
   CategoricalStatsResults

Continuous Variable Stats
-------------------------

.. autosummary::
   :toctree: generated/

   bioenv

Examples
--------
Load a 4x4 distance matrix and grouping vector denoting 2 groups of objects.
Note that these statistical methods require symmetric distances:

>>> from skbio.core.distance import DistanceMatrix
>>> dm = DistanceMatrix([[0, 1, 1, 4],
...                      [1, 0, 3, 2],
...                      [1, 3, 0, 3],
...                      [4, 2, 3, 0]],
...                     ['s1', 's2', 's3', 's4'])
>>> grouping = ['Group1', 'Group1', 'Group2', 'Group2']

Create an ANOSIM instance and run the method with 99 permutations:

>>> import numpy as np
>>> np.random.seed(0) # Make output deterministic; not necessary for normal use
>>> from skbio.math.stats.distance import ANOSIM
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

>>> results.statistic
0.25
>>> print results.p_value
None
>>> results.permutations
0

You can also provide a ``pandas.DataFrame`` and a column denoting the grouping
instead of a grouping vector. The following data frame's ``Group`` column
specifies the same grouping as the vector we used in all of the previous
examples:

>>> np.random.seed(0) # Make output deterministic; not necessary for normal use
>>> import pandas as pd
>>> df = pd.DataFrame.from_dict(
...     {'Group': {'s2': 'Group1', 's3': 'Group2', 's4': 'Group2',
...                's5': 'Group3', 's1': 'Group1'}})
>>> anosim = ANOSIM(dm, df, column='Group')
>>> results = anosim(99)
>>> print results
Method name  Sample size  Number of groups  R statistic  p-value  \
Number of permutations
     ANOSIM            4                 2         0.25     0.67  \
                    99
<BLANKLINE>

The results match the results we saw in the first example above.

Note that when providing a data frame, the ordering of rows and/or columns does
not affect the grouping vector that is extracted. The data frame must be
indexed by the distance matrix IDs (i.e., the row labels must be distance
matrix IDs).

If IDs (rows) are present in the data frame but not in the distance matrix,
they are ignored (the previous example's ``s5`` ID illustrates this behavior).
Thus, the data frame can be a superset of the distance matrix IDs. Note that
the reverse is not true: IDs in the distance matrix *must* be present in the
data frame or an error will be raised.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .base import CategoricalStatsResults, bioenv
from .anosim import ANOSIM
from .permanova import PERMANOVA
from ._mantel import mantel

__all__ = ['ANOSIM', 'PERMANOVA', 'CategoricalStatsResults', 'bioenv',
           'mantel']

from numpy.testing import Tester
test = Tester().test
