"""
Distance matrices and distance-based statistics (:mod:`skbio.stats.distance`)
=============================================================================

.. currentmodule:: skbio.stats.distance

This subpackage provides functionality for serializing, deserializing, and
manipulating dissimilarity and distance matrices in memory. It also contains
various statistical methods that operate on distance matrices, often relating
distances (e.g., community distances) to categorical and/or continuous
variables of interest (e.g., gender or age). Methods are also provided for
comparing distance matrices (e.g., computing the correlation between two or
more distance matrices using the Mantel test).

Data Structures: DissimilarityMatrix and DistanceMatrix
-------------------------------------------------------

This package provides two matrix classes, `DissimilarityMatrix` and
`DistanceMatrix`. Both classes can store measures of difference/distinction
between objects. A dissimilarity/distance matrix includes both a matrix of
dissimilarities/distances (floats) between objects, as well as unique IDs
(object labels; strings) identifying each object in the matrix.

`DissimilarityMatrix` can be used to store measures of dissimilarity between
objects, and does not require that the dissimilarities are symmetric (e.g.,
dissimilarities obtained using the *Gain in PD* measure [1]_).
`DissimilarityMatrix` is a more general container to store differences than
`DistanceMatrix`.

`DistanceMatrix` has the additional requirement that the differences it
stores are symmetric (e.g., Euclidean or Hamming distances).

.. note:: `DissimilarityMatrix` can be used to store distances, but it is
   recommended to use `DistanceMatrix` to store this type of data as it
   provides an additional check for symmetry. A distance matrix *is a*
   dissimilarity matrix; this is modeled in the class design by having
   `DistanceMatrix` subclass `DissimilarityMatrix`.

Classes
^^^^^^^

.. autosummary::
   :toctree: generated/

   DissimilarityMatrix
   DistanceMatrix

Functions
^^^^^^^^^

.. autosummary::
   :toctree: generated/

   randdm

Exceptions
^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   DissimilarityMatrixError
   DistanceMatrixError
   MissingIDError

Examples
^^^^^^^^
Assume we have the following delimited text file storing distances between
three objects with IDs ``a``, ``b``, and ``c``::

    \\ta\\tb\\tc
    a\\t0.0\\t0.5\\t1.0
    b\\t0.5\\t0.0\\t0.75
    c\\t1.0\\t0.75\\t0.0

Load a distance matrix from the file:

>>> from StringIO import StringIO
>>> from skbio import DistanceMatrix
>>> dm_fh = StringIO("\\ta\\tb\\tc\\n"
...                  "a\\t0.0\\t0.5\\t1.0\\n"
...                  "b\\t0.5\\t0.0\\t0.75\\n"
...                  "c\\t1.0\\t0.75\\t0.0\\n")
>>> dm = DistanceMatrix.read(dm_fh)
>>> print(dm)
3x3 distance matrix
IDs:
'a', 'b', 'c'
Data:
[[ 0.    0.5   1.  ]
 [ 0.5   0.    0.75]
 [ 1.    0.75  0.  ]]

Access the distance (scalar) between objects ``'a'`` and ``'c'``:

>>> dm['a', 'c']
1.0

Get a row vector of distances between object ``'b'`` and all other objects:

>>> dm['b']
array([ 0.5 ,  0.  ,  0.75])

numpy indexing/slicing also works as expected. Extract the third column:

>>> dm[:, 2]
array([ 1.  ,  0.75,  0.  ])

Serialize the distance matrix to delimited text file:

>>> out_fh = StringIO()
>>> dm.write(out_fh)
>>> out_fh.getvalue() == dm_fh.getvalue()
True

A distance matrix object can also be created from an existing ``numpy.array``
(or an array-like object, such as a nested Python list):

>>> import numpy as np
>>> data = np.array([[0.0, 0.5, 1.0],
...                  [0.5, 0.0, 0.75],
...                  [1.0, 0.75, 0.0]])
>>> ids = ["a", "b", "c"]
>>> dm_from_np = DistanceMatrix(data, ids)
>>> print(dm_from_np)
3x3 distance matrix
IDs:
'a', 'b', 'c'
Data:
[[ 0.    0.5   1.  ]
 [ 0.5   0.    0.75]
 [ 1.    0.75  0.  ]]
>>> dm_from_np == dm
True

IDs may be omitted when constructing a dissimilarity/distance matrix.
Monotonically-increasing integers (cast as strings) will be automatically used:

>>> dm = DistanceMatrix(data)
>>> dm.ids
('0', '1', '2')

Distance-based statistics
-------------------------

In addition to the data structures described above, this package provides the
following distance-based statistical methods.

Categorical Variable Stats
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ANOSIM
   PERMANOVA
   CategoricalStatsResults

Continuous Variable Stats
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   bioenv

Distance Matrix Comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   mantel
   pwmantel

Examples
^^^^^^^^
Load a 4x4 distance matrix and grouping vector denoting 2 groups of objects.
Note that these statistical methods require symmetric distances:

>>> from skbio import DistanceMatrix
>>> dm = DistanceMatrix([[0, 1, 1, 4],
...                      [1, 0, 3, 2],
...                      [1, 3, 0, 3],
...                      [4, 2, 3, 0]],
...                     ['s1', 's2', 's3', 's4'])
>>> grouping = ['Group1', 'Group1', 'Group2', 'Group2']

Create an ANOSIM instance and run the method with 99 permutations:

>>> import numpy as np
>>> np.random.seed(0) # Make output deterministic; not necessary for normal use
>>> from skbio.stats.distance import ANOSIM
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

References
----------
.. [1] Faith, D. P. (1992). "Conservation evaluation and phylogenetic
   diversity".

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._base import (DissimilarityMatrixError, DistanceMatrixError,
                    MissingIDError, DissimilarityMatrix, DistanceMatrix,
                    CategoricalStatsResults, randdm)
from ._bioenv import bioenv
from ._anosim import ANOSIM
from ._permanova import PERMANOVA
from ._mantel import mantel, pwmantel

__all__ = ['DissimilarityMatrixError', 'DistanceMatrixError', 'MissingIDError',
           'DissimilarityMatrix', 'DistanceMatrix', 'randdm', 'ANOSIM',
           'PERMANOVA', 'CategoricalStatsResults', 'bioenv', 'mantel',
           'pwmantel']

from numpy.testing import Tester
test = Tester().test
