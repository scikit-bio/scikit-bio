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
   :toctree:

   DissimilarityMatrix
   DistanceMatrix

Functions
^^^^^^^^^

.. autosummary::
   :toctree:

   randdm

Exceptions
^^^^^^^^^^

.. autosummary::
   :toctree:

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

>>> from io import StringIO
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
>>> _ = dm.write(out_fh)
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
   :toctree:

   anosim
   permanova
   permdisp

Continuous Variable Stats
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree:

   bioenv

Distance Matrix Comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree:

   mantel
   pwmantel

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
                    randdm)
from ._bioenv import bioenv
from ._anosim import anosim
from ._permanova import permanova
from ._mantel import mantel, pwmantel
from ._permdisp import permdisp

__all__ = ['DissimilarityMatrixError', 'DistanceMatrixError', 'MissingIDError',
           'DissimilarityMatrix', 'DistanceMatrix', 'randdm', 'anosim',
           'permanova', 'bioenv', 'mantel', 'pwmantel', 'permdisp']
