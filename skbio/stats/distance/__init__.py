r"""Distance matrix-based statistics (:mod:`skbio.stats.distance`)
==============================================================

.. currentmodule:: skbio.stats.distance

This module provides functionality for storing, manipulating, and analyzing pairwise
relationships, especially distances, between biological objects. It contains various
statistical methods that operate on distance matrices, often relating distances to
categorical and/or numeric variables of interest (e.g., gender and/or age). Methods
are also provided for comparing distance matrices (e.g., computing the correlation
between two or more distance matrices using the Mantel test).


Data structures
---------------

This package provides three hierarchical matrix classes: ``PairwiseMatrix`` >
``SymmetricMatrix`` > ``DistanceMatrix``, to store measures of pairwise relationships
between objects. A matrix object stores a 2-D float array of pairwise measures and a
tuple of unique string labels identifying each object in the matrix.

``PairwiseMatrix`` can be used to store any measures of pairwise relationships. Its
subclass ``SymmetricMatrix`` requires that the relationships are **symmetric** (i.e.,
:math:`D(i, j) = D(j, i)`). ``DistanceMatrix``, the most specific class of the three,
further requires that the matrix is **hollow** (i.e., :math:`D(i, i) = 0`), usually
representing the difference/distinction between objects (e.g., Euclidean or Hamming
distances). However, ``DistanceMatrix`` does not enforce non-negativity or triangle
inequality, two properties of a metric distance, therefore it is suitable for more
general measures (e.g., Bray-Curtis dissimilarity).

Classes
^^^^^^^

.. autosummary::
   :toctree:

   PairwiseMatrix
   SymmetricMatrix
   DistanceMatrix

Functions
^^^^^^^^^

.. autosummary::
   :toctree:

   randdm

Exceptions
^^^^^^^^^^

.. autosummary::

   PairwiseMatrixError
   SymmetricMatrixError
   DistanceMatrixError
   MissingIDError


Distance-based statistics
-------------------------

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


Examples
--------
Assume we have the following delimited text file storing distances between
three objects with IDs ``a``, ``b``, and ``c``::

    \ta\tb\tc
    a\t0.0\t0.5\t1.0
    b\t0.5\t0.0\t0.75
    c\t1.0\t0.75\t0.0

Load a distance matrix from the file:

>>> from io import StringIO
>>> from skbio.stats.distance import DistanceMatrix
>>> dm_fh = StringIO("\ta\tb\tc\n"
...                  "a\t0.0\t0.5\t1.0\n"
...                  "b\t0.5\t0.0\t0.75\n"
...                  "c\t1.0\t0.75\t0.0\n")
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

IDs may be omitted when constructing a pairwise/symmetric/distance matrix.
Monotonically-increasing integers (cast as strings) will be automatically used:

>>> dm = DistanceMatrix(data)
>>> dm.ids
('0', '1', '2')

"""  # noqa: D407, D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._base import (
    PairwiseMatrixError,
    DissimilarityMatrixError,
    SymmetricMatrixError,
    DistanceMatrixError,
    MissingIDError,
    PairwiseMatrix,
    DissimilarityMatrix,
    SymmetricMatrix,
    DistanceMatrix,
    randdm,
)
from ._bioenv import bioenv
from ._anosim import anosim
from ._permanova import permanova
from ._mantel import mantel, pwmantel
from ._permdisp import permdisp

__all__ = [
    "PairwiseMatrixError",
    "DissimilarityMatrixError",
    "SymmetricMatrixError",
    "DistanceMatrixError",
    "MissingIDError",
    "PairwiseMatrix",
    "DissimilarityMatrix",
    "SymmetricMatrix",
    "DistanceMatrix",
    "randdm",
    "anosim",
    "permanova",
    "bioenv",
    "mantel",
    "pwmantel",
    "permdisp",
]
