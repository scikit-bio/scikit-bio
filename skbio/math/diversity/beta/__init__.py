"""
Beta diversity measures (:mod:`skbio.math.diversity.beta`)
==========================================================

.. currentmodule:: skbio.math.diversity.beta

This package contains helper functions for working with scipy's pairwise
distance (``pdist``) functions in scikit-bio, and will eventually be expanded
to contain pairwise distance/dissimilarity methods that are not implemented
(or planned to be implemented) in scipy.

The functions in this package currently support applying ``pdist`` functions
to all pairs of samples in a sample by observation count or abundance matrix
and returning an ``skbio.DistanceMatrix`` object. This application is
illustrated below for a few different forms of input.

Functions
---------

.. autosummary::
   :toctree: generated/

    pw_distances
    pw_distances_from_table

Examples
--------
Create a table object containing 7 OTUs and 6 samples.

>>> from skbio.math.diversity.beta import pw_distances
>>> import numpy as np
>>> data = [[23, 64, 14, 0, 0, 3, 1],
...         [0, 3, 35, 42, 0, 12, 1],
...         [0, 5, 5, 0, 40, 40, 0],
...         [44, 35, 9, 0, 1, 0, 0],
...         [0, 2, 8, 0, 35, 45, 1],
...         [0, 0, 25, 35, 0, 19, 0]]
>>> ids = list('ABCDEF')

Compute Bray-Curtis distances between all pairs of samples and return a
DistanceMatrix object.

>>> bc_dm = pw_distances(data, ids, "braycurtis")
>>> print(bc_dm)
6x6 distance matrix
IDs:
A, B, C, D, E, F
Data:
[[ 0.          0.78787879  0.86666667  0.30927835  0.85714286  0.81521739]
 [ 0.78787879  0.          0.78142077  0.86813187  0.75        0.1627907 ]
 [ 0.86666667  0.78142077  0.          0.87709497  0.09392265  0.71597633]
 [ 0.30927835  0.86813187  0.87709497  0.          0.87777778  0.89285714]
 [ 0.85714286  0.75        0.09392265  0.87777778  0.          0.68235294]
 [ 0.81521739  0.1627907   0.71597633  0.89285714  0.68235294  0.        ]]

Compute Jaccard distances between all pairs of samples and return a
DistanceMatrix object.

>>> j_dm = pw_distances(data, ids, "jaccard")
>>> print(j_dm)
6x6 distance matrix
IDs:
A, B, C, D, E, F
Data:
[[ 0.          0.83333333  1.          1.          0.83333333  1.        ]
 [ 0.83333333  0.          1.          1.          0.83333333  1.        ]
 [ 1.          1.          0.          1.          1.          1.        ]
 [ 1.          1.          1.          0.          1.          1.        ]
 [ 0.83333333  0.83333333  1.          1.          0.          1.        ]
 [ 1.          1.          1.          1.          1.          0.        ]]

Determine if the resulting distance matrices are significantly correlated
by computing the Mantel correlation between them. Then determine if the p-value
is significant based on an alpha of 0.05.

>>> from skbio.math.stats.distance import mantel
>>> r, p_value = mantel(j_dm, bc_dm)
>>> print(r)
-0.209362157621
>>> print(p_value < 0.05)
False

Compute PCoA for both distance matrices, and then find the Procrustes
M-squared value that results from comparing the coordinate matrices.

>>> from skbio.math.stats.ordination import PCoA
>>> bc_pc = PCoA(bc_dm).scores()
>>> j_pc = PCoA(j_dm).scores()
>>> from skbio.math.stats.spatial import procrustes
>>> print procrustes(bc_pc.site, j_pc.site)[2]
0.466134984787

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .base import (
    pw_distances, pw_distances_from_table)

__all__ = ["pw_distances", "pw_distances_from_table"]

from numpy.testing import Tester
test = Tester().test
