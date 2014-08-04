r"""
Ordination methods (:mod:`skbio.stats.ordination`)
==================================================

.. currentmodule:: skbio.stats.ordination

This module contains several ordination methods, including Principal
Coordinate Analysis, Correspondence Analysis, Redundancy Analysis and
Canonical Correspondence Analysis.

Classes
-------

.. autosummary::
   :toctree: generated/

   PCoA
   CA
   RDA
   CCA
   OrdinationResults

Examples
--------

This is an artificial dataset (table 11.3 in [1]_) that represents fish
abundance in different sites (`Y`, the response variables) and
environmental variables (`X`, the explanatory variables).

>>> import numpy as np
>>> X = np.array([[1.0, 0.0, 1.0, 0.0],
...               [2.0, 0.0, 1.0, 0.0],
...               [3.0, 0.0, 1.0, 0.0],
...               [4.0, 0.0, 0.0, 1.0],
...               [5.0, 1.0, 0.0, 0.0],
...               [6.0, 0.0, 0.0, 1.0],
...               [7.0, 1.0, 0.0, 0.0],
...               [8.0, 0.0, 0.0, 1.0],
...               [9.0, 1.0, 0.0, 0.0],
...               [10.0, 0.0, 0.0, 1.0]])
>>> Y = np.array([[1, 0, 0, 0, 0, 0, 2, 4, 4],
...               [0, 0, 0, 0, 0, 0, 5, 6, 1],
...               [0, 1, 0, 0, 0, 0, 0, 2, 3],
...               [11, 4, 0, 0, 8, 1, 6, 2, 0],
...               [11, 5, 17, 7, 0, 0, 6, 6, 2],
...               [9, 6, 0, 0, 6, 2, 10, 1, 4],
...               [9, 7, 13, 10, 0, 0, 4, 5, 4],
...               [7, 8, 0, 0, 4, 3, 6, 6, 4],
...               [7, 9, 10, 13, 0, 0, 6, 2, 0],
...               [5, 10, 0, 0, 2, 4, 0, 1, 3]])

We can now create a CCA object to perform canonical correspondence
analysis. Matrix `X` contains a continuous variable (depth) and a
categorical one (substrate type) encoded using a one-hot encoding. We
explicitly need to avoid perfect collinearity, so we'll drop one of
the substrate types (the last column of `X`). We also expect to
increase pandas integration to ease analyses.

>>> from skbio.stats.ordination import CCA
>>> ordination_result = CCA(Y, X[:, :-1],
...                         ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
...                          'Site5', 'Site6', 'Site7', 'Site8', 'Site9'],
...                         ['Species0', 'Species1', 'Species2', 'Species3',
...                          'Species4', 'Species5', 'Species6', 'Species7',
...                          'Species8'])

Exploring the results we see that the first three axes explain about
80% of all the variance.

>>> sc_2 = ordination_result.scores(scaling=2)
>>> print sc_2.proportion_explained
[ 0.46691091  0.23832652  0.10054837  0.10493671  0.04480535  0.02974698
  0.01263112  0.00156168  0.00053235]

References
----------

.. [1] Legendre P. and Legendre L. 1998. Numerical Ecology. Elsevier,
   Amsterdam.

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .correspondence_analysis import CA
from .redundancy_analysis import RDA
from .canonical_correspondence_analysis import CCA
from .principal_coordinate_analysis import PCoA
from .base import OrdinationResults

__all__ = ['CA', 'RDA', 'CCA', 'PCoA', 'OrdinationResults']

from numpy.testing import Tester
test = Tester().test
