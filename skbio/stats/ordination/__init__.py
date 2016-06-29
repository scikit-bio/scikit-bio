r"""
Ordination methods (:mod:`skbio.stats.ordination`)
==================================================

.. currentmodule:: skbio.stats.ordination

This module contains several ordination methods, including Principal
Coordinate Analysis, Correspondence Analysis, Redundancy Analysis and
Canonical Correspondence Analysis.


Functions
---------

.. autosummary::
   :toctree: generated/

   ca
   pcoa
   cca
   rda
   mean_and_std
   corr
   scale
   svd_rank
   e_matrix
   f_matrix

Classes
-------

.. autosummary::
   :toctree: generated/

   OrdinationResults

Examples
--------

This is an artificial dataset (table 11.3 in [1]_) that represents fish
abundance in different sites (`Y`, the response variables) and
environmental variables (`X`, the explanatory variables).

>>> import numpy as np
>>> import pandas as pd

First we need to construct our explanatory variable dataset `X`.

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
>>> transects = ['depth', 'substrate_coral', 'substrate_sand',
...              'substrate_other']
>>> sites = ['site1', 'site2', 'site3', 'site4', 'site5', 'site6', 'site7',
...          'site8', 'site9', 'site10']
>>> X = pd.DataFrame(X, sites, transects)

Then we need to create a dataframe with the information about the species
observed at different sites.

>>> species = ['specie1', 'specie2', 'specie3', 'specie4', 'specie5',
...            'specie6', 'specie7', 'specie8', 'specie9']
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
>>> Y = pd.DataFrame(Y, sites, species)

We can now perform canonical correspondence analysis. Matrix `X` contains a
continuous variable (depth) and a categorical one (substrate type) encoded
using a one-hot encoding.
>>> from skbio.stats.ordination import cca

We explicitly need to avoid perfect collinearity, so we'll drop one of the
substrate types (the last column of `X`).

>>> del X['substrate_other']
>>> ordination_result = cca(Y, X, scaling=2)

Exploring the results we see that the first three axes explain about
80% of all the variance.

>>> ordination_result.proportion_explained
CCA1    0.466911
CCA2    0.238327
CCA3    0.100548
CCA4    0.104937
CCA5    0.044805
CCA6    0.029747
CCA7    0.012631
CCA8    0.001562
CCA9    0.000532
dtype: float64

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

from skbio.util import TestRunner

from ._redundancy_analysis import rda
from ._correspondence_analysis import ca
from ._canonical_correspondence_analysis import cca
from ._principal_coordinate_analysis import pcoa
from ._ordination_results import OrdinationResults
from ._utils import (mean_and_std, scale, svd_rank, corr, e_matrix, f_matrix)

__all__ = ['ca', 'rda', 'cca', 'pcoa', 'OrdinationResults',
           'mean_and_std', 'scale', 'svd_rank', 'corr',
           'e_matrix', 'f_matrix']

test = TestRunner(__file__).test
