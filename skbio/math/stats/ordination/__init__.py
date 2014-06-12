r"""
Ordination methods (:mod:`skbio.math.stats.ordination`)
=======================================================

.. currentmodule:: skbio.math.stats.ordination

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

Let's import a few utility functions to load the data:

>>> from skbio.util.testing import get_data_path
>>> import os
>>> from functools import partial
>>> path = partial(get_data_path,
...                subfolder=os.path.join('skbio', 'math', 'stats',
...                                       'ordination', 'tests', 'data'))

It's an artificial dataset (table 11.3 in [1]_) that represents fish
abundance in different sites (`Y`, the response variables) and
environmental variables (`X`, the explanatory variables).

>>> import numpy as np
>>> X = np.loadtxt(path('example3_X'))
>>> Y = np.loadtxt(path('example3_Y'))

We can now create a CCA object to perform canonical correspondence
analysis. Matrix `X` contains a continuous variable (depth) and a
categorical one (substrate type) encoded using a one-hot encoding. We
explicitly need to avoid perfect collinearity, so we'll drop one of
the substrate types (the last column of `X`). We also expect to
increase pandas integration to ease analyses.

>>> from skbio.math.stats.ordination import CCA
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
