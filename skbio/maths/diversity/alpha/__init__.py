#!/usr/bin/env python
"""
Alpha diversity indices (:mod:`skbio.maths.diversity.alpha`)
============================================================

.. currentmodule:: skbio.maths.diversity.alpha

This package provides implementations of various alpha diversity indices. Some
functions also provide the option to generate confidence intervals (CIs).

Functions
---------

.. autosummary::
   :toctree: generated/

   ace
   observed_species

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .ace import ace
from .chao1 import chao1, chao1_confidence
from .base import observed_species
from .gini import gini_index
from .lladser import lladser_pe, lladser_ci

# Need a better way to handle __all__ than updating by hand
#__all__ = ['ace', 'observed_species']

from numpy.testing import Tester
test = Tester().test
