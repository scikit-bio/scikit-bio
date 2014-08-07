#!/usr/bin/env python

"""
Visualizations (:mod:`skbio.draw`)
==================================

.. currentmodule:: skbio.draw.distributions

This module provides functionality for visualization of data.

Distribution visualizations
---------------------------

Functions
---------

.. autosummary::
   :toctree: generated/

   boxplots
   grouped_distributions

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._distributions import boxplots, grouped_distributions

__all__ = ['boxplots', 'grouped_distributions']

from numpy.testing import Tester
test = Tester().test
