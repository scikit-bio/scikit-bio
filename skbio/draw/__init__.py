"""
Visualizations (:mod:`skbio.draw`)
==================================

.. currentmodule:: skbio.draw

This module provides functionality for visualization of data.

Distribution visualizations
---------------------------

Functions
^^^^^^^^^

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

from numpy.testing import Tester

from ._distributions import boxplots, grouped_distributions

__all__ = ['boxplots', 'grouped_distributions']

test = Tester().test
