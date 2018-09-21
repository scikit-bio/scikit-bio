"""
Statistics (:mod:`skbio.stats`)
===============================

.. currentmodule:: skbio.stats

This package contains various statistical methods, including ordination
techniques and distance matrix-based statistics.

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   distance
   evolve
   ordination
   gradient
   power
   composition

Functions
---------

.. autosummary::
   :toctree: generated/

   subsample_counts
   isubsample

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._subsample import subsample_counts, isubsample

__all__ = ['subsample_counts', 'isubsample']
