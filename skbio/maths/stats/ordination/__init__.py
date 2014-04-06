r"""
Ordination methods (:mod:`skbio.maths.stats.ordination`)
========================================================

.. currentmodule:: skbio.maths.stats.ordination

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

__all__ = ['CA', 'RDA', 'CCA', 'PCoA']
