"""
Evolutionary statistics (:mod:`skbio.math.stats.evolve`)
========================================================

.. currentmodule:: skbio.math.stats.evolve

This package contains statistics pertaining to phylogenies and evolution. 

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   hommola

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .hommola import hommola_cospeciation

from numpy.testing import Tester
test = Tester().test
