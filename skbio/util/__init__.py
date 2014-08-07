"""
Utility functionality (:mod:`skbio.util`)
=========================================

.. currentmodule:: skbio.util

This package provides general exception/warning definitions used throughout
scikit-bio, as well as several subpackages containing various utility
functionality, including I/O and unit-testing convenience functions.

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   misc
   testing

Warnings
--------

.. autosummary::
   :toctree: generated/

   EfficiencyWarning

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from ._warning import EfficiencyWarning

__all__ = ['EfficiencyWarning']

from numpy.testing import Tester
test = Tester().test
