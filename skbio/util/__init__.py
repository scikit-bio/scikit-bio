"""
Utility functionality (:mod:`skbio.util`)
=========================================

.. currentmodule:: skbio.util

This package provides several subpackages containing various utility
functionality, including custom exception/warning definitions and I/O and
unit-testing convenience functions.

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   exception
   warning
   io
   misc
   testing

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from numpy.testing import Tester
test = Tester().test
