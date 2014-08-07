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

   io
   misc
   testing

Exceptions
----------

.. autosummary::
   :toctree: generated/

   FileFormatError
   RecordError
   FieldError

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

from ._exception import FileFormatError, RecordError, FieldError
from ._warning import EfficiencyWarning

__all__ = ['FileFormatError', 'RecordError', 'FieldError', 'EfficiencyWarning']

from numpy.testing import Tester
test = Tester().test
