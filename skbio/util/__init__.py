"""
Utility functionality (:mod:`skbio.util`)
=========================================

.. currentmodule:: skbio.util

This package provides general exception/warning definitions used throughout
scikit-bio, as well as various utility functionality, including I/O and
unit-testing convenience functions.

Testing functionality
---------------------

Common functionality to support testing in skbio.

.. autosummary::
   :toctree: generated/

   get_data_path

Miscellaneous functionality
---------------------------

Generally useful functions that don't fit in more specific locations.

.. autosummary::
   :toctree: generated/

   create_dir
   flatten
   remove_files
   safe_md5
   is_casava_v180_or_later

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
from ._misc import (create_dir, flatten, remove_files, safe_md5,
                    is_casava_v180_or_later)
from ._testing import get_data_path

__all__ = ['FileFormatError', 'RecordError', 'FieldError', 'EfficiencyWarning',
           'create_dir', 'flatten', 'remove_files', 'safe_md5',
           'is_casava_v180_or_later', 'get_data_path']

from numpy.testing import Tester
test = Tester().test
