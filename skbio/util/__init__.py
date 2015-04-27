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
   _skip_if_no_matplotlib
   TestRunner

Miscellaneous functionality
---------------------------

Generally useful functions that don't fit in more specific locations.

.. autosummary::
   :toctree: generated/

   cardinal_to_ordinal
   create_dir
   find_duplicates
   flatten
   is_casava_v180_or_later
   remove_files
   safe_md5

Exceptions
----------

.. autosummary::
   :toctree: generated/

   TestingUtilError

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
# ----------------------------------------------------------------------------

from ._warning import EfficiencyWarning
from ._exception import TestingUtilError
from ._misc import (cardinal_to_ordinal, create_dir, find_duplicates, flatten,
                    is_casava_v180_or_later, remove_files, safe_md5)
from ._testing import (get_data_path, _skip_if_no_matplotlib,
                       _not_has_matplotlib, TestRunner)

__all__ = ['EfficiencyWarning', 'TestingUtilError',
           'cardinal_to_ordinal', 'create_dir', 'find_duplicates', 'flatten',
           'is_casava_v180_or_later', 'remove_files', 'safe_md5',
           'get_data_path', 'TestRunner', '_skip_if_no_matplotlib',
           '_not_has_matplotlib']

test = TestRunner(__file__).test
