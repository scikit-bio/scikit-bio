r"""Utilities for Developers (:mod:`skbio.util`)
============================================

.. currentmodule:: skbio.util

This package provides general exception/warning definitions used throughout
scikit-bio, as well as various utility functionality, including I/O and
unit-testing convenience functions.


User configuration utilities
----------------------------

Configuration options for functional input and output.

.. autosummary::
   :toctree: generated/

   config
   get_config
   set_config


Testing utilities
-----------------

Common functionality to support testing in skbio.

.. autosummary::
   :toctree: generated/

   get_data_path
   assert_ordination_results_equal
   assert_ordination_results_equal_np
   assert_data_frame_almost_equal


Plotting utilities
------------------

.. autosummary::
   :toctree: generated/

   PlottableMixin


Decorators
----------

.. autosummary::
   :toctree: generated/

   overrides
   classproperty
   classonlymethod
   deprecated
   aliased
   register_aliases
   params_aliased


Miscellaneous utilities
-----------------------

Generally useful functionality that doesn't fit in more specific locations.

.. autosummary::
   :toctree: generated/

   get_rng
   cardinal_to_ordinal
   find_duplicates
   safe_md5

"""  # noqa: D412, D416, D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._misc import cardinal_to_ordinal, find_duplicates, safe_md5, get_rng
from ._testing import (
    get_data_path,
    assert_ordination_results_equal,
    assert_ordination_results_equal_np,
    assert_data_frame_almost_equal,
    pytestrunner,
)
from ._decorator import (
    overrides,
    classproperty,
    classonlymethod,
    deprecated,
    aliased,
    register_aliases,
    params_aliased,
)
from ._plotting import PlottableMixin
from .config._config import get_config, set_config

__all__ = [
    "cardinal_to_ordinal",
    "find_duplicates",
    "safe_md5",
    "get_rng",
    "get_data_path",
    "assert_ordination_results_equal",
    "assert_ordination_results_equal_np",
    "assert_data_frame_almost_equal",
    "pytestrunner",
    "PlottableMixin",
    "overrides",
    "classproperty",
    "classonlymethod",
    "deprecated",
    "aliased",
    "register_aliases",
    "params_aliased",
    "get_config",
    "set_config",
]
