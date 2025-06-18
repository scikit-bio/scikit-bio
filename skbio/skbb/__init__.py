r"""Optimized libraries from scikit-bio-binaries
============================================

.. currentmodule:: skbio.skbb

This module provides an interface to the optimized functions
that are avaialble in the scikit-bio-binaries package, 
via the libskbb shared library.

Since scikit-bio-binaries is not guaragteed to be installed
the module provides ways to check for the presence of the 
needed shared library.

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2025--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._util import (
    skbb_available,
    skbb_get_api_version,
    skbb_set_random_seed,
)
from ._ordination import (
    skbb_pcoa_fsvd,
)

__all__ = [
    "skbb_available",
    "skbb_get_api_version",
    "skbb_set_random_seed",
    "skbb_pcoa_fsvd",
]
