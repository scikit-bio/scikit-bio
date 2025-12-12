r"""Optimized libraries from scikit-bio-binaries
============================================

.. currentmodule:: skbio.binaries

This module provides an interface to the optimized functions
that are available in the scikit-bio-binaries package,
via the libskbb shared library.

Since scikit-bio-binaries is not guaranteed to be installed
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
    available,
    get_api_version,
    set_random_seed,
    py_to_bin_random_seed,
)
from ._ordination import (
    pcoa_fsvd_available,
    pcoa_fsvd,
)

from ._distance import (
    permanova_available,
    permanova,
)

__all__ = [
    "available",
    "get_api_version",
    "set_random_seed",
    "py_to_bin_random_seed",
    "pcoa_fsvd_available",
    "pcoa_fsvd",
    "permanova_available",
    "permanova",
]
