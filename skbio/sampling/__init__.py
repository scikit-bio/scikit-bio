r"""Sampling methods (:mod:`skbio.sampling`)
========================================

.. currentmodule:: skbio.sampling

This module provides methods for sampling from frequency data.

Functions
---------

.. autosummary::
   :toctree: generated/

   subsample_counts
   isubsample

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._subsample import subsample_counts, isubsample


__all__ = ["subsample_counts", "isubsample"]
