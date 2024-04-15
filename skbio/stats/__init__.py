r"""Multivariate Statistics (:mod:`skbio.stats`)
============================================

.. currentmodule:: skbio.stats

This module provides various statistical methods to support the analyses of
high-dimensional biological data to uncover the relationships among samples,
features and metadata. Examples include distance matrix-based statistics,
ordination methods, composition statistics, and data subsampling techniques.


Distance matrix statistics
--------------------------

.. autosummary::
   :toctree: generated/

   distance


Ordination methods
------------------

.. autosummary::
   :toctree: generated/

   ordination


Composition statistics
----------------------

.. autosummary::
   :toctree: generated/

   composition


Data subsampling
----------------

.. autosummary::
   :toctree: generated/

   subsample_counts
   isubsample


Other statistical methods
-------------------------

.. autosummary::
   :toctree: generated/

   evolve
   gradient
   power

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
