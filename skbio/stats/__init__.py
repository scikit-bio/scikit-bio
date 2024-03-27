"""Statistics (:mod:`skbio.stats`)
===============================

.. currentmodule:: skbio.stats

This module contains various statistical methods, such as distance matrix-based
statistics, ordination techniques, composition statistics, and data subsampling
techniques.


Distance matrix-based statistics
--------------------------------

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
