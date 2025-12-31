# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Multimodal analysis (:mod:`skbio.stats.multimodal`)
====================================================

.. currentmodule:: skbio.stats.multimodal

This module provides methods for multimodal data analysis, including techniques
for learning joint embeddings across multiple data modalities (e.g., microbes
and metabolites).


Functions
---------

.. autosummary::
   :toctree:

   mmvec


Classes
-------

.. autosummary::
   :toctree:

   MMvecResults

"""

from ._mmvec import mmvec, MMvecResults

__all__ = [
    "mmvec",
    "MMvecResults",
]
