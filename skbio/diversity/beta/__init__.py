"""Beta diversity measures (:mod:`skbio.diversity.beta`)
=====================================================

.. currentmodule:: skbio.diversity.beta

This package provides implementations of beta diversity measures for computing sample
dissimilarity. Users should also explore SciPy's :func:`~scipy.spatial.distance.pdist`,
which provides implementations of additional beta diversity metrics with interfaces
similar to those provided here, and can be directly used in the driver function
:func:`~skbio.diversity.beta_diversity`.

Functions
---------

.. autosummary::
   :toctree:

    unweighted_unifrac
    weighted_unifrac

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._unifrac import unweighted_unifrac, weighted_unifrac

__all__ = ["unweighted_unifrac", "weighted_unifrac"]
