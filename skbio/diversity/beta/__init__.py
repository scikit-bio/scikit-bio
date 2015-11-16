"""
Beta diversity measures (:mod:`skbio.diversity.beta`)
=====================================================

.. currentmodule:: skbio.diversity.beta

This package contains helper functions for working with scipy's pairwise
distance (``pdist``) functions in scikit-bio, and will eventually be expanded
to contain pairwise distance/dissimilarity methods that are not implemented
(or planned to be implemented) in scipy.

The functions in this package currently support applying ``pdist`` functions
to all pairs of samples in a sample by observation count or abundance matrix
and returning an ``skbio.DistanceMatrix`` object. This application is
illustrated below for a few different forms of input.

Functions
---------

.. autosummary::
   :toctree: generated/

    unweighted_unifrac
    weighted_unifrac

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import TestRunner

from ._unifrac import unweighted_unifrac, weighted_unifrac

__all__ = ["unweighted_unifrac", "weighted_unifrac"]

test = TestRunner(__file__).test
