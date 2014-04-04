#!/usr/bin/env python
r"""
Alpha diversity measures (:mod:`skbio.maths.diversity.alpha`)
=============================================================

.. currentmodule:: skbio.maths.diversity.alpha

This module provides implementations of various alpha diversity measures.

Functions
---------

.. autosummary::
   :toctree: generated/

   observed_species

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np


def berger_parker_d(counts):
    """Fraction of the sample that belongs to the most abundant species.

    References
    ----------
    .. [1] Berger & Parker 1970, by way of SDR-IV online help.

    """
    return counts.max() / counts.sum()

def observed_species(counts):
    """Calculate number of distinct species."""
    return (counts != 0).sum()

def osd(counts):
    """Calculate observed, singles and doubles from counts."""
    return observed_species(counts), (counts == 1).sum(), (counts == 2).sum()
