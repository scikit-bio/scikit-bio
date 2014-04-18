#!/usr/bin/env python
"""
Alpha diversity indices (:mod:`skbio.maths.diversity.alpha`)
============================================================

.. currentmodule:: skbio.maths.diversity.alpha

This package provides implementations of various alpha diversity indices. Some
functions also provide the option to generate confidence intervals (CIs).

Functions
---------

.. autosummary::
   :toctree: generated/

   ace
   berger_parker_d
   brillouin_d
   chao1
   chao1_confidence
   dominance
   doubles
   enspie
   equitability
   esty_ci
   fisher_alpha
   gini_index
   goods_coverage
   heip_e
   kempton_taylor_q
   lladser_ci
   lladser_pe
   margalef
   mcintosh_d
   mcintosh_e
   menhinick
   michaelis_menten_fit
   observed_species
   osd
   robbins
   shannon
   simpson
   simpson_e
   simpson_reciprocal
   singles
   strong

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .ace import ace
from .chao1 import chao1, chao1_confidence
from .base import (
    berger_parker_d, brillouin_d, dominance, doubles, enspie, equitability,
    esty_ci, fisher_alpha, goods_coverage, heip_e, kempton_taylor_q, margalef,
    mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, observed_species,
    osd, robbins, shannon, simpson, simpson_e, simpson_reciprocal, singles,
    strong)
from .gini import gini_index
from .lladser import lladser_pe, lladser_ci

__all__ = [m for m in dir() if not m.startswith('_')]

from numpy.testing import Tester
test = Tester().test
