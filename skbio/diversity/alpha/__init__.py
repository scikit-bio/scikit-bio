"""
Alpha diversity measures (:mod:`skbio.diversity.alpha`)
=======================================================

.. currentmodule:: skbio.diversity.alpha

This package provides implementations of alpha diversity measures, including
measures of richness, dominance, and evenness. Some functions generate
confidence intervals (CIs). These functions have the suffix ``_ci``.

Functions
---------

.. autosummary::
   :toctree:

   ace
   berger_parker_d
   brillouin_d
   chao1
   chao1_ci
   dominance
   doubles
   enspie
   esty_ci
   faith_pd
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
   observed_otus
   osd
   pielou_e
   robbins
   shannon
   simpson
   simpson_e
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

from ._ace import ace
from ._chao1 import chao1, chao1_ci
from ._faith_pd import faith_pd
from ._base import (
    berger_parker_d, brillouin_d, dominance, doubles, enspie,
    esty_ci, fisher_alpha, goods_coverage, heip_e, kempton_taylor_q,
    margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit,
    observed_otus, osd, pielou_e, robbins, shannon, simpson, simpson_e,
    singles, strong)
from ._gini import gini_index
from ._lladser import lladser_pe, lladser_ci


__all__ = ['ace', 'chao1', 'chao1_ci', 'berger_parker_d',
           'brillouin_d', 'dominance', 'doubles', 'enspie', 'esty_ci',
           'faith_pd', 'fisher_alpha', 'gini_index', 'goods_coverage',
           'heip_e', 'kempton_taylor_q', 'margalef', 'mcintosh_d',
           'mcintosh_e', 'menhinick', 'michaelis_menten_fit', 'observed_otus',
           'osd', 'pielou_e', 'robbins', 'shannon', 'simpson', 'simpson_e',
           'singles', 'strong', 'lladser_pe', 'lladser_ci']
