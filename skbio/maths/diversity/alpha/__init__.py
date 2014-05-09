#!/usr/bin/env python
"""
Alpha diversity measures (:mod:`skbio.maths.diversity.alpha`)
=============================================================

.. currentmodule:: skbio.maths.diversity.alpha

This package provides implementations of various alpha diversity measures,
including measures of richness, dominance, and evenness. Some functions also
have the ability to generate confidence intervals (CIs).

All alpha diversity measures accept a vector of counts within a single sample,
where each count is the number of individuals (e.g., biological sequences) that
were observed for a particular "species". We use the term "species" here very
loosely, as these could be counts of any type of feature/observation
(e.g., OTUs). We'll refer to this vector as the *counts vector* or simply
*counts* throughout the documentation.

The counts vector must be one-dimensional and contain integers representing the
number of individuals seen (or *counted*) for a particular species. Negative
values are not allowed; the counts vector may only contain integers greater
than or equal to zero.

The counts vector is `array_like`: anything that can be converted into a 1-D
numpy array is acceptable input. For example, you can provide a numpy array or
a native Python list and the results should be identical.

If the input to an alpha diversity measure does not meet the above
requirements, the function will raise either a ``ValueError`` or a
``TypeError``, depending on the condition that is violated.

.. note:: There are different ways that samples are represented in the
   ecological literature and in related software. The alpha diversity measures
   provided here *always* assume that the input contains abundance data: each
   count represents the number of individuals seen for a particular species in
   the sample. For example, if you have two species, where 3 individuals were
   observed in one of the species and only a single individual was observed in
   the other, you could represent this data in the following forms (among
   others):

   As a vector of counts. This is the expected type of input for the alpha
   diversity measures in this module. There are 3 individuals for the species
   at index 0, and 1 individual for the species at index 1:

   >>> counts = [3, 1]

   As a vector of indices. The species at index 0 appears 3 times, while the
   species at index 1 appears 1 time:

   >>> indices = [0, 0, 0, 1]

   As a vector of frequencies. We have 1 species that is a singleton and 1
   species that is a tripleton. We do not have any 0-tons or doubletons:

   >>> frequencies = [0, 1, 0, 1]

   Always use the first representation (a counts vector) with this module.

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
   singles
   strong

Examples
--------

Import a few of the alpha diversity measures we'll be using in the examples:

>>> import numpy as np

Assume we have the following abundance data for a sample, represented as a
counts vector:

>>> counts = [1, 0, 0, 4, 1, 2, 3, 0]

We can count the number of "species":

>>> observed_species(counts)
5

Note that the species with counts of zero are ignored.

In the previous example, we provided a Python list as input. We can also
provide other types of input that are `array_like`:

>>> observed_species((1, 0, 0, 4, 1, 2, 3, 0)) # tuple
5
>>> observed_species(np.array([1, 0, 0, 4, 1, 2, 3, 0])) # numpy array
5

All of the alpha diversity measures work in this manner.

Let's see how many singletons and doubletons there are in the sample:

>>> singles(counts)
2
>>> doubles(counts)
1

Let's calculate Menhinick's richness index:

>>> menhinick(counts)
1.507556722888818

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
    osd, robbins, shannon, simpson, simpson_e, singles, strong)
from .gini import gini_index
from .lladser import lladser_pe, lladser_ci

__all__ = [m for m in dir() if not m.startswith('_')]

from numpy.testing import Tester
test = Tester().test
