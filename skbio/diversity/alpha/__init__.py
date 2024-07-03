r"""Alpha diversity measures (:mod:`skbio.diversity.alpha`)
=======================================================

.. currentmodule:: skbio.diversity.alpha

This package provides implementations of various alpha diversity [1]_ metrics,
including measures of richness, diversity, evenness, dominance, and coverage.

Some functions generate confidence intervals (CIs). These functions have the
suffix ``_ci``.


Richness metrics
----------------

**Richness** [2]_ measures the number of species (taxa) in a community.

Due to incomplete sampling, the number of observed species (``sobs``) in a
sample is usually lower than the true number of species in the community.
Metrics have been proposed to estimate the latter based on the distribution
of observed species in the sample.

.. autosummary::
   :toctree:

   ace
   chao1
   chao1_ci
   doubles
   faith_pd
   margalef
   menhinick
   michaelis_menten_fit
   observed_features
   observed_otus
   osd
   singles
   sobs


Diversity metrics
-----------------

**Diversity** [3]_ measures the number and relative abundances of species
(taxa) in a community. It combines richness and evenness.

Some diversity metrics describe the effective number of species (a.k.a., true
diversity) -- the number of equally-abundant species that produce the same
diversity measurement.

.. autosummary::
   :toctree:

   brillouin_d
   enspie
   fisher_alpha
   hill
   inv_simpson
   kempton_taylor_q
   phydiv
   renyi
   shannon
   simpson
   tsallis


Evenness metrics
----------------

**Evenness** [4]_ (or equitability) measures the closeness of species (taxa) in a
community in terms of abundance (number of individuals within the species). The
calculation of evenness involves the relative abundances of species.

.. autosummary::
   :toctree:

   heip_e
   mcintosh_e
   pielou_e
   simpson_e


Dominance metrics
-----------------

**Dominance** [5]_ (or concentration) measures the degree that one or a few
most abundant species (taxa) represent the great majority of a community. It
can be considered as a measure of community unevenness.

It should be noted that higher dominance corresponds to lower biodiversity.

.. autosummary::
   :toctree:

   berger_parker_d
   dominance
   gini_index
   mcintosh_d
   simpson_d
   strong


Coverage metrics
----------------

**Coverage** [6]_ measures the proportion of individuals of a community that
have been observed (or unobserved) in a sample. It describes the completeness
of sampling.

.. autosummary::
   :toctree:

   esty_ci
   goods_coverage
   lladser_ci
   lladser_pe
   robbins


References
----------
.. [1] https://en.wikipedia.org/wiki/Alpha_diversity

.. [2] https://en.wikipedia.org/wiki/Species_richness

.. [3] https://en.wikipedia.org/wiki/Species_diversity

.. [4] https://en.wikipedia.org/wiki/Species_evenness

.. [5] https://en.wikipedia.org/wiki/Dominance_%28ecology%29

.. [6] Good, I. J. (1953). The population frequencies of species and the
   estimation of population parameters. Biometrika, 40(3-4), 237-264.

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._base import (
    berger_parker_d,
    brillouin_d,
    dominance,
    doubles,
    enspie,
    esty_ci,
    fisher_alpha,
    goods_coverage,
    heip_e,
    hill,
    inv_simpson,
    kempton_taylor_q,
    margalef,
    mcintosh_d,
    mcintosh_e,
    menhinick,
    michaelis_menten_fit,
    observed_features,
    observed_otus,
    osd,
    pielou_e,
    renyi,
    robbins,
    shannon,
    simpson,
    simpson_d,
    simpson_e,
    singles,
    sobs,
    strong,
    tsallis,
)
from ._ace import ace
from ._chao1 import chao1, chao1_ci
from ._gini import gini_index
from ._lladser import lladser_pe, lladser_ci
from ._pd import faith_pd, phydiv


__all__ = [
    "ace",
    "chao1",
    "chao1_ci",
    "berger_parker_d",
    "brillouin_d",
    "dominance",
    "doubles",
    "enspie",
    "esty_ci",
    "faith_pd",
    "fisher_alpha",
    "gini_index",
    "goods_coverage",
    "heip_e",
    "hill",
    "inv_simpson",
    "kempton_taylor_q",
    "lladser_pe",
    "lladser_ci",
    "margalef",
    "mcintosh_d",
    "mcintosh_e",
    "menhinick",
    "michaelis_menten_fit",
    "observed_features",
    "observed_otus",
    "osd",
    "phydiv",
    "pielou_e",
    "renyi",
    "robbins",
    "shannon",
    "simpson",
    "simpson_d",
    "simpson_e",
    "singles",
    "sobs",
    "strong",
    "tsallis",
]
