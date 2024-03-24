r"""Composition Statistics (:mod:`skbio.composition`)
=================================================

.. currentmodule:: skbio.composition

This module provides functions for compositional data analysis.

Many omics datasets are inherently compositional -- meaning that they are best
interpreted as proportions or percentages rather than absolute counts.

Formally, sample :math:`x` is a composition if :math:`\sum_{i=0}^D x_{i} = c`
and :math:`x_{i} > 0`, :math:`1 \leq i \leq D` and :math:`c` is a real-valued
constant and there are :math:`D` components (features) for this composition.
In this module :math:`c=1`. Compositional data can be analyzed using
**Aitchison geometry** [1]_.

However, in this framework, standard real Euclidean operations such as addition
and multiplication no longer apply. Only operations such as perturbation and
power can be used to manipulate this data.

This module allows two styles of manipulation of compositional data.
Compositional data can be analyzed using perturbation and power operations,
which can be useful for simulation studies. The alternative strategy is to
transform compositional data into the real space. Right now, the centre log
ratio transform (clr) and the isometric log ratio transform (ilr) [2]_ can be
used to accomplish this. This transform can be useful for performing standard
statistical methods such as parametric hypothesis testing, regression and more.

The major caveat of using this framework is dealing with zeros. In Aitchison
geometry, only compositions with non-zero components can be considered.
The multiplicative replacement technique [3]_ can be used to substitute these
zeros with small pseudocounts without introducing major distortions to the
data.

Functions
---------

.. autosummary::
   :toctree: generated/

   closure
   multi_replace
   multiplicative_replacement
   perturb
   perturb_inv
   power
   inner
   clr
   clr_inv
   ilr
   ilr_inv
   alr
   alr_inv
   centralize
   vlr
   pairwise_vlr
   tree_basis
   ancom
   sbp_basis
   dirmult_ttest

References
----------
.. [1] V. Pawlowsky-Glahn, J. J. Egozcue, R. Tolosana-Delgado (2015),
   Modeling and Analysis of Compositional Data, Wiley, Chichester, UK

.. [2] J. J. Egozcue.,  "Isometric Logratio Transformations for
   Compositional Data Analysis" Mathematical Geology, 35.3 (2003)

.. [3] J. A. Martin-Fernandez,  "Dealing With Zeros and Missing Values in
   Compositional Data Sets Using Nonparametric Imputation",
   Mathematical Geology, 35.3 (2003)

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._composition import (
    closure,
    multi_replace,
    multiplicative_replacement,
    perturb,
    perturb_inv,
    power,
    inner,
    clr,
    clr_inv,
    ilr,
    ilr_inv,
    alr,
    alr_inv,
    centralize,
    vlr,
    pairwise_vlr,
    tree_basis,
    ancom,
    sbp_basis,
    dirmult_ttest,
)


__all__ = [
    "closure",
    "multi_replace",
    "multiplicative_replacement",
    "perturb",
    "perturb_inv",
    "power",
    "inner",
    "clr",
    "clr_inv",
    "ilr",
    "ilr_inv",
    "alr",
    "alr_inv",
    "centralize",
    "vlr",
    "pairwise_vlr",
    "tree_basis",
    "ancom",
    "sbp_basis",
    "dirmult_ttest",
]
