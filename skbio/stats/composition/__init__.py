r"""Composition Statistics (:mod:`skbio.stats.composition`)
=======================================================

.. currentmodule:: skbio.stats.composition

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


Differential abundance
----------------------

Statistical tests for the differential abundance (DA) of components among groups of
compositions.

.. autosummary::
   :toctree:

   ancom
   ancombc
   dirmult_ttest
   dirmult_lme
   struc_zero

.. note::
   Differential abundance tests will be moved to a separate module ``differential`` in
   the next release of scikit-bio. The current location will be kept as an alias.


Arithmetic operations
---------------------

Manipulate compositional data within the Aitchison space.

.. autosummary::
   :toctree:

   centralize
   closure
   inner
   perturb
   perturb_inv
   power


Log-ratio transformation
------------------------

Convert compositional data into log-ratio space to enable subsequent comparison
and statistical analysis.

.. autosummary::
   :toctree:

   alr
   alr_inv
   clr
   clr_inv
   rclr
   ilr
   ilr_inv

.. note::
   Arithmetic operations and log-ratio transformations support array formats compliant
   with the `Python array API standard <https://data-apis.org/array-api/latest/>`_
   without transition through NumPy. For example, they can directly consume and return
   GPU-resident PyTorch tensors.


Correlation analysis
--------------------

Measure the pairwise relationships of compositional data.

.. autosummary::
   :toctree:

   vlr
   pairwise_vlr


Zero handling
-------------

Replace zero values in compositional data with positive values, which is
necessary prior to logarithmic operations.

.. autosummary::
   :toctree:

   multi_replace


Basis construction
------------------

Generate basis vectors for compositional data via hierarchical partitioning, to
allow for decomposition and transformation, such as ilr transform.

.. autosummary::
   :toctree:

   sbp_basis
   tree_basis


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

from ._base import (
    centralize,
    closure,
    inner,
    perturb,
    perturb_inv,
    power,
    alr,
    alr_inv,
    clr,
    clr_inv,
    rclr,
    ilr,
    ilr_inv,
    vlr,
    pairwise_vlr,
    multi_replace,
    sbp_basis,
    tree_basis,
)
from ._ancom import ancom
from ._ancombc import ancombc, struc_zero
from ._dirmult import dirmult_ttest, dirmult_lme

__all__ = [
    "centralize",
    "closure",
    "inner",
    "perturb",
    "perturb_inv",
    "power",
    "alr",
    "alr_inv",
    "clr",
    "clr_inv",
    "rclr",
    "ilr",
    "ilr_inv",
    "vlr",
    "pairwise_vlr",
    "multi_replace",
    "sbp_basis",
    "tree_basis",
    "ancom",
    "ancombc",
    "struc_zero",
    "dirmult_ttest",
    "dirmult_lme",
]
