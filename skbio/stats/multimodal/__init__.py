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

Multimodal analysis addresses a common challenge in biological studies: how to
integrate measurements from different molecular layers (e.g., 16S rRNA gene
sequencing, metabolomics, proteomics) collected from the same samples. By
learning shared latent representations, these methods can reveal relationships
between different types of biological entities.


MMvec (Microbe-Metabolite Vectors)
----------------------------------

The primary method in this module is MMvec, which learns a co-occurrence model
between two sets of features. While originally developed for microbe-metabolite
interactions, the methodology is generic and can be applied to **any two omics
modalities representable as compositional (count-based) data**, such as:

- Microbes and metabolites (the original use case)
- Microbes and host transcripts
- Proteins and metabolites
- Any pair of feature tables sharing the same samples

The core idea is to model the conditional probability of observing a feature
from one modality given the presence of a feature from another:

.. math::

    P(\text{feature}_j^{(2)} | \text{feature}_i^{(1)}) \propto
    \exp(U_i \cdot V_j + b_{U_i} + b_{V_j})

where :math:`U` and :math:`V` are low-dimensional embedding matrices and
:math:`b_U`, :math:`b_V` are bias terms. This neural network-inspired approach
enables:

- **Dimensionality reduction**: Features from both modalities are embedded in
  the same latent space, enabling visualization and downstream analysis.
- **Prediction**: Given a composition from one modality, predict expected
  profiles in the other modality.
- **Interpretability**: The learned conditional probabilities (ranks) indicate
  which features from one modality are associated with features from the other.

Example usage::

    >>> import numpy as np
    >>> import pandas as pd
    >>> from skbio.stats.multimodal import mmvec

    >>> # Microbe counts (samples x microbes)
    >>> microbes = pd.DataFrame(
    ...     np.random.randint(0, 100, size=(50, 10)),
    ...     columns=[f'OTU_{i}' for i in range(10)]
    ... )
    >>> # Metabolite intensities (samples x metabolites)
    >>> metabolites = pd.DataFrame(
    ...     np.random.randint(0, 100, size=(50, 20)),
    ...     columns=[f'metabolite_{i}' for i in range(20)]
    ... )

    >>> # Fit model
    >>> result = mmvec(microbes, metabolites, n_components=3)

    >>> # Access embeddings
    >>> result.microbe_embeddings.shape
    (10, 4)
    >>> result.metabolite_embeddings.shape
    (20, 4)

    >>> # Get conditional probabilities
    >>> probs = result.probabilities()

    >>> # Predict metabolites for new samples
    >>> predictions = result.predict(microbes.iloc[:5])

    >>> # Evaluate on held-out data
    >>> q2 = result.score(microbes.iloc[-10:], metabolites.iloc[-10:])


References
----------
.. [1] Morton, J.T., et al. "Learning representations of microbe-metabolite
   interactions." Nature Methods 16, 1306-1314 (2019).
   https://doi.org/10.1038/s41592-019-0616-3


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
