# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Simulation utilities for testing MMvec."""

import numpy as np
import pandas as pd

from skbio.stats.composition import clr_inv as softmax
from skbio.stats.composition import ilr_inv
from skbio.util import get_rng


def random_multimodal(
    n_microbes=20,
    n_metabolites=100,
    n_samples=100,
    latent_dim=3,
    low=-1,
    high=1,
    microbe_total=10,
    metabolite_total=100,
    uB=0,
    sigmaB=2,
    sigmaQ=0.1,
    uU=0,
    sigmaU=1,
    uV=0,
    sigmaV=1,
    seed=0,
):
    """Generate synthetic microbe-metabolite co-occurrence data.

    Parameters
    ----------
    n_microbes : int
        Number of microbial species to simulate.
    n_metabolites : int
        Number of metabolites to simulate.
    n_samples : int
        Number of samples to generate.
    latent_dim : int
        Number of latent dimensions for the embedding.
    low : float
        Lower bound of the sample gradient.
    high : float
        Upper bound of the sample gradient.
    microbe_total : int
        Total microbial counts per sample.
    metabolite_total : int
        Total metabolite counts per sample.
    uB : float
        Mean of regression coefficient distribution.
    sigmaB : float
        Standard deviation of regression coefficient distribution.
    sigmaQ : float
        Standard deviation of noise in microbe composition.
    uU : float
        Mean of microbe embedding distribution.
    sigmaU : float
        Standard deviation of microbe embedding distribution.
    uV : float
        Mean of metabolite embedding distribution.
    sigmaV : float
        Standard deviation of metabolite embedding distribution.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    microbe_counts : pd.DataFrame
        Count table of microbial abundances (n_samples, n_microbes).
    metabolite_counts : pd.DataFrame
        Count table of metabolite abundances (n_samples, n_metabolites).
    design : np.ndarray
        Design matrix with intercept and gradient (n_samples, 2).
    beta : np.ndarray
        Regression coefficients (2, n_microbes).
    U : np.ndarray
        True microbe embedding matrix (n_microbes, latent_dim).
    Ubias : np.ndarray
        True microbe bias vector (n_microbes, 1).
    V : np.ndarray
        True metabolite embedding matrix (latent_dim, n_metabolites - 1).
    Vbias : np.ndarray
        True metabolite bias vector (1, n_metabolites - 1).

    """
    rng = get_rng(seed)

    # Regression coefficients for gradient model
    beta = rng.normal(uB, sigmaB, size=(2, n_microbes))

    # Design matrix: intercept + gradient
    design = np.vstack(
        (np.ones(n_samples), np.linspace(low, high, n_samples))
    ).T

    # Generate microbe compositions from ILR with noise
    microbes = ilr_inv(
        rng.multivariate_normal(
            mean=np.zeros(n_microbes - 1),
            cov=np.diag([sigmaQ] * (n_microbes - 1)),
            size=n_samples,
        )
    )

    # Generate latent embeddings
    Umain = rng.normal(uU, sigmaU, size=(n_microbes, latent_dim))
    Vmain = rng.normal(uV, sigmaV, size=(latent_dim, n_metabolites - 1))

    Ubias = rng.normal(uU, sigmaU, size=(n_microbes, 1))
    Vbias = rng.normal(uV, sigmaV, size=(1, n_metabolites - 1))

    # Augmented matrices for computing probabilities
    U_ = np.hstack((np.ones((n_microbes, 1)), Ubias, Umain))
    V_ = np.vstack((Vbias, np.ones((1, n_metabolites - 1)), Vmain))

    # Compute conditional probabilities P(metabolite | microbe)
    phi = np.hstack((np.zeros((n_microbes, 1)), U_ @ V_))
    probs = softmax(phi)

    # Generate count data
    microbe_counts = np.zeros((n_samples, n_microbes))
    metabolite_counts = np.zeros((n_samples, n_metabolites))

    n1 = microbe_total
    n2 = metabolite_total // microbe_total

    for n in range(n_samples):
        # Draw microbe counts
        otu = rng.multinomial(n1, microbes[n, :])
        # For each microbe, draw metabolites conditional on that microbe
        for i in range(n_microbes):
            ms = rng.multinomial(otu[i] * n2, probs[i, :])
            metabolite_counts[n, :] += ms
        microbe_counts[n, :] += otu

    # Create DataFrames with meaningful IDs
    otu_ids = [f"OTU_{d}" for d in range(microbe_counts.shape[1])]
    ms_ids = [f"metabolite_{d}" for d in range(metabolite_counts.shape[1])]
    sample_ids = [f"sample_{d}" for d in range(metabolite_counts.shape[0])]

    microbe_counts = pd.DataFrame(
        microbe_counts, index=sample_ids, columns=otu_ids
    )
    metabolite_counts = pd.DataFrame(
        metabolite_counts, index=sample_ids, columns=ms_ids
    )

    return (
        microbe_counts,
        metabolite_counts,
        design,
        beta,
        Umain,
        Ubias,
        Vmain,
        Vbias,
    )
