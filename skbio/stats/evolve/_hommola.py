# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from scipy.stats import pearsonr

from skbio import DistanceMatrix
from skbio.util._decorator import experimental


@experimental(as_of="0.4.0")
def hommola_cospeciation(host_dist, par_dist, interaction, permutations=999):
    """Perform Hommola et al (2009) host/parasite cospeciation test.

    This test for host/parasite cospeciation is as described in [1]_. This test
    is a modification of a Mantel test, expanded to accept the case where
    multiple hosts map to a single parasite (and vice versa).

    For a basic Mantel test, the distance matrices being compared must have the
    same number of values. To determine the significance of the correlations
    between distances in the two matrices, the correlation coefficient of those
    distances is calculated and compared to the correlation coefficients
    calculated from a set of matrices in which rows and columns have been
    permuted.

    In this test, rather than comparing host-host to parasite-parasite
    distances directly (requiring one host per parasite), the distances are
    compared for each interaction edge between host and parasite. Thus, a host
    interacting with two different parasites will be represented in two
    different edges, with the host-host distance for the comparison between
    those edges equal to zero, and the parasite-parasite distance equal to the
    distance between those two parasites. Like in the Mantel test, significance
    of the interaction is assessed by permutation, in this case permutation of
    the host-symbiont interaction links.

    Note that the null hypothesis being tested here is that the hosts and
    parasites have evolved independently of one another. The alternative to
    this is a somewhat weaker case than what is often implied with the term
    'cospeciation,' which is that each incidence of host speciation is
    recapitulated in an incidence of symbiont speciation (strict
    co-cladogenesis). Although there may be many factors that could contribute
    to non-independence of host and symbiont phylogenies, this loss of
    explanatory specificity comes with increased robustness to phylogenetic
    uncertainty. Thus, this test may be especially useful for cases where host
    and/or symbiont phylogenies are poorly resolved, or when simple correlation
    between host and symbiont evolution is of more interest than strict
    co-cladogenesis.

    This test requires pairwise distance matrices for hosts and symbionts, as
    well as an interaction matrix specifying links between hosts (in columns)
    and symbionts (in rows). This interaction matrix should have the same
    number of columns as the host distance matrix, and the same number of rows
    as the symbiont distance matrix. Interactions between hosts and symbionts
    should be indicated by values of ``1`` or ``True``, with non-interactions
    indicated by values of ``0`` or ``False``.

    Parameters
    ----------
    host_dist : 2-D array_like or DistanceMatrix
        Symmetric matrix of m x m pairwise distances between hosts.
    par_dist : 2-D array_like or DistanceMatrix
        Symmetric matrix of n x n pairwise distances between parasites.
    interaction : 2-D array_like, bool
        n x m binary matrix of parasite x host interactions. Order of hosts
        (columns) should be identical to order of hosts in `host_dist`, as
        should order of parasites (rows) be identical to order of parasites in
        `par_dist`.
    permutations : int, optional
        Number of permutations used to compute p-value. Must be greater than or
        equal to zero. If zero, statistical significance calculations will be
        skipped and the p-value will be ``np.nan``.

    Returns
    -------
    corr_coeff : float
        Pearson correlation coefficient of host : parasite association.
    p_value : float
        Significance of host : parasite association computed using
        `permutations` and a one-sided (greater) alternative hypothesis.
    perm_stats : 1-D numpy.ndarray, float
        Correlation coefficients observed using permuted host : parasite
        interactions. Length will be equal to the number of permutations used
        to compute p-value (see `permutations` parameter above).

    See Also
    --------
    skbio.stats.distance.mantel
    scipy.stats.pearsonr

    Notes
    -----
    It is assumed that the ordering of parasites in `par_dist` and hosts in
    `host_dist` are identical to their ordering in the rows and columns,
    respectively, of the interaction matrix.

    This code is loosely based on the original R code from [1]_.

    References
    ----------
    .. [1] Hommola K, Smith JE, Qiu Y, Gilks WR (2009) A Permutation Test of
       Host-Parasite Cospeciation. Molecular Biology and Evolution, 26,
       1457-1468.

    Examples
    --------
    >>> from skbio.stats.evolve import hommola_cospeciation

    Create arrays for host distances, parasite distances, and their
    interactions (data taken from example in [1]_):

    >>> hdist = [[0,3,8,8,9], [3,0,7,7,8], [8,7,0,6,7], [8,7,6,0,3],
    ...          [9,8,7,3,0]]
    >>> pdist = [[0,5,8,8,8], [5,0,7,7,7], [8,7,0,4,4], [8,7,4,0,2],
    ...          [8,7,4,2,0]]
    >>> interaction = [[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0],
    ...                [0,0,0,1,1]]

    Run the cospeciation test with 99 permutations. Note that the correlation
    coefficient for the observed values counts against the final reported
    p-value:

    >>> corr_coeff, p_value, perm_stats = hommola_cospeciation(
    ...     hdist, pdist, interaction, permutations=99)
    >>> print("%.3f" % corr_coeff)
    0.832

    In this case, the host distances have a fairly strong positive correlation
    with the symbiont distances. However, this may also reflect structure
    inherent in the phylogeny, and is not itself indicative of significance.

    >>> p_value <= 0.05
    True

    After permuting host : parasite interactions, we find that the observed
    correlation is indeed greater than we would expect by chance.

    """
    host_dist = DistanceMatrix(host_dist)
    par_dist = DistanceMatrix(par_dist)
    interaction = np.asarray(interaction, dtype=bool)

    num_hosts = host_dist.shape[0]
    num_pars = par_dist.shape[0]

    if num_hosts < 3 or num_pars < 3:
        raise ValueError("Distance matrices must be a minimum of 3x3 in size.")
    if num_hosts != interaction.shape[1]:
        raise ValueError("Number of interaction matrix columns must match "
                         "number of hosts in `host_dist`.")
    if num_pars != interaction.shape[0]:
        raise ValueError("Number of interaction matrix rows must match "
                         "number of parasites in `par_dist`.")
    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if interaction.sum() < 3:
        raise ValueError("Must have at least 3 host-parasite interactions in "
                         "`interaction`.")

    # shortcut to eliminate nested for-loops specifying pairwise interaction
    # partners as randomizeable indices
    pars, hosts = np.nonzero(interaction)
    pars_k_labels, pars_t_labels = _gen_lists(pars)
    hosts_k_labels, hosts_t_labels = _gen_lists(hosts)

    # get a vector of pairwise distances for each interaction edge
    x = _get_dist(hosts_k_labels, hosts_t_labels, host_dist.data,
                  np.arange(num_hosts))
    y = _get_dist(pars_k_labels, pars_t_labels, par_dist.data,
                  np.arange(num_pars))

    # calculate the observed correlation coefficient for these hosts/symbionts
    corr_coeff = pearsonr(x, y)[0]

    # now do permutatitons. initialize index lists of the appropriate size
    mp = np.arange(num_pars)
    mh = np.arange(num_hosts)

    # initialize list of shuffled correlation vals
    perm_stats = np.empty(permutations)

    if permutations == 0 or np.isnan(corr_coeff):
        p_value = np.nan
        perm_stats.fill(np.nan)
    else:
        for i in range(permutations):
            # generate a shuffled list of indexes for each permutation. this
            # effectively randomizes which host is associated with which
            # symbiont, but maintains the distribution of genetic distances
            np.random.shuffle(mp)
            np.random.shuffle(mh)

            # get pairwise distances in shuffled order
            y_p = _get_dist(pars_k_labels, pars_t_labels, par_dist.data, mp)
            x_p = _get_dist(hosts_k_labels, hosts_t_labels, host_dist.data, mh)

            # calculate shuffled correlation coefficient
            perm_stats[i] = pearsonr(x_p, y_p)[0]

        p_value = ((perm_stats >= corr_coeff).sum() + 1) / (permutations + 1)

    return corr_coeff, p_value, perm_stats


def _get_dist(k_labels, t_labels, dists, index):
    """Subset a distance matrix using a set of (randomizable) index labels.

    Parameters
    ----------
    k_labels : numpy.array
        index labels specifying row-wise member of pairwise interaction
    t_labels : numpy.array
        index labels specifying column-wise member of pairwise interaction
    dists : numpy.array
        pairwise distance matrix
    index : numpy.array of int
        permutable indices for changing order in pairwise distance matrix

    Returns
    -------
    vec : list of float
        List of distances associated with host:parasite edges.

    """
    return dists[index[k_labels], index[t_labels]]


def _gen_lists(labels):
    """Generate matched lists of row and column index labels.

    Shortcut function for generating matched lists of row and col index
    labels for the set of pairwise comparisons specified by the list of those
    indices recovered using ``np.nonzero(interaction)``.

    Reproduces values of iterated indices from the nested for-loops contained
    in ``get_dist`` function in original code from [1]_.

    Parameters
    ----------
    labels : numpy.array
        array containing the indices of nonzero elements in one dimension of an
        interaction matrix

    Returns
    -------
    k_labels : numpy.array
        index labels specifying row-wise member of pairwise interaction
    t_labels : numpy.array
        index labels specifying column-wise member of pairwise interaction

    References
    ----------
    .. [1] Hommola K, Smith JE, Qiu Y, Gilks WR (2009) A Permutation Test of
       Host-Parasite Cospeciation. Molecular Biology and Evolution, 26,
       1457-1468.

    """
    i_array, j_array = np.transpose(np.tri(len(labels)-1)).nonzero()
    j_array = j_array + 1
    return labels[i_array], labels[j_array]
