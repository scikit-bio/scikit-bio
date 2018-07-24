# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn

import pandas as pd
import numpy as np
from scipy.linalg import eigh

from skbio.stats.distance import DistanceMatrix
from skbio.util._decorator import experimental
from ._ordination_results import OrdinationResults
from ._utils import e_matrix, f_matrix, scale

# - In cogent, after computing eigenvalues/vectors, the imaginary part
#   is dropped, if any. We know for a fact that the eigenvalues are
#   real, so that's not necessary, but eigenvectors can in principle
#   be complex (see for example
#   http://math.stackexchange.com/a/47807/109129 for details) and in
#   that case dropping the imaginary part means they'd no longer be
#   so, so I'm not doing that.


@experimental(as_of="0.4.0")
def pcoa(distance_matrix):
    r"""Perform Principal Coordinate Analysis.

    Principal Coordinate Analysis (PCoA) is a method similar to PCA
    that works from distance matrices, and so it can be used with
    ecologically meaningful distances like UniFrac for bacteria.

    In ecology, the euclidean distance preserved by Principal
    Component Analysis (PCA) is often not a good choice because it
    deals poorly with double zeros (Species have unimodal
    distributions along environmental gradients, so if a species is
    absent from two sites at the same site, it can't be known if an
    environmental variable is too high in one of them and too low in
    the other, or too low in both, etc. On the other hand, if an
    species is present in two sites, that means that the sites are
    similar.).

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        A distance matrix.

    Returns
    -------
    OrdinationResults
        Object that stores the PCoA results, including eigenvalues, the
        proportion explained by each of them, and transformed sample
        coordinates.

    See Also
    --------
    OrdinationResults

    Notes
    -----
    It is sometimes known as metric multidimensional scaling or
    classical scaling.

    .. note::

       If the distance is not euclidean (for example if it is a
       semimetric and the triangle inequality doesn't hold),
       negative eigenvalues can appear. There are different ways
       to deal with that problem (see Legendre & Legendre 1998, \S
       9.2.3), but none are currently implemented here.

       However, a warning is raised whenever negative eigenvalues
       appear, allowing the user to decide if they can be safely
       ignored.
    """
    distance_matrix = DistanceMatrix(distance_matrix)

    E_matrix = e_matrix(distance_matrix.data)

    # If the used distance was euclidean, pairwise distances
    # needn't be computed from the data table Y because F_matrix =
    # Y.dot(Y.T) (if Y has been centred).
    F_matrix = f_matrix(E_matrix)

    # If the eigendecomposition ever became a bottleneck, it could
    # be replaced with an iterative version that computes the
    # largest k eigenvectors.
    eigvals, eigvecs = eigh(F_matrix)

    # eigvals might not be ordered, so we order them (at least one
    # is zero). cogent makes eigenvalues positive by taking the
    # abs value, but that doesn't seem to be an approach accepted
    # by L&L to deal with negative eigenvalues. We raise a warning
    # in that case. First, we make values close to 0 equal to 0.
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0
    if np.any(eigvals < 0):
        warn(
            "The result contains negative eigenvalues."
            " Please compare their magnitude with the magnitude of some"
            " of the largest positive eigenvalues. If the negative ones"
            " are smaller, it's probably safe to ignore them, but if they"
            " are large in magnitude, the results won't be useful. See the"
            " Notes section for more details. The smallest eigenvalue is"
            " {0} and the largest is {1}.".format(eigvals.min(),
                                                  eigvals.max()),
            RuntimeWarning
            )
    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]

    # Scale eigenvalues to have lenght = sqrt(eigenvalue). This
    # works because np.linalg.eigh returns normalized
    # eigenvectors. Each row contains the coordinates of the
    # objects in the space of principal coordinates. Note that at
    # least one eigenvalue is zero because only n-1 axes are
    # needed to represent n points in an euclidean space.

    # If we return only the coordinates that make sense (i.e., that have a
    # corresponding positive eigenvalue), then Jackknifed Beta Diversity
    # won't work as it expects all the OrdinationResults to have the same
    # number of coordinates. In order to solve this issue, we return the
    # coordinates that have a negative eigenvalue as 0
    num_positive = (eigvals >= 0).sum()
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)

    coordinates = eigvecs * np.sqrt(eigvals)
    proportion_explained = eigvals / eigvals.sum()

    axis_labels = ['PC%d' % i for i in range(1, eigvals.size + 1)]
    return OrdinationResults(
        short_method_name='PCoA',
        long_method_name='Principal Coordinate Analysis',
        eigvals=pd.Series(eigvals, index=axis_labels),
        samples=pd.DataFrame(coordinates, index=distance_matrix.ids,
                             columns=axis_labels),
        proportion_explained=pd.Series(proportion_explained,
                                       index=axis_labels))


@experimental(as_of="0.5.3")
def pcoa_biplot(ordination, y):
    """Compute the projection of descriptors into a PCoA matrix

    This implementation is as described in Chapter 9 of Legendre & Legendre,
    Numerical Ecology 3rd edition.

    Parameters
    ----------
    ordination: OrdinationResults
        The computed principal coordinates analysis of dimensions (n, c) where
        the matrix ``y`` will be projected onto.
    y: DataFrame
        Samples by features table of dimensions (n, m). These can be
        environmental features or abundance counts. This table should be
        normalized in cases of dimensionally heterogenous physical variables.

    Returns
    -------
    OrdinationResults
        The modified input object that includes projected features onto the
        ordination space in the ``features`` attribute.
    """

    # acknowledge that most saved ordinations lack a name, however if they have
    # a name, it should be PCoA
    if (ordination.short_method_name != '' and
       ordination.short_method_name != 'PCoA'):
        raise ValueError('This biplot computation can only be performed in a '
                         'PCoA matrix.')

    if set(y.index) != set(ordination.samples.index):
        raise ValueError('The eigenvectors and the descriptors must describe '
                         'the same samples.')

    eigvals = ordination.eigvals
    coordinates = ordination.samples
    N = coordinates.shape[0]

    # align the descriptors and eigenvectors in a sample-wise fashion
    y = y.reindex(coordinates.index)

    # S_pc from equation 9.44
    # Represents the covariance matrix between the features matrix and the
    # column-centered eigenvectors of the pcoa.
    spc = (1 / (N - 1)) * y.values.T.dot(scale(coordinates, ddof=1))

    # U_proj from equation 9.55, is the matrix of descriptors to be projected.
    #
    # Only get the power of non-zero values, otherwise this will raise a
    # divide by zero warning. There shouldn't be negative eigenvalues(?)
    Uproj = np.sqrt(N - 1) * spc.dot(np.diag(np.power(eigvals, -0.5,
                                                      where=eigvals > 0)))

    ordination.features = pd.DataFrame(data=Uproj,
                                       index=y.columns.copy(),
                                       columns=coordinates.columns.copy())

    return ordination
