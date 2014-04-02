#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import warnings

import numpy as np

from .base import Ordination, OrdinationResults
from skbio.core.distance import SymmetricDistanceMatrix

# - In cogent, after computing eigenvalues/vectors, the imaginary part
#   is dropped, if any. We know for a fact that the eigenvalues are
#   real, so that's not necessary, but eigenvectors can in principle
#   be complex (see for example
#   http://math.stackexchange.com/a/47807/109129 for details) and in
#   that case dropping the imaginary part means they'd no longer be
#   so, so I'm not doing that.

# - The rest of the ordination files works from a data table (sites x
#   species), but PCoA works from a distance matrix, so the ordination
#   results from the former (i.e., site scores, species scores, etc)
#   don't map very well to PCoA ordination results (i.e., "object"
#   scores). See also base.py


class PCoA(Ordination):
    r"""Perform Principal Coordinate Analysis.

    Principal Coordinate Analysis (PCoA) is a method similar to PCA
    that works from distance matrices, and so it can be used with
    ecologically meaningful distances like unifrac for bacteria.

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
    ==========
    distance_matrix : SymmetricDistanceMatrix
        A distance matrix.

    Notes
    =====
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
    short_method_name = 'PCoA'
    long_method_name = 'Principal Coordinate Analysis'

    def __init__(self, distance_matrix):
        if isinstance(distance_matrix, SymmetricDistanceMatrix):
            self.dm = np.asarray(distance_matrix.data, dtype=np.float64)
        else:
            raise TypeError("Input must be a SymmetricDistanceMatrix.")
        self._pcoa()

    def _pcoa(self):
        E_matrix = self._E_matrix(self.dm)

        # If the used distance was euclidean, pairwise distances
        # needn't be computed from the data table Y because F_matrix =
        # Y.dot(Y.T) (if Y has been centred).
        F_matrix = self._F_matrix(E_matrix)

        # If the eigendecomposition ever became a bottleneck, it could
        # be replaced with an iterative version that computes the
        # largest k eigenvectors.
        eigvals, eigvecs = np.linalg.eigh(F_matrix)

        # eigvals might not be ordered, so we order them (at least one
        # is zero). cogent makes eigenvalues positive by taking the
        # abs value, but that doesn't seem to be an approach accepted
        # by L&L to deal with negative eigenvalues. We raise a warning
        # in that case. First, we make values close to 0 equal to 0.
        negative_close_to_zero = np.isclose(eigvals, 0)
        eigvals[negative_close_to_zero] = 0
        if np.any(eigvals < 0):
            warnings.warn(
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
        self.eigvals = eigvals[idxs_descending]
        self.eigvecs = eigvecs[:, idxs_descending]

    def scores(self):
        # Scale eigenvalues to have lenght = sqrt(eigenvalue). This
        # works because np.linalg.eigh returns normalized
        # eigenvectors. Each row contains the coordinates of the
        # objects in the space of principal coordinates. Note that at
        # least one eigenvalue is zero because only n-1 axes are
        # needed to represent n points in an euclidean space.

        # We only return coordinates that make sense (i.e., that have
        # a corresponding positive eigenvalue)
        num_positive = (self.eigvals > 0).sum()
        eigvecs = self.eigvecs[:, :num_positive]
        eigvals = self.eigvals[:num_positive]

        coordinates = eigvecs * np.sqrt(eigvals)

        perc_expl = (eigvals / eigvals.sum()) * 100

        # TODO: Improve OrdinationResults to better cope with PCoA
        return OrdinationResults(eigvals=eigvals, species=coordinates,
                                 perc_expl=perc_expl)

    @staticmethod
    def _E_matrix(distance_matrix):
        """Compute E matrix from a distance matrix.

        Squares and divides by -2 the input elementwise. Eq. 9.20 in
        Legendre & Legendre 1998."""
        return distance_matrix * distance_matrix / -2

    @staticmethod
    def _F_matrix(E_matrix):
        """Compute F matrix from E matrix.

        Centring step: for each element, the mean of the corresponding
        row and column are substracted, and the mean of the whole
        matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
        row_means = E_matrix.mean(axis=1, keepdims=True)
        col_means = E_matrix.mean(axis=0, keepdims=True)
        matrix_mean = E_matrix.mean()
        return E_matrix - row_means - col_means + matrix_mean
