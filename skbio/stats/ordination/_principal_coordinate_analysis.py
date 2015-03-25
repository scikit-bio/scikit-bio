# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from warnings import warn

import numpy as np
import scipy as sp

from skbio.stats.distance import DistanceMatrix
from ._base import Ordination, OrdinationResults

# - In cogent, after computing eigenvalues/vectors, the imaginary part
#   is dropped, if any. We know for a fact that the eigenvalues are
#   real, so that's not necessary, but eigenvectors can in principle
#   be complex (see for example
#   http://math.stackexchange.com/a/47807/109129 for details) and in
#   that case dropping the imaginary part means they'd no longer be
#   so, so I'm not doing that.


class PCoABase(Ordination):
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
    distance_matrix : DistanceMatrix
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
    # defined in subclass
    short_method_name = None
    long_method_name = None
    eig_methods = {}

    def __init__(self, distance_matrix, eig_f=None, **eig_kwargs):
        if isinstance(distance_matrix, DistanceMatrix):
            self.dm = np.asarray(distance_matrix.data, dtype=np.float64)
            self.ids = distance_matrix.ids
        else:
            raise TypeError("Input must be a DistanceMatrix.")

        self._eig_kwargs = eig_kwargs
        if eig_f is None:
            self._eig_f = np.linalg.eigh
            self._eig_name = 'eigh'
        else:
            if eig_f not in self._eig_methods:
                eig_methods = ', '.join(self._eig_methods)
                raise KeyError("Unknown method: %s. The available methods "
                               "are: %s" % (eig_f, eig_methods))
            else:
                self._eig_f = self._eig_methods[eig_f]
                self._eig_name = eig_f

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
        eigvals, eigvecs = self._eig_f(F_matrix, **self._eig_kwargs)

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
        self.eigvals = eigvals[idxs_descending]
        self.eigvecs = eigvecs[:, idxs_descending]

    def scores(self):
        """Compute coordinates in transformed space.

        Returns
        -------
        OrdinationResults
            Object that stores the computed eigenvalues, the
            proportion explained by each of them (per unit) and
            transformed coordinates, etc.

        See Also
        --------
        OrdinationResults
        """
        raise NotImplementedError("Must be defined in a subclass")

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


class PCoA(PCoABase):
    short_method_name = 'PCoA'
    long_method_name = 'Principal Coordinate Analysis'
    eig_methods = {'eigh': np.linalg.eigh,
                   'eigsh': sp.sparse.linalg.eigsh}

    def scores(self):
        """Compute coordinates in transformed space.

        Returns
        -------
        OrdinationResults
            Object that stores the computed eigenvalues, the
            proportion explained by each of them (per unit) and
            transformed coordinates, etc.

        See Also
        --------
        OrdinationResults
        """
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
        num_positive = (self.eigvals >= 0).sum()
        eigvecs = self.eigvecs
        eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
        eigvals = self.eigvals
        eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)

        coordinates = eigvecs * np.sqrt(eigvals)

        proportion_explained = eigvals / eigvals.sum()

        return OrdinationResults(eigvals=eigvals, site=coordinates,
                                 proportion_explained=proportion_explained,
                                 site_ids=self.ids)
