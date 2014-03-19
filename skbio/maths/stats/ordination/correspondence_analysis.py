#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from .base import Ordination, OrdinationResults
from .utils import svd_rank


class CA(Ordination):
    r"""Compute correspondence analysis, a multivariate statistical
    technique for ordination.

    In general, rows in the data table will correspond to sites and
    columns to species, but the method is symmetric. In order to
    measure the correspondence between rows and columns, the
    :math:`\chi^2` distance is used, and those distances are preserved
    in the transformed space. The :math:`\chi^2` distance doesn't take
    double zeros into account, and so it is expected to produce better
    ordination that PCA when the data has lots of zero values.

    It is related to Principal Component Analysis (PCA) but it should
    be preferred in the case of steep or long gradients, that is, when
    there are many zeros in the input data matrix.

    Parameters
    ----------
    X : array_like
        Contingency table. It can be applied to different kinds of
        data tables but data must be non-negative and dimensionally
        homogeneous (quantitative or binary).

    Notes
    -----
    The algorithm is based on Legendre & Legendre (1998) 9.4.1. and is
    expected to give the same results as ``cca(X)`` in R's package
    vegan.

    See Also
    --------
    CCA
    """
    short_method_name = 'CA'
    long_method_name = 'Canonical Analysis'

    def __init__(self, X):
        self.X = np.asarray(X, dtype=np.float64)
        self._ca()

    def _ca(self):
        X = self.X
        r, c = X.shape

        if X.min() < 0:
            raise ValueError("Input matrix elements must be non-negative.")

        # Step 1 (similar to Pearson chi-square statistic)
        grand_total = X.sum()
        Q = X / grand_total

        column_marginals = Q.sum(axis=0)
        row_marginals = Q.sum(axis=1)
        # Let's store them since they're needed to compute scores
        self.column_marginals = column_marginals
        self.row_marginals = row_marginals

        # Formula 9.32 in Lagrange & Lagrange (1998). Notice that it's
        # an scaled version of the contribution of each cell towards
        # Pearson chi-square statistic.
        expected = np.outer(row_marginals, column_marginals)
        Q_bar = (Q - expected) / np.sqrt(expected)  # Eq. 9.32

        # Step 2 (Singular Value Decomposition)
        U_hat, W, Ut = np.linalg.svd(Q_bar, full_matrices=False)
        # Due to the centering, there are at most min(r, c) - 1 non-zero
        # eigenvalues (which are all positive)
        rank = svd_rank(Q_bar.shape, W)
        assert rank <= min(r, c) - 1
        self.U_hat = U_hat[:, :rank]
        self.W = W[:rank]
        self.U = Ut[:rank].T

    def scores(self, scaling):
        r"""Compute site and species scores for different scalings.

        Parameters
        ----------
        scaling : int

            For a more detailed explanation of the interpretation, check
            Legendre & Legendre 1998, section 9.4.3. The notes that
            follow are quick recommendations.

            Scaling type 1 maintains :math:`\chi^2` distances between
            rows (sites): in the transformed space, the euclidean
            distances between rows are equal to the :math:`\chi^2`
            distances between rows in the original space. It should be
            used when studying the ordination of sites. Rows (sites)
            that are near a column (species) have high contributions
            from it.

            Scaling type 2 preserves :math:`\chi^2` distances between
            columns (species), so euclidean distance between columns
            after transformation is equal to :math:`\chi^2` distance
            between columns in the original space. It is best used
            when we are interested in the ordination of species. A
            column (species) that is next to a row (site) means that
            it is more abundant there.

            Other types of scalings are currently not implemented, as
            they're less used by ecologists (Legendre & Legendre 1998,
            p. 456).

            In general, species appearing far from the center of the
            biplot and far from its edges will probably exhibit better
            relationships than species either in the center (may be
            multimodal species, not related to the shown ordination
            axes...) or the edges (sparse species...).
        """

        if scaling not in {1, 2}:
            raise NotImplementedError(
                "Scaling {0} not implemented.".format(scaling))
        # Both scalings are a bit intertwined, so we'll compute both and
        # then choose
        V = self.column_marginals[:, None]**-0.5 * self.U
        V_hat = self.row_marginals[:, None]**-0.5 * self.U_hat
        F = V_hat * self.W
        # According to Formula 9.43, this should hold
        #assert np.allclose(F, (row_marginals**-1)[:, None] * Q.dot(V))
        # but it doesn't (notice that W**2==Lambda):
        # (9.43a) F = V_hat W = D(p_i+)^{-1/2} U_hat W
        #           = D(p_i+)^{-1/2} Q_bar U W^{-1} W  (substituting 9.38)
        #           = D(p_i+)^{-1/2} Q_bar U
        # (9.43b) F = D(p_i+)^{-1} Q V
        #           = D(p_i+)^{-1} Q D(p_+j)^{-1/2} U  (substituting 9.41)
        #           = D(p_i+)^{-1/2} D(p_i+)^{-1/2} Q D(p_+j)^{-1/2} U
        #           = D(p_i+)^{-1/2} Q_tilde U         (using 9.40)
        # It holds if we replace Q in 9.43b with Q after centering, ie
        #assert np.allclose(
        #    F,
        #    (row_marginals**-1)[:, None] * (Q - expected).dot(V))
        # Comparing results with vegan and the examples in the book, 9.43a
        # is the right one. The same issue happens in 9.44, where also
        # 9.44a is the one that matches vegan's output.
        # (9.44a) F_hat = V W = D(p_+j)^{-1/2} U W
        #               = D(p_+j)^{-1/2} Q_bar' U_hat W^{-1} W (using 9.39)
        #               = D(p_+j)^{-1/2} Q_bar' U_hat
        # (9.44b) F_hat = D(p_+j)^{-1} Q' V_hat
        #               = D(p_+j)^{-1/2} Q_tilde' U_hat (using 9.40 and 9.42)
        F_hat = V * self.W

        # Eigenvalues
        eigvals = self.W**2

        # Species scores
        species_scores = [V, F_hat][scaling - 1]
        # Site scores (weighted averages of species scores)
        site_scores = [F, V_hat][scaling - 1]
        return OrdinationResults(eigvals=eigvals, species=species_scores,
                                 site=site_scores)
