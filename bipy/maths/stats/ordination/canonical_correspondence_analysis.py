#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from .base import Ordination, OrdinationResults
from .utils import corr, svd_rank, scale


class CCA(Ordination):
    short_method_name = 'CCA'
    long_method_name = 'Canonical Correspondence Analysis'

    def __init__(self, Y, X):
        r"""Compute constrained (also known as canonical) correspondence
        analysis.

        Canonical (or constrained) correspondence analysis is a
        multivariate ordination technique.

        In exploratory data analysis, ordination (or multivariate
        gradient analysis) complements clustering by arranging objects
        (species, samples...) along gradients so that similar ones are
        closer and dissimilar ones are further. There's a good overview
        of the available techniques in
        http://ordination.okstate.edu/overview.htm.

        Parameters
        ----------
        Y : array_like
            Community data matrix of shape (n, m): a contingency table
            for m species at n sites.
        X : array_like
            Constraining matrix of shape (n, q): q quantitative
            environmental variables at n sites.

        Notes
        -----
        There are a couple of techniques to search for multivariate
        relationships between two datasets with very similar names, so
        this can get confusing:

        - Canonical correlation analysis is a statistical tool that,
          given two vectors of random variables, finds linear
          combinations that have maximum correlation with each other. In
          some sense, assumes linear responses of *species* to
          *environmental variables*. This is not implemented here.

        - Canonical (or constrained) correspondence analysis appeared in
          community ecology [1]_ and relates community composition to
          the variation in the environment (or in other factors). It
          works from data on abundances or counts of individuals and
          environmental variables, and outputs ordination axes that
          maximize niche separation among species. It is better suited
          to extract the niches of taxa than linear multivariate methods
          like canonical correlation analysis because it assumes
          unimodal response curves (habitat preferences are often
          unimodal functions of habitat variables [2]_).

        References
        ----------

        .. [1] Cajo J. F. Ter Braak, "Canonical Correspondence
            Analysis: A New Eigenvector Technique for Multivariate
            Direct Gradient Analysis", Ecology 67.5 (1986),
            pp. 1167-1179.

        .. [2] Cajo J.F. Braak and Piet F.M. Verdonschot, "Canonical
            correspondence analysis and related multivariate methods
            in aquatic ecology", Aquatic Sciences 57.3 (1995),
            pp. 255-289.

        """

        self.Y = np.asarray(Y, dtype=np.float64)
        self.X = np.asarray(X, dtype=np.float64)
        self._cca()

    def _cca(self):
        X, Y = self.X, self.Y
        if X.shape[0] != Y.shape[0]:
            raise ValueError("Contingency and environmental tables must have"
                             " the same number of rows (sites). X has {0}"
                             " rows but Y has {1}.".format(X.shape[0],
                                                           Y.shape[0]))
        if Y.min() < 0:
            raise ValueError("Contingency table must be nonnegative")
        row_max = Y.max(axis=1)
        if np.any(row_max <= 0):
            # Or else the lstsq call to compute Y_hat breaks
            raise ValueError("Contingency table cannot contain row of only 0s")

        # Step 1 (similar to Pearson chi-square statistic)
        grand_total = Y.sum()
        Q = Y / grand_total  # Relative frequencies of X (contingency table)

        # Species and site weights (marginal totals)
        column_marginals = Q.sum(axis=0)
        row_marginals = Q.sum(axis=1)

        # Formula 9.32 in Lagrange & Lagrange (1998). Notice that it's an
        # scaled version of the contribution of each cell towards Pearson
        # chi-square statistic.
        expected = np.outer(row_marginals, column_marginals)
        Q_bar = (Q - expected) / np.sqrt(expected)

        # Step 2. Standardize columns of Y with respect to site weights,
        # using the maximum likelyhood variance estimator (Legendre &
        # Legendre 1998, p. 595)
        X = scale(X, weights=row_marginals, ddof=0)

        # Step 3. Weighted multiple regression.
        X_weighted = row_marginals[:, None]**0.5 * X
        B, _, rank_lstsq, _ = np.linalg.lstsq(X_weighted, Q_bar)
        Y_hat = X_weighted.dot(B)
        Y_res = Q_bar - Y_hat

        # Step 4. Eigenvalue decomposition
        u, s, vt = np.linalg.svd(Y_hat, full_matrices=False)
        rank = svd_rank(Y_hat.shape, s)
        s = s[:rank]
        u = u[:, :rank]
        vt = vt[:rank]
        U = vt.T

        # Step 5. Eq. 9.38
        U_hat = Q_bar.dot(U) * s**-1

        # Residuals analysis
        u_res, s_res, vt_res = np.linalg.svd(Y_res, full_matrices=False)
        rank = svd_rank(Y_res.shape, s_res)
        s_res = s_res[:rank]
        u_res = u_res[:, :rank]
        vt_res = vt_res[:rank]

        U_res = vt_res.T
        U_hat_res = Q_bar.dot(U_res) * s_res**-1

        # Storing values needed to compute scores
        iter_ = (('column_marginals', column_marginals),
                 ('row_marginals', row_marginals),
                 ('U', U),
                 ('U_res', U_res),
                 ('U_hat', U_hat),
                 ('U_hat_res', U_hat_res),
                 ('u', u), ('Y_hat', Y_hat),
                 ('s', s), ('s_res', s_res),
                 ('X_weighted', X_weighted[:, :rank_lstsq]))
        for val_name, val in iter_:
            setattr(self, val_name, val)

        self.eigenvalues = np.r_[s, s_res]**2

    def scores(self, scaling):
        if scaling not in {1, 2}:
            raise NotImplementedError(
                "Scaling {0} not implemented.".format(scaling))
        # In this case scores are also a bit intertwined, so we'll
        # almost compute them both and then choose.

        ## Scalings (p. 596 L&L 1998)
        # Species scores, scaling 1
        V = (self.column_marginals**-0.5)[:, None] * self.U

        # Site scores, scaling 2
        V_hat = (self.row_marginals**-0.5)[:, None] * self.U_hat

        # Site scores, scaling 1
        F = V_hat * self.s

        # Species scores, scaling 2
        F_hat = V * self.s

        # Site scores which are linear combinations of environmental
        # variables
        Z_scaling1 = ((self.row_marginals**-0.5)[:, None] *
                      self.Y_hat.dot(self.U))
        Z_scaling2 = Z_scaling1 * self.s**-1

        # Species residual scores, scaling 1
        V_res = (self.column_marginals**-0.5)[:, None] * self.U_res

        # Site residual scores, scaling 2
        V_hat_res = (self.row_marginals**-0.5)[:, None] * self.U_hat_res

        # Site residual scores, scaling 1
        F_res = V_hat_res * self.s_res

        # Species residual scores, scaling 2
        F_hat_res = V_res * self.s_res

        eigvals = self.eigenvalues
        species_scores = np.hstack([(V, V_res),
                                    (F_hat, F_hat_res)][scaling - 1])
        site_scores = np.hstack([(F, F_res),
                                 (V_hat, V_hat_res)][scaling - 1])
        site_constraints = np.hstack([(Z_scaling1, F_res),
                                      (Z_scaling2, V_hat_res)][scaling - 1])
        biplot_scores = corr(self.X_weighted, self.u)
        return OrdinationResults(eigvals=eigvals,
                                 species=species_scores,
                                 site=site_scores,
                                 biplot=biplot_scores,
                                 site_constraints=site_constraints)
