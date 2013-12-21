#! /usr/bin/env python
from __future__ import print_function, absolute_import
from collections import namedtuple
try:  # py2 compatibility
    from itertools import izip as zip, imap as map
except ImportError:
    pass
from itertools import count
import numpy as np
import matplotlib.pyplot as plt

from .utils import corr, svd_rank, scale


class OrdinationResults(namedtuple('OrdinationResults',
                                   ('eigvals', 'species', 'site', 'biplot',
                                    'site_constraints'))):
    __slots__ = ()  # To avoid creating a dict, as a namedtuple
                    # doesn't have it

    def __new__(cls, eigvals, species, site, biplot=None, site_constraints=None):
        return super(OrdinationResults, cls).__new__(cls, eigvals, species, site,
                                                     biplot, site_constraints)


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'

    def biplot(self, scaling=2, choices=[0, 1]):
        choices = list(choices)

        def scatter_with_names(ax, data, names, marker, color, label):
            """Helper to plot scatter points and their names."""
            xdata, ydata = data[:, choices].T
            ax.scatter(xdata, ydata, marker=marker, color=color, label=label)
            for namei, xi, yi in zip(names, xdata, ydata):
                ax.annotate(namei, (xi, yi), xytext=(0, 5),
                            textcoords='offset points', ha='center')

        def make_names(prefix):
            prefix += "{n}"
            return map(lambda i: prefix.format(n=i), count())

        # Extract scores
        scores = self.scores(scaling)
        eigvals = scores.eigvals
        species_scores = scores.species
        site_scores = scores.site
        biplot_scores = scores.biplot  # Can be None
        # Actual plotting
        fig, ax = plt.subplots()
        plt.axvline(color='k')
        plt.axhline(color='k')
        scatter_with_names(ax, site_scores, make_names("R"), marker='o',
                           color='b', label='Sites (rows)')
        scatter_with_names(ax, species_scores, make_names("C"),
                           marker='s', color='r', label='Species (columns)')
        if biplot_scores is not None:
            for (xi, yi) in biplot_scores[:, choices]:
                ax.arrow(0, 0, xi, yi, head_width=0.1, facecolor='none')

        for i, set_label in enumerate([ax.set_xlabel, ax.set_ylabel]):
            set_label('{ax_name}{ax_n} - {proportion:.3%}'.format(
                ax_name=self.__class__.short_method_name,
                ax_n=choices[i],
                proportion=(eigvals / np.sum(eigvals))[choices[i]]))
        ax.legend(loc='best').draggable()
        ax.set_title('{0} (scaling = {1})'.format(
            self.__class__.long_method_name, scaling))
        ax.set_aspect('equal')
        plt.show()


class CA(Ordination):
    short_method_name = 'CA'
    long_method_name = 'Canonical Analysis'

    def __init__(self, X):
        r"""Compute correspondence analysis, a multivariate statistical
        technique for ordination.

        In general, rows in the data table will correspond to sites and
        columns to species, but the method is symmetric. In order to
        measure the correspondence between rows and columns, the
        :math:`\chi^2` distance is used, and those distances are
        preserved in the transformed space.

        It is related to Principal Component Analysis (PCA) but it
        should be preferred in the case of steep or long gradients, that
        is, when there are many zeros in the input data matrix.

        Parameters
        ----------
        X : array_like
            Contingency table. Data must be non-negative and
            dimensionally homogeneous.

        Notes
        -----
        The algorithm is based on Legendre & Legendre (1998) 9.4.1. and
        is expected to give the same results as `cca(X)` in R's package
        vegan.
        """
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
            rows (sites) and should be used when checking the ordination
            of sites. Rows (sites) that are near a column (species) have
            high contributions from it.

            Scaling type 2 preserves :math:`\chi^2` distances between
            columns (species) and is best used when interested in the
            ordination of species. If a column (species) is next to a
            row (site) it means that it's more abundant there.

            Scaling type 3 is currently not implemented, as it's less
            used by ecologists (Legendre & Legendre 1998, p. 456).

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
        #assert np.allclose(F,
        #                   (row_marginals**-1)[:, None] * (Q - expected).dot(V))
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


class RDA(Ordination):
    short_method_name = 'RDA'
    long_method_name = 'Redundancy Analysis'

    def __init__(self, Y, X, scale_Y=False):
        r"""Compute redundancy analysis, a type of canonical analysis.

        It is related to PCA and multiple regression because the
        explained variables `Y` are fitted to the explanatory variables
        `X` and PCA is then performed on the fitted values. A similar
        process is performed on the residuals.

        RDA should be chosen if the studied gradient is small, and CCA
        when it's large, so that the contingency table is sparse.

        Parameters
        ----------
        Y : array_like
            :math:`n \times p` response matrix. Its columns need be
            dimensionally homogeneous (or you can set `scale_Y=True`).
        X : array_like
            :math:`n \times m, n \geq m` matrix of explanatory
            variables. Its columns need not be standardize, but doing so
            turns regression coefficients into standard regression
            coefficients.
        scale_Y : bool, optional
            Controls whether the response matrix columns are scaled to
            have unit standard deviation. Defaults to `False`.
        """
        self.Y = np.asarray(Y, dtype=np.float64)
        self.X = np.asarray(X, dtype=np.float64)
        self._rda(scale_Y)

    def _rda(self, scale_Y):
        n, p = self.Y.shape
        n_, m = self.X.shape
        if n != n_:
            raise ValueError(
                "Both data matrices must have the same number of rows.")
        if n < m:
            # Mmm actually vegan is able to do this case, too
            raise ValueError(
                "Explanatory variables cannot have less rows than columns.")

        # Centre response variables (they must be dimensionally
        # homogeneous)
        Y = scale(self.Y, with_std=scale_Y)
        # Centre explanatory variables
        X = scale(self.X, with_std=False)

        # Distribution of variables should be examined and transformed
        # if necessary (see paragraph 4 in p. 580 L&L 1998)

        # Compute Y_hat (fitted values by multivariate linear
        # regression, that is, linear least squares). Formula 11.6 in
        # L&L 1998 involves solving the normal equations, but that fails
        # when cond(X) ~ eps**(-0.5). A more expensive but much more
        # stable solution (fails when cond(X) ~ eps**-1) is computed
        # using the QR decomposition of X = QR:
        # (11.6) Y_hat = X [X' X]^{-1} X' Y
        #              = QR [R'Q' QR]^{-1} R'Q' Y
        #              = QR [R' R]^{-1} R'Q' Y
        #              = QR R^{-1} R'^{-1} R' Q' Y
        #              = Q Q' Y
        # and B (matrix of regression coefficients)
        # (11.4) B = [X' X]^{-1} X' Y
        #          = R^{-1} R'^{-1} R' Q' Y
        #          = R^{-1} Q'
        #Q, R = np.linalg.qr(X)
        #Y_hat = Q.dot(Q.T).dot(Y)
        #B = scipy.linalg.solve_triangular(R, Q.T.dot(Y))
        # This works provided X has full rank. When not, you can still
        # fix it using R's pseudoinverse or partitioning R. To avoid any
        # issues, like the numerical instability when trying to
        # reproduce an example in L&L where X was rank-deficient, we'll
        # just use `np.linalg.lstsq`, which uses the SVD decomposition
        # under the hood and so it's also more expensive.
        B, _, rank_X, _ = np.linalg.lstsq(X, Y)
        Y_hat = X.dot(B)
        # Now let's perform PCA on the fitted values from the multiple
        # regression
        u, s, vt = np.linalg.svd(Y_hat, full_matrices=False)
        # vt are the right eigenvectors, which is what we need to
        # perform PCA. That is, we're changing points in Y_hat from the
        # canonical basis to the orthonormal basis given by the right
        # eigenvectors of Y_hat (or equivalently, the eigenvectors of
        # the covariance matrix Y_hat.T.dot(Y_hat))
        # See 3) in p. 583 in L&L 1998
        rank = svd_rank(Y_hat.shape, s)
        assert rank <= min(p, m, n - 1), (
            "There are {} non-zero eigenvalues, but the maximum is {}.".format(
                rank, min(p, m, n - 1)))

        U = vt[:rank].T  # U as in Fig. 11.2

        # Ordination in the space of response variables. Its columns are
        # site scores. (Eq. 11.12)
        F = Y.dot(U)
        # Ordination in the space of explanatory variables. Its columns
        # are fitted site scores. (Eq. 11.13)
        Z = Y_hat.dot(U)

        # Canonical coefficients (formula 11.14)
        #C = B.dot(U)  # Not used

        Y_res = Y - Y_hat
        # PCA on the residuals
        u_res, s_res, vt_res = np.linalg.svd(Y_res, full_matrices=False)
        # See 9) in p. 587 in L&L 1998
        rank_res = svd_rank(Y_res.shape, s_res)
        assert rank_res <= min(p, n - 1), (
            "There are {} non-zero eigenvalues, but the maximum is {}.".format(
                rank_res, min(p, m, n - 1)))

        U_res = vt_res[:rank_res].T
        F_res = Y_res.dot(U_res)  # Ordination in the space of residuals

        # Storing values needed to compute scores
        for val_name, val in (('U', U), ('U_res', U_res), ('F', F),
                              ('F_res', F_res), ('Z', Z), ('u', u[:, :rank])):
            setattr(self, val_name, val)

        self.eigenvalues = np.r_[s[:rank], s_res[:rank_res]]

    def scores(self, scaling):
        """Compute site, species and biplot scores for different scalings.

        Parameters
        ----------
        scaling : int

            TODO: improve explanation based on p. 403 L&L.

            Scaling type 1 produces a distance biplot. It focuses on the
            ordination of sites. Especially interesting when most
            explanatory variables are binary. Not yet implemented.

            Scaling type 2 produces a correlation biplot. It focuses
            on the relationships among explained variables (`Y`).
        """
        if scaling not in {1, 2}:
            raise NotImplementedError("Only scalings 1, 2 available for RDA.")
        # According to the vegan-FAQ.pdf, the scaling factor for scores
        # is (notice that L&L 1998 says in p. 586 that such scaling
        # doesn't affect the interpretation of a biplot)
        eigvals = self.eigenvalues
        const = np.sum(eigvals**2)**0.25
        if scaling == 1:
            scaling_factor = const
        elif scaling == 2:
            scaling_factor = eigvals / const
        species_scores = np.hstack((self.U, self.U_res)) * scaling_factor
        site_scores = np.hstack((self.F, self.F_res)) / scaling_factor
        # TODO not yet used/displayed
        site_constraints = np.hstack((self.Z, self.F_res)) / scaling_factor
        # vegan seems to compute them as corr(self.X[:, :rank_X],
        # self.u) but I don't think that's a good idea. In fact, if
        # you take the example shown in Figure 11.3 in L&L 1998 you
        # can see that there's an arrow for each of the 4
        # environmental variables (depth, coral, sand, other) even if
        # other = not(coral or sand)
        biplot_scores = corr(self.X, self.u)
        # The "Correlations of environmental variables with site
        # scores" from table 11.4 are quite similar to vegan's biplot
        # scores, but they're computed like this:
        #corr(self.X, self.F))
        return OrdinationResults(eigvals=eigvals,
                                 species=species_scores,
                                 site=site_scores,
                                 biplot=biplot_scores,
                                 site_constraints=site_constraints)


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
          some sense, assumes linear responses of ``species'' to
          ``environmental variables''. This is not implemented here.

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
                 correspondence analysis and related multivariate
                 methods in aquatic ecology", Aquatic Sciences 57.3
                 (1995), pp. 255-289.
        """

        self.Y = np.asarray(Y, dtype=np.float64)
        self.X = np.asarray(X, dtype=np.float64)
        self._cca()

    def _cca(self):
        X, Y = self.X, self.Y
        if X.shape[0] != Y.shape[0]:
            raise ValueError("Contingency and environmental tables must have"
                             " the same number of rows (sites)")
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
        # Legendre 1998, p. 612)
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

        self.column_marginals = column_marginals
        self.row_marginals = row_marginals
        self.U = U
        # Storing values needed to compute scores (better way to do this?)
        for val_name, val in (('column_marginals', column_marginals), ('row_marginals', row_marginals), ('U', U), ('U_res', U_res), ('U_hat', U_hat), ('U_hat_res', U_hat_res), ('u', u), ('Y_hat', Y_hat), ('s', s), ('s_res', s_res), ('X_weighted', X_weighted[:, :rank_lstsq])):
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
        Z_scaling1 = (self.row_marginals**-0.5)[:, None] * self.Y_hat.dot(self.U)
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


if __name__ == '__main__':
    import os
    path = os.path.dirname(os.path.abspath(__file__))

    def get_path(fn):
        return os.path.join(path, 'test', 'data', fn)

    X = np.loadtxt(get_path('L&L_CA_data'))
    ordint = CA(X)
    ordint.biplot(1)
    ordint.biplot(2)

    Y = np.loadtxt(get_path('example2_Y'))
    X = np.loadtxt(get_path('example2_X')).reshape(-1, 4, order='F')
    ordint = RDA(Y, X)
    ordint.biplot()

    Y = np.loadtxt(get_path('example3_Y'))
    X = np.loadtxt(get_path('example3_X')).reshape(-1, 4, order='F')
    ordint = CCA(Y, X)
    ordint.biplot()
