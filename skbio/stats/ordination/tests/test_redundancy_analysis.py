# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
from unittest import TestCase, main

from skbio import OrdinationResults
from skbio.stats.ordination import rda
from skbio.util import get_data_path, assert_ordination_results_equal


class TestRDAErrors(TestCase):
    def setUp(self):
        pass

    def test_shape(self):
        for n, p, n_, m in [(3, 4, 2, 1), (3, 4, 3, 10)]:
            Y = pd.DataFrame(np.random.randn(n, p))
            X = pd.DataFrame(np.random.randn(n_, m))
            yield npt.assert_raises, ValueError, rda, Y, X, None, None


class TestRDAResults(TestCase):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there are no written
    # results in L&L.
    def setUp(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        self.sample_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
                           'Site5', 'Site6', 'Site7', 'Site8', 'Site9']
        self.feature_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                            'Species4', 'Species5']
        self.env_ids = list(map(str, range(4)))
        self.pc_ids = ['RDA1', 'RDA2', 'RDA3', 'RDA4', 'RDA5', 'RDA6', 'RDA7']

        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example2_Y')),
            index=self.sample_ids, columns=self.feature_ids)

        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example2_X')),
            index=self.sample_ids, columns=self.env_ids)

    def test_scaling1(self):

        scores = rda(self.Y, self.X, scaling=1)

        sample_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('example2_sample_constraints_scaling1')))

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_species_scaling1_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_site_scaling1_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_sample_constraints_scaling1')),
            index=self.sample_ids,
            columns=self.pc_ids)
        mat = np.loadtxt(get_data_path(
            'example2_biplot_scaling1'))
        cropped_pc_ids = self.pc_ids[:mat.shape[1]]
        biplot_scores = pd.DataFrame(mat,
                                     index=self.env_ids,
                                     columns=cropped_pc_ids)

        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)

        eigvals = pd.Series([25.897954, 14.982578, 8.937841, 6.139956,
                             1.680705, 0.577350, 0.275984],
                            index=self.pc_ids)

        exp = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_directionality=True,
                                        decimal=6)

    def test_scaling2(self):

        scores = rda(self.Y, self.X, scaling=2)
        mat = np.loadtxt(get_data_path('example2_biplot_scaling2'))
        cropped_pc_ids = self.pc_ids[:mat.shape[1]]
        biplot_scores = pd.DataFrame(mat,
                                     index=self.env_ids,
                                     columns=cropped_pc_ids)

        sample_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('example2_sample_constraints_scaling2')))

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_species_scaling2_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_site_scaling2_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_sample_constraints_scaling2')),
            index=self.sample_ids,
            columns=self.pc_ids)

        mat = np.loadtxt(get_data_path(
            'example2_biplot_scaling2'))
        cropped_pc_ids = self.pc_ids[:mat.shape[1]]
        biplot_scores = pd.DataFrame(mat,
                                     index=self.env_ids,
                                     columns=cropped_pc_ids)

        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)

        eigvals = pd.Series([25.897954, 14.982578, 8.937841, 6.139956,
                             1.680705, 0.577350, 0.275984],
                            index=self.pc_ids)

        exp = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_directionality=True,
                                        decimal=6)


class TestRDAResults_biplot_score(TestCase):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there are no written
    # results in L&L.
    def setUp(self):
        """varespec and varechem from Väre etal. 1995 DOI: 10.2307/3236351"""

        self.Y = pd.read_csv(get_data_path('varespec.csv'), index_col=0)
        self.X = pd.read_csv(get_data_path('varechem.csv'), index_col=0)
        self.Y.index.name = None
        self.X.index.name = None

    def test_biplot_score(self):

        rda_ = rda(y=self.Y, x=self.X, scale_Y=False, scaling=1)

        # Load data as computed with vegan 2.4-3:
        # library(vegan)
        # data(varechem)
        # data(varespec)
        # rda_ = rda(X=varespec, Y=varechem, scale=FALSE)
        # write.table(summary(rda_, scaling=1)$biplot,
        #             'vare_rda_biplot_from_vegan.csv', sep=',')
        # write.table(summary(rda_, scaling=1)$sites,
        #                     'vare_rda_sites_from_vegan.csv', sep=',')
        # write.table(summary(rda_, scaling=1)$species,
        #                     'vare_rda_species_from_vegan.csv', sep=',')
        # write.table(summary(rda_, scaling=1)$constraints, #
        #                     'vare_rda_constraints_from_vegan.csv', sep=',')
        # write.table(summary(rda_, scaling=1)$cont$importance[2, ],
        #                     'vare_rda_propexpl_from_vegan.csv', sep=',')
        # write.table(summary(rda_, scaling=1)$cont$importance[1, ],
        #                     'vare_rda_eigvals_from_vegan.csv', sep=',')

        vegan_features = pd.read_csv(
            get_data_path('vare_rda_species_from_vegan.csv'))
        vegan_samples = pd.read_csv(
            get_data_path('vare_rda_sites_from_vegan.csv'))
        vegan_biplot = pd.read_csv(
            get_data_path('vare_rda_biplot_from_vegan.csv'))
        vegan_constraints = pd.read_csv(
            get_data_path('vare_rda_constraints_from_vegan.csv'))
        vegan_propexpl = pd.read_csv(
            get_data_path('vare_rda_propexpl_from_vegan.csv'))
        vegan_propexpl = pd.Series(
            vegan_propexpl.x.values, index=rda_.eigvals.index)
        vegan_eigvals = pd.read_csv(
            get_data_path('vare_rda_eigvals_from_vegan.csv'))
        vegan_eigvals = pd.Series(
            vegan_eigvals.x.values, index=rda_.eigvals.index)

        # scikit-bio returns singular values, whereas vegan returns eigenvalues
        vegan_eigvals = np.sqrt(vegan_eigvals*vegan_eigvals.shape[0])
        vegan_propexpl = vegan_eigvals/vegan_eigvals.sum()

        # transform the output of rda_ to match column selection of vegan
        res_samples = rda_.samples.iloc[:, 0:6]
        res_features = rda_.features.iloc[:, 0:6]

        rda_ = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=res_samples,
            features=res_features,
            sample_constraints=rda_.sample_constraints.iloc[:, 0:6],
            biplot_scores=rda_.biplot_scores.iloc[:, 0:6],
            proportion_explained=rda_.proportion_explained,
            eigvals=rda_.eigvals)

        exp = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=vegan_constraints,
            biplot_scores=vegan_biplot,
            proportion_explained=vegan_propexpl,
            eigvals=vegan_eigvals)

        # This scaling constant is required to make skbio comparable to vegan.
        scaling = (rda_.eigvals.iloc[0] / rda_.eigvals.iloc[:6])
        exp.biplot_scores *= scaling
        assert_ordination_results_equal(
            rda_, exp,
            ignore_directionality=True,
            decimal=6)


if __name__ == '__main__':
    main()
