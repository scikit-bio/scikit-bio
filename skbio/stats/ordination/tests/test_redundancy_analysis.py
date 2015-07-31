# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


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
        self.env_ids = map(str, range(4))
        self.pc_ids = ['RDA1', 'RDA2', 'RDA3', 'RDA4', 'RDA5', 'RDA6', 'RDA7']

        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example2_Y')),
            index=self.sample_ids, columns=self.feature_ids)

        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example2_X')),
            index=self.sample_ids, columns=self.env_ids)

    def test_scaling1(self):

        scores = rda(self.Y, self.X, scaling=1)

        biplot_scores = pd.DataFrame(np.loadtxt(
            get_data_path('example2_biplot_scaling1')))

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

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_biplot_scaling1')))

        # These are wrong. See issue #1002
        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)
        # These are wrong. See issue #1002
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
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)

    def test_scaling2(self):

        scores = rda(self.Y, self.X, scaling=2)

        biplot_scores = pd.DataFrame(np.loadtxt(
            get_data_path('example2_biplot_scaling2')))

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

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_biplot_scaling2')))

        # These are wrong. See issue #1002
        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)
        # These are wrong. See issue #1002
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
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)


if __name__ == '__main__':
    main()
