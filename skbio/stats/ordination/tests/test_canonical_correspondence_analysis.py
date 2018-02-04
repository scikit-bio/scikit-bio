# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
from unittest import TestCase, main

from skbio import OrdinationResults
from skbio.stats.ordination import cca
from skbio.util import get_data_path, assert_ordination_results_equal


class TestCCAErrors(TestCase):
    def setUp(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        self.Y = pd.DataFrame(np.loadtxt(get_data_path('example3_Y')))
        self.X = pd.DataFrame(np.loadtxt(get_data_path('example3_X')))

    def test_shape(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(ValueError):
            cca(Y, X[:-1])

    def test_Y_values(self):
        X, Y = self.X, self.Y
        Y[0, 0] = -1
        with npt.assert_raises(ValueError):
            cca(Y, X)
        Y[0] = 0
        with npt.assert_raises(ValueError):
            cca(Y, X)

    def test_scaling(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(NotImplementedError):
            cca(Y, X, 3)

    def test_all_zero_row(self):
        X, Y = pd.DataFrame(np.zeros((3, 3))), pd.DataFrame(np.zeros((3, 3)))
        with npt.assert_raises(ValueError):
            cca(X, Y)


class TestCCAResults1(TestCase):
    def setUp(self):
        """Data from table 11.3 in Legendre & Legendre 1998
        (p. 590). Loaded results as computed with vegan 2.0-8 and
        compared with table 11.5 if also there."""
        self.feature_ids = ['Feature0', 'Feature1', 'Feature2', 'Feature3',
                            'Feature4', 'Feature5', 'Feature6', 'Feature7',
                            'Feature8']
        self.sample_ids = ['Sample0', 'Sample1', 'Sample2', 'Sample3',
                           'Sample4', 'Sample5', 'Sample6', 'Sample7',
                           'Sample8', 'Sample9']
        self.env_ids = ['Constraint0', 'Constraint1',
                        'Constraint2']
        self.pc_ids = ['CCA1', 'CCA2', 'CCA3', 'CCA4', 'CCA5', 'CCA6', 'CCA7',
                       'CCA8', 'CCA9']
        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example3_Y')),
            columns=self.feature_ids,
            index=self.sample_ids)
        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example3_X'))[:, :-1],
            columns=self.env_ids,
            index=self.sample_ids
            )

    def test_scaling1(self):
        scores = cca(self.Y, self.X, scaling=1)

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_species_scaling1_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_site_scaling1_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_sample_constraints_scaling1')),
            index=self.sample_ids,
            columns=self.pc_ids)
        mat = np.loadtxt(get_data_path(
            'example3_biplot_scaling1'))
        cropped_pcs = self.pc_ids[:mat.shape[1]]
        biplot_scores = pd.DataFrame(mat,
                                     index=self.env_ids,
                                     columns=cropped_pcs)

        proportion_explained = pd.Series([0.466911, 0.238327, 0.100548,
                                          0.104937, 0.044805, 0.029747,
                                          0.012631, 0.001562, 0.000532],
                                         index=self.pc_ids)
        eigvals = pd.Series([0.366136, 0.186888, 0.078847, 0.082288,
                             0.035135, 0.023327, 0.009905, 0.001225,
                             0.000417], index=self.pc_ids)

        exp = OrdinationResults(
            'CCA', 'Canonical Correspondence Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        decimal=6)

    def test_scaling2(self):
        scores = cca(self.Y, self.X, scaling=2)

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_species_scaling2_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_site_scaling2_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_sample_constraints_scaling2')),
            index=self.sample_ids,
            columns=self.pc_ids)

        mat = np.loadtxt(get_data_path(
            'example3_biplot_scaling2'))

        cropped_pc_ids = self.pc_ids[:mat.shape[1]]
        biplot_scores = pd.DataFrame(mat,
                                     index=self.env_ids,
                                     columns=cropped_pc_ids)

        proportion_explained = pd.Series([0.466911, 0.238327, 0.100548,
                                          0.104937, 0.044805, 0.029747,
                                          0.012631, 0.001562, 0.000532],
                                         index=self.pc_ids)
        eigvals = pd.Series([0.366136, 0.186888, 0.078847, 0.082288,
                             0.035135, 0.023327, 0.009905, 0.001225,
                             0.000417], index=self.pc_ids)

        exp = OrdinationResults(
            'CCA', 'Canonical Correspondence Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        decimal=6)


if __name__ == '__main__':
    main()
