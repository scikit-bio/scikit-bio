# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np
import pandas as pd
import numpy.testing as npt

from skbio import OrdinationResults
from skbio.io import OrdinationFormatError
from skbio.io.format.ordination import (
    _ordination_to_ordination_results, _ordination_results_to_ordination,
    _ordination_sniffer)
from skbio.util import get_data_path, assert_ordination_results_equal


class OrdinationTestData(TestCase):
    def setUp(self):
        self.valid_fps = map(
            get_data_path,
            ['ordination_L&L_CA_data_scores', 'ordination_example3_scores',
             'ordination_PCoA_sample_data_3_scores',
             'ordination_example2_scores'])

        # Store filepath, regex for matching the error message that should be
        # raised when reading the file, and whether the file should be matched
        # by the sniffer (True) or not (False).
        self.invalid_fps = map(lambda e: (get_data_path(e[0]), e[1], e[2]), [
            ('empty', r'end of file.*Eigvals header', False),
            ('whitespace_only', r'Eigvals header not found', False),
            ('ordination_error1', r'Eigvals header not found', False),
            ('ordination_error2',
             r'Proportion explained header not found', False),
            ('ordination_error3', r'Species header not found', True),
            ('ordination_error4', r'Site header not found', True),
            ('ordination_error5', r'Biplot header not found', True),
            ('ordination_error6', r'Site constraints header not found', True),
            ('ordination_error7', r'empty line', False),
            ('ordination_error8', r'9.*Proportion explained.*8', True),
            ('ordination_error9', r'2 values.*1 in row 1', True),
            ('ordination_error10', r'2 values.*1 in row 1', True),
            ('ordination_error11', r'Site constraints ids and site ids', True),
            ('ordination_error12', r'9.*Eigvals.*8', True),
            ('ordination_error13', r'9.*Proportion explained.*8', True),
            ('ordination_error14', r'Site is 0: 9 x 0', True),
            ('ordination_error15', r'9 values.*8 in row 1', True),
            ('ordination_error16', r'Biplot is 0: 3 x 0', True),
            ('ordination_error17', r'3 values.*2 in row 1', True),
            ('ordination_error18',
             r'proportion explained.*eigvals: 8 != 9', True),
            ('ordination_error19',
             r'coordinates.*species.*eigvals: 1 != 2', True),
            ('ordination_error20',
             r'coordinates.*site.*eigvals: 1 != 2', True),
            ('ordination_error21', r'one eigval', False),
            ('ordination_error22', r'end of file.*blank line', False),
            ('ordination_error23',
             r'end of file.*Proportion explained section',
             True),
            ('ordination_error24',
             r'end of file.*row 2.*Species section', True)
        ])


class OrdinationResultsReaderWriterTests(OrdinationTestData):
    def setUp(self):
        super(OrdinationResultsReaderWriterTests, self).setUp()

        # define in-memory results, one for each of the valid files in
        # self.valid_fps

        # CA results
        axes_ids = ['CA1', 'CA2']
        species_ids = ['Species1', 'Species2', 'Species3']
        site_ids = ['Site1', 'Site2', 'Site3']
        eigvals = pd.Series([0.0961330159181, 0.0409418140138], axes_ids)
        species = pd.DataFrame([[0.408869425742, 0.0695518116298],
                                [-0.1153860437, -0.299767683538],
                                [-0.309967102571, 0.187391917117]],
                               index=species_ids, columns=axes_ids)
        site = pd.DataFrame([[-0.848956053187, 0.882764759014],
                             [-0.220458650578, -1.34482000302],
                             [1.66697179591, 0.470324389808]],
                            index=site_ids, columns=axes_ids)
        biplot = None
        site_constraints = None
        prop_explained = None
        ca_scores = OrdinationResults(
            'CA', 'Correspondence Analysis', eigvals=eigvals, features=species,
            samples=site, biplot_scores=biplot,
            sample_constraints=site_constraints,
            proportion_explained=prop_explained)

        # CCA results
        axes_ids = ['CCA%d' % i for i in range(1, 10)]
        species_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                       'Species4', 'Species5', 'Species6', 'Species7',
                       'Species8']
        site_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4', 'Site5',
                    'Site6', 'Site7', 'Site8', 'Site9']

        eigvals = pd.Series([0.366135830393, 0.186887643052, 0.0788466514249,
                             0.082287840501, 0.0351348475787, 0.0233265839374,
                             0.0099048981912, 0.00122461669234,
                             0.000417454724117], axes_ids)
        species = pd.DataFrame(np.loadtxt(
            get_data_path('ordination_exp_Ordination_CCA_species')),
            index=species_ids, columns=axes_ids)
        site = pd.DataFrame(
            np.loadtxt(get_data_path('ordination_exp_Ordination_CCA_site')),
            index=site_ids, columns=axes_ids)
        biplot = pd.DataFrame(
            [[-0.169746767979, 0.63069090084, 0.760769036049],
             [-0.994016563505, 0.0609533148724, -0.0449369418179],
             [0.184352565909, -0.974867543612, 0.0309865007541]],
            columns=axes_ids[:3])
        site_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('ordination_exp_Ordination_CCA_site_constraints')),
            index=site_ids, columns=axes_ids)
        prop_explained = None
        cca_scores = OrdinationResults('CCA',
                                       'Canonical Correspondence Analysis',
                                       eigvals=eigvals, features=species,
                                       samples=site, biplot_scores=biplot,
                                       sample_constraints=site_constraints,
                                       proportion_explained=prop_explained)

        # PCoA results
        axes_ids = ['PC%d' % i for i in range(1, 10)]
        species_ids = None
        site_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']
        eigvals = pd.Series([0.512367260461, 0.300719094427, 0.267912066004,
                             0.208988681078, 0.19169895326, 0.16054234528,
                             0.15017695712, 0.122457748167, 0.0], axes_ids)
        species = None
        site = pd.DataFrame(
            np.loadtxt(get_data_path('ordination_exp_Ordination_PCoA_site')),
            index=site_ids, columns=axes_ids)
        biplot = None
        site_constraints = None
        prop_explained = pd.Series([0.267573832777, 0.15704469605,
                                    0.139911863774, 0.109140272454,
                                    0.100111048503, 0.0838401161912,
                                    0.0784269939011, 0.0639511763509, 0.0],
                                   axes_ids)
        pcoa_scores = OrdinationResults('PCoA',
                                        'Principal Coordinate Analysis',
                                        eigvals=eigvals, features=species,
                                        samples=site, biplot_scores=biplot,
                                        sample_constraints=site_constraints,
                                        proportion_explained=prop_explained)

        # RDA results
        axes_ids = ['RDA%d' % i for i in range(1, 8)]
        species_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                       'Species4', 'Species5']
        site_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4', 'Site5',
                    'Site6', 'Site7', 'Site8', 'Site9']
        eigvals = pd.Series([25.8979540892, 14.9825779819, 8.93784077262,
                             6.13995623072, 1.68070536498, 0.57735026919,
                             0.275983624351], axes_ids)
        species = pd.DataFrame(np.loadtxt(
            get_data_path('ordination_exp_Ordination_RDA_species')),
            index=species_ids, columns=axes_ids)
        site = pd.DataFrame(
            np.loadtxt(get_data_path('ordination_exp_Ordination_RDA_site')),
            index=site_ids, columns=axes_ids)
        biplot = pd.DataFrame(
            [[0.422650019179, -0.559142585857, -0.713250678211],
             [0.988495963777, 0.150787422017, -0.0117848614073],
             [-0.556516618887, 0.817599992718, 0.147714267459],
             [-0.404079676685, -0.9058434809, -0.127150316558]],
            columns=axes_ids[:3])
        site_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('ordination_exp_Ordination_RDA_site_constraints')),
            index=site_ids, columns=axes_ids)
        prop_explained = None
        rda_scores = OrdinationResults(
            'RDA', 'Redundancy Analysis', eigvals=eigvals, features=species,
            samples=site, biplot_scores=biplot,
            sample_constraints=site_constraints,
            proportion_explained=prop_explained)

        self.ordination_results_objs = [ca_scores, cca_scores, pcoa_scores,
                                        rda_scores]

    def test_read_valid_files(self):
        for fp, obj in zip(self.valid_fps, self.ordination_results_objs):
            obs = _ordination_to_ordination_results(fp)
            assert_ordination_results_equal(
                obs, obj, ignore_method_names=True,
                ignore_axis_labels=True)

    def test_read_invalid_files(self):
        for invalid_fp, error_msg_regexp, _ in self.invalid_fps:
            with self.assertRaisesRegex(OrdinationFormatError,
                                        error_msg_regexp):
                _ordination_to_ordination_results(invalid_fp)

    def test_write(self):
        for fp, obj in zip(self.valid_fps, self.ordination_results_objs):
            fh = io.StringIO()
            _ordination_results_to_ordination(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with io.open(fp) as fh:
                exp = fh.read()

            npt.assert_equal(obs, exp)

    def test_roundtrip_read_write(self):
        for fp in self.valid_fps:
            # Read.
            obj1 = _ordination_to_ordination_results(fp)

            # Write.
            fh = io.StringIO()
            _ordination_results_to_ordination(obj1, fh)
            fh.seek(0)

            # Read.
            obj2 = _ordination_to_ordination_results(fh)
            fh.close()

            assert_ordination_results_equal(obj1, obj2)


class SnifferTests(OrdinationTestData):
    def setUp(self):
        super(SnifferTests, self).setUp()

    def test_matches_and_nonmatches(self):
        # Sniffer should match all valid files, and will match some invalid
        # ones too because it doesn't exhaustively check the entire file.
        for fp in self.valid_fps:
            self.assertEqual(_ordination_sniffer(fp), (True, {}))

        for fp, _, expected_sniffer_match in self.invalid_fps:
            self.assertEqual(_ordination_sniffer(fp),
                             (expected_sniffer_match, {}))


if __name__ == '__main__':
    main()
