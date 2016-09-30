# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import io

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata

from skbio.io.format.gff3 import (_gff3_sniffer,
                                  _gff3_to_interval_metadata,
                                  _gff3_to_generator,
                                  _interval_metadata_to_gff3)


class GFF3IOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('gff3_multi_record')
        self.single_fp = get_data_path('gff3_single_record')

        intvls = [{'bounds': [(0, 4641652)],
                   'metadata': {'SOURCE': 'European Nucleotide Archive',
                                'TYPE': 'chromosome',
                                'SCORE': '.',
                                'STRAND': '.',
                                'PHASE': '.',
                                'ATTR': ('ID=chromosome:Chromosome;'
                                         'Alias=U00096.3;Is_circular=true')}},
                  {'bounds': [(147, 148)],
                   'metadata': {'SOURCE': 'regulondb_feature',
                                'TYPE': 'biological_region',
                                'SCORE': '.',
                                'STRAND': '+',
                                'PHASE': '.',
                                'ATTR': ('external_name=Promoter thrLp '
                                         '(RegulonDB:ECK120010236);'
                                         'logic_name=regulondb_promoter')}},
                  {'bounds': [(2, 98)],
                   'metadata': {'SOURCE': 'Prodigal_v2.60',
                                'TYPE': 'CDS',
                                'SCORE': '1.8',
                                'STRAND': '+',
                                'PHASE': '0',
                                'ATTR': ('ID=1_1;partial=10;start_type=Edge;'
                                         'rbs_motif=None;rbs_spacer=None;'
                                         'gc_cont=0.427;conf=60.39;score=1.84;'
                                         'cscore=-0.88;sscore=2.72;rscore=0.00;'
                                         'uscore=0.00;tscore=3.22;')}},
                  {'bounds': [(336, 2799)],
                   'metadata': {'SOURCE': 'Prodigal_v2.60',
                                'TYPE': 'CDS',
                                'SCORE': '333.8',
                                'STRAND': '+',
                                'PHASE': '0',
                                'ATTR': ('ID=1_2;partial=00;start_type=ATG;'
                                         'rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;'
                                         'gc_cont=0.531;conf=99.99;score=333.83;'
                                         'cscore=323.42;sscore=10.41;rscore=9.55;'
                                         'uscore=1.51;tscore=0.00;')}}]

        self.upper_bound1 = 4641652
        self.imd1 = IntervalMetadata(self.upper_bound1)
        self.imd1.add(**intvls[0])
        self.imd1.add(**intvls[1])

        self.upper_bound2 = 2799
        self.imd2 = IntervalMetadata(self.upper_bound2)
        self.imd2.add(**intvls[2])
        self.imd2.add(**intvls[3])


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = map(get_data_path, [
            'gff3_multi_record',
            'gff3_single_record'])
        self.negative_fps = map(get_data_path, [
            'empty',
            'whitespace_only',
            'gff3_bad_missing_directive'])

    def test_positive(self):
        for fp in self.positive_fps:
            self.assertEqual(_gff3_sniffer(fp), (True, {}))
        for fp in self.negative_fps:
            self.assertEqual(_gff3_sniffer(fp), (False, {}))


class ReaderTests(GFF3IOTests):
    def test_gff3_to_interval_metadata(self):
        obs = _gff3_to_interval_metadata(
            self.single_fp, upper_bound=4641652, rec_num=1)

        self.assertEqual(obs, self.imd1)

    def test_gff3_to_generator(self):
        obss = _gff3_to_generator(
            self.multi_fp, upper_bounds=[self.upper_bound1, self.upper_bound2])
        for obs, exp in zip(obss, [self.imd1, self.imd2]):
            self.assertEqual(obs, exp)


class WriterTests(GFF3IOTests):
    def test_interval_metadata_to_gff3(self):
        for _im, _fp in zip([self.simple_exp_im, self.double_exp_im],
                            [self.simple_fp, self.double_fp]):

            fh = io.StringIO()
            _IntervalMetadata_to_gff3(_im, fh)
            obs = fh.getvalue().splitlines()

            with io.open(_fp) as fh:
                exp = fh.read().splitlines()
            self.assertSetEqual(set(exp), set(obs))


class RoundtripTests(GFF3IOTests):
    def test_roundtrip_generator(self):
        ''''''


if __name__ == '__main__':
    main()
