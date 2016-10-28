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
from skbio.io import GFF3FormatError
from skbio.io.format.gff3 import (_yield_record,
                                  _gff3_sniffer,
                                  _gff3_to_interval_metadata,
                                  _gff3_to_generator,
                                  _interval_metadata_to_gff3,
                                  _serialize_imd_gff3)


class GFF3IOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('gff3_multi_record')
        self.single_fp = get_data_path('gff3_single_record')

        intvls = [{'bounds': [(0, 4641652)],
                   'metadata': {'source': 'European Nucleotide Archive',
                                'type': 'chromosome',
                                'score': '.',
                                'strand': '.',
                                'ID': 'chromosome:Chromosome',
                                'Alias': 'U00096.3',
                                'Is_circular': 'true'}},
                  {'bounds': [(147, 148)],
                   'metadata': {'source': 'regulondb_feature',
                                'type': 'biological_region',
                                'score': '.',
                                'strand': '+',
                                'external_name':
                                'Promoter thrLp (RegulonDB:ECK120010236)',
                                'logic_name': 'regulondb_promoter'}},
                  {'bounds': [(336, 2799)],
                   'metadata': {'source': 'Prodigal_v2.60',
                                'type': 'gene',
                                'score': '1.8',
                                'strand': '+',
                                'phase': 0,
                                'ID': '1_1',
                                'gc_cont': '0.427'}},
                  {'bounds': [(336, 2799)],
                   'metadata': {'source': 'Prodigal_v2.60',
                                'type': 'CDS',
                                'score': '333.8',
                                'strand': '+',
                                'phase': 0,
                                'ID': '1_2',
                                'Parent': '1_1',
                                'rbs_motif': 'GGAG/GAGG',
                                'rbs_spacer': '5-10bp'}},
                  {'bounds': [(0, 50), (55, 100)],
                   'metadata': {'source': 'Prodigal_v2.60',
                                'type': 'gene',
                                'score': '1.8',
                                'strand': '+',
                                'phase': 0,
                                'ID': '1_1',
                                'gene': 'FXR receptor'}}]

        self.upper_bound1 = 4641652
        self.imd1 = IntervalMetadata(self.upper_bound1)
        self.imd1.add(**intvls[0])
        self.imd1.add(**intvls[1])

        self.upper_bound2 = 2799
        self.imd2 = IntervalMetadata(self.upper_bound2)
        self.imd2.add(**intvls[2])
        self.imd2.add(**intvls[3])

        self.upper_bound3 = 200
        self.imd3 = IntervalMetadata(self.upper_bound3)
        self.imd3.add(**intvls[4])

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
    def test_yield_record(self):
        obs = [('seqid1', ['seqid1\txxx', 'seqid1\tyyy']),
               ('seqid2', ['seqid2\tzzz'])]
        s = ('seqid1\txxx\n'
             'seqid1\tyyy\n'
             'seqid2\tzzz\n')
        fh = io.StringIO(s)
        for i, j in zip(_yield_record(fh), obs):
            self.assertEqual(i, j)

    def test_gff3_to_interval_metadata(self):
        obs = _gff3_to_interval_metadata(
            self.single_fp,
            interval_metadata=IntervalMetadata(self.upper_bound1),
            rec_num=1)

        self.assertEqual(obs, self.imd1)

    def test_gff3_to_interval_metadata_bad(self):
        with self.assertRaisesRegex(GFF3FormatError,
                                    'do not have 9 columns in this line'):
            _gff3_to_interval_metadata(
                get_data_path('gff3_bad_wrong_columns'),
                interval_metadata=IntervalMetadata(self.upper_bound1))

    def test_gff3_to_generator(self):
        imd = {'Chromosome': IntervalMetadata(self.upper_bound1),
               'gi|556503834|ref|NC_000913.3|':
               IntervalMetadata(self.upper_bound2)}
        obss = _gff3_to_generator(self.multi_fp, interval_metadata_dict=imd)
        for obs, exp in zip(obss, [self.imd1, self.imd2]):
            self.assertEqual(obs, exp)


class WriterTests(GFF3IOTests):
    def test_interval_metadata_to_gff3(self):
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(self.imd1, fh, seq_id='Chromosome')
            # only compare the uncommented lines because the comments are not
            # stored in IntervalMetadata
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        with open(self.single_fp) as f:
            exp = [i[:-1] for i in f.readlines() if not i.startswith('#')]

        self.assertEqual(obs, exp)


    def test_interval_metadata_to_gff3_missing_field(self):
        exp = 'ctg123\t.\tgene\t1\t9\t.\t.\t.\tID=gene00001;Name=EDEN'
        imd = IntervalMetadata(9)
        imd.add([(0, 9)], metadata={
            'type': 'gene', 'ID': 'gene00001', 'Name': 'EDEN'})
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(imd, fh, seq_id='ctg123')
            # only compare the uncommented lines because the comments are not
            # stored in IntervalMetadata
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        self.assertEqual([exp], obs)


    def test_interval_metadata_to_gff3_sub_region(self):
        with open(self.multi_fp) as f:
            lines = [i.strip() for i in f if not i.startswith('#')]

        with io.StringIO() as fh:
            _serialize_imd_gff3(self.imd3, seq_id='NC_7', fh=fh, skip=False)
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]
        exp = lines[-3:]
        self.assertEqual(exp, obs)

        with io.StringIO() as fh:
            _serialize_imd_gff3(self.imd3, seq_id='NC_7', fh=fh)
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]
        exp = [lines[-3]]
        self.assertEqual(exp, obs)


class RoundtripTests(GFF3IOTests):
    def test_roundtrip_interval_metadata(self):
        ''''''
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(
                _gff3_to_interval_metadata(
                    self.single_fp,
                    interval_metadata=IntervalMetadata(self.upper_bound1)),
                fh,
                seq_id='Chromosome')
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        with open(self.single_fp) as f:
            exp = [i[:-1] for i in f.readlines() if not i.startswith('#')]

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
