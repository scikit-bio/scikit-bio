# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import io

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata
from skbio import DNA, Sequence
from skbio.io import GFF3FormatError
from skbio.io.format.gff3 import (_yield_record,
                                  _parse_record,
                                  _parse_attr,
                                  _gff3_sniffer,
                                  _gff3_to_interval_metadata,
                                  _interval_metadata_to_gff3,
                                  _gff3_to_generator,
                                  _generator_to_gff3,
                                  _gff3_to_sequence,
                                  _sequence_to_gff3,
                                  _gff3_to_dna,
                                  _dna_to_gff3,
                                  _serialize_interval_metadata)


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

        self.upper_bound = 4641652
        self.imd1 = IntervalMetadata(self.upper_bound)
        self.imd1.add(**intvls[0])
        self.imd1.add(**intvls[1])

        self.imd2 = IntervalMetadata(None)
        self.imd2.add(**intvls[2])
        self.imd2.add(**intvls[3])

        self.imd3 = IntervalMetadata(None)
        self.imd3.add(**intvls[4])

        self.seq_fp = get_data_path('gff3_dna')
        self.seq = Sequence('ATGCATGCATGC',
                            metadata={'id': 'NC_1',
                                      'description': 'species X'})
        self.seq.interval_metadata.add(
            [(0, 9)],
            metadata={'source': 'Prodigal_v2.60',
                      'type': 'gene',
                      'score': '.',
                      'strand': '+',
                      'phase': 0,
                      'ID': 'gene1',
                      'Name': 'FXR'})
        self.dna = DNA(self.seq)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = map(get_data_path, [
            'gff3_multi_record',
            'gff3_single_record',
            'gff3_dna'])
        self.negative_fps = map(get_data_path, [
            'empty',
            'whitespace_only',
            'gff3_bad_missing_directive'])

    def test_positive(self):
        for fp in self.positive_fps:
            self.assertEqual(_gff3_sniffer(fp), (True, {}))

    def test_negative(self):
        for fp in self.negative_fps:
            self.assertEqual(_gff3_sniffer(fp), (False, {}))


class ReaderTests(GFF3IOTests):
    def test_parse_attr(self):
        s = 'Dbxref=GO:000152,GO:001234;Note=fooo'
        obs = _parse_attr(s)
        exp = {'db_xref': 'GO:000152,GO:001234', 'note': 'fooo'}
        self.assertEqual(exp, obs)

    def test_yield_record(self):
        obs = [('data', 'seqid1', ['seqid1\txxx', 'seqid1\tyyy']),
               ('data', 'seqid2', ['seqid2\tzzz'])]
        s = ('seqid1\txxx\n'
             'seqid1\tyyy\n'
             'seqid2\tzzz\n')
        fh = io.StringIO(s)
        for i, j in zip(_yield_record(fh), obs):
            self.assertEqual(i, j)

    def test_parse_record_raise(self):
        chars = 'abc?!'
        for char in chars:
            lines = [
                'ctg123\t.\tgene\t1000\t9000\t.\t+\t%s\tID=gene00001' % char]
            with self.assertRaisesRegex(
                    GFF3FormatError,
                    r"unknown value for phase column: '%s'" % char):
                _parse_record(lines, 10000)

    def test_yield_record_raise(self):
        s = '##gff-version   3\nseq_1 . gene 1 3 . + . ID=gene01\n'
        with io.StringIO(s) as fh:
            with self.assertRaises(GFF3FormatError):
                list(_yield_record(fh))

    def test_gff3_to_interval_metadata(self):
        obs = _gff3_to_interval_metadata(
            self.single_fp, seq_id='Chromosome')

        self.assertEqual(obs, self.imd1)

    def test_gff3_to_interval_metadata_empty(self):
        exp = IntervalMetadata(None)
        obs = _gff3_to_interval_metadata(
            # the seq id does not exist
            self.single_fp, seq_id='foo')
        self.assertEqual(obs, exp)

    def test_gff3_to_interval_metadata_bad(self):
        with self.assertRaisesRegex(GFF3FormatError,
                                    r'do not have 9 columns in this line'):
            _gff3_to_interval_metadata(
                get_data_path('gff3_bad_wrong_columns'),
                seq_id='Chromosome')

    def test_gff3_to_generator(self):
        exps = [('Chromosome', self.imd1),
                ('gi|556503834|ref|NC_000913.3|', self.imd2)]
        obss = _gff3_to_generator(self.multi_fp)
        for obs, exp in zip(obss, exps):
            self.assertEqual(obs, exp)

    def test_gff3_to_generator_empty(self):
        empty_fps = map(get_data_path, ['empty', 'whitespace_only'])
        for empty_fp in empty_fps:
            obs = list(_gff3_to_generator(empty_fp))
            self.assertEqual(obs, [])

    def test_gff3_to_sequence(self):
        obs = _gff3_to_sequence(self.seq_fp)
        self.assertEqual(obs, self.seq)

    def test_gff3_to_dna(self):
        obs = _gff3_to_dna(self.seq_fp)
        self.assertEqual(obs, self.dna)


class WriterTests(GFF3IOTests):
    def test_interval_metadata_to_gff3(self):
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(self.imd1, fh, seq_id='Chromosome')
            # only compare the uncommented lines because the comments are not
            # stored in IntervalMetadata
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        with open(self.single_fp) as f:
            exp = [i.rstrip() for i in f.readlines() if not i.startswith('#')]

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

    def test_interval_metadata_to_gff3_escape(self):
        # test escape of reserved char in GFF3
        exp = 'ctg123\t.\tgene\t1\t9\t.\t.\t.\tID=a%3B%3D%26%2Cb'
        imd = IntervalMetadata(9)
        imd.add([(0, 9)], metadata={
            'type': 'gene', 'ID': 'a;=&,b'})
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(imd, fh, seq_id='ctg123')
            # only compare the uncommented lines because the comments are not
            # stored in IntervalMetadata
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        self.assertEqual([exp], obs)

    def test_interval_metadata_to_gff3_multiple_values(self):
        # test multiple values of db_xref are correctly serialized
        exp = 'ctg123\t.\tgene\t1\t9\t.\t.\t.\tDbxref=GO:000152,GO:001234'
        imd = IntervalMetadata(9)
        imd.add([(0, 9)], metadata={
            'type': 'gene', 'db_xref': ['GO:000152', 'GO:001234']})
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(imd, fh, seq_id='ctg123')
            # only compare the uncommented lines because the comments are not
            # stored in IntervalMetadata
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        self.assertEqual([exp], obs)

    def test_interval_metadata_to_gff3_empty(self):
        imd = IntervalMetadata(None)
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(imd, fh, seq_id='foo')
            obs = fh.getvalue()
        self.assertEqual(obs, '##gff-version 3\n')

    def test_interval_metadata_to_gff3_sub_region(self):
        seq_id = 'NC 7'
        with open(self.multi_fp) as f:
            exp = [i.strip() for i in f if i.startswith(seq_id)]

        with io.StringIO() as fh:
            _serialize_interval_metadata(
                self.imd3, seq_id=seq_id, fh=fh, skip_subregion=False)
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]
        self.assertEqual(exp, obs)

        with io.StringIO() as fh:
            _serialize_interval_metadata(self.imd3, seq_id=seq_id, fh=fh)
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]
        # all the rest lines except the 1st are sub-region lines, so only
        # compare the first line from exp
        self.assertEqual(exp[:1], obs)

    def test_sequence_to_gff3(self):
        with io.StringIO() as fh:
            _sequence_to_gff3(self.seq, fh)
            obs = fh.getvalue()

        with open(self.seq_fp) as fh:
            exp = fh.read()

        self.assertEqual(exp, obs)

    def test_dna_to_gff3(self):
        with io.StringIO() as fh:
            _dna_to_gff3(self.dna, fh)
            obs = fh.getvalue()

        with open(self.seq_fp) as fh:
            exp = fh.read()

        self.assertEqual(exp, obs)

    def test_raise_subregion(self):
        im = IntervalMetadata(None)
        im.add([(0, 3), (7, 9)], metadata={'type': 'gene'})
        with io.StringIO() as fh:
            with self.assertRaises(GFF3FormatError):
                _serialize_interval_metadata(
                    im, seq_id='a', fh=fh, skip_subregion=False)


class RoundtripTests(GFF3IOTests):
    def test_roundtrip_interval_metadata(self):
        ''''''
        with io.StringIO() as fh:
            _interval_metadata_to_gff3(
                _gff3_to_interval_metadata(
                    self.single_fp,
                    seq_id='Chromosome'),
                fh,
                seq_id='Chromosome')
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        with open(self.single_fp) as f:
            exp = [i.rstrip() for i in f.readlines() if not i.startswith('#')]

        self.assertEqual(obs, exp)

    def test_roundtrip_interval_metadata_generator(self):
        with io.StringIO() as fh:
            _generator_to_gff3(
                _gff3_to_generator(self.multi_fp), fh, skip_subregion=False)
            obs = [i for i in fh.getvalue().splitlines()
                   if not i.startswith('#')]

        with open(self.multi_fp) as f:
            exp = [i.rstrip() for i in f.readlines() if not i.startswith('#')]

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
