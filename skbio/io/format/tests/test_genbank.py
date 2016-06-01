# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import numpy.testing as npt
from unittest import TestCase, main

from skbio import Protein, DNA, RNA, Sequence
from skbio.util import get_data_path
from skbio.metadata import Feature
from skbio.io import GenBankFormatError
from skbio.io.format.genbank import (
    _genbank_sniffer,
    _genbank_to_generator, _genbank_to_sequence,
    _genbank_to_dna, _genbank_to_rna, _genbank_to_protein,
    _parse_locus, _parse_reference,
    _parse_interval,
    _parse_section_default,
    _generator_to_genbank, _sequence_to_genbank,
    _protein_to_genbank, _rna_to_genbank, _dna_to_genbank,
    _serialize_locus)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'genbank_5_blanks_start_of_file',
            'genbank_single_record_upper',
            'genbank_single_record_lower',
            'genbank_multi_records']))

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only',
            'genbank_6_blanks_start_of_file',
            'genbank_w_beginning_whitespace',
            'genbank_missing_locus_name']))

    def test_positives(self):
        for fp in self.positive_fps:
            self.assertEqual(_genbank_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negative_fps:
            self.assertEqual(_genbank_sniffer(fp), (False, {}))


class GenBankIOTests(TestCase):
    # parent class to set up test data for the child class
    def setUp(self):
        # test locus line
        self.locus = (
            (['LOCUS       NC_005816   9609 bp   '
              'DNA   circular   CON   07-FEB-2015'],
             {'division': 'CON', 'mol_type': 'DNA', 'shape': 'circular',
              'locus_name': 'NC_005816', 'date': '07-FEB-2015',
              'unit': 'bp', 'size': 9609}),
            (['LOCUS       SCU49845   5028 bp   '
              'DNA      PLN   21-JUN-1999'],
             {'division': 'PLN', 'mol_type': 'DNA', 'shape': None,
             'locus_name': 'SCU49845', 'date': '21-JUN-1999',
              'unit': 'bp', 'size': 5028}),
            (['LOCUS       NP_001832   360 aa      '
              'linear   PRI   18-DEC-2001'],
             {'division': 'PRI', 'mol_type': None, 'shape': 'linear',
              'locus_name': 'NP_001832', 'date': '18-DEC-2001',
              'unit': 'aa', 'size': 360}))

        # test single record and read uppercase sequence
        self.single_upper_fp = get_data_path('genbank_single_record_upper')
        self.single_lower_fp = get_data_path('genbank_single_record_lower')
        self.single = (
            'GSREILDFK',
            {'LOCUS': {'date': '23-SEP-1994',
                       'division': 'BCT',
                       'locus_name': 'AAB29917',
                       'mol_type': None,
                       'shape': 'linear',
                       'size': 9,
                       'unit': 'aa'},
             'id': 'AAB29917'
             },
            None,
            Protein)

        self.single_rna_fp = get_data_path('genbank_single_record')
        self.single_rna = (
            'gugaaacaaagcacuauugcacuggcugucuuaccguuacuguuuaccccugugacaaaagcc',
            {'ACCESSION': 'M14399',
             'COMMENT': 'Original source text: E.coli, cDNA to mRNA.',
             'DEFINITION': u"alkaline phosphatase signal mRNA, 5' end.",
             'KEYWORDS': 'alkaline phosphatase; signal peptide.',
             'LOCUS': {'date': '26-APR-1993',
                       'division': 'BCT',
                       'locus_name': 'ECOALKP',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp'},
             'id': 'ECOALKP',
             'SOURCE': {'ORGANISM': 'Escherichia coli',
                        'taxonomy': 'Bacteria; Proteobacteria; '
                        'Gammaproteobacteria; Enterobacteriales; '
                        'Enterobacteriaceae; Escherichia.'},
             'VERSION': 'M14399.1  GI:145229'},
            {
                Feature(db_xref='"taxon:562"',
                        left_partial_=False,
                        location='1..63',
                        mol_type='"mRNA"',
                        organism='"Escherichia coli"',
                        rc_=False,
                        right_partial_=False,
                        type_='source'): [(0, 63)],
                Feature(codon_start='1',
                        db_xref=('"GI:145230"',
                                 '"taxon:562"',
                                 '"taxon:561"'),
                        left_partial_=False,
                        location='1..>63',
                        note='"alkaline phosphatase signal peptide"',
                        protein_id='"AAA23431.1"',
                        rc_=False,
                        right_partial_=True,
                        transl_table='11',
                        translation='"MKQSTIALAVLPLLFTPVTKA"',
                        type_='CDS'): [(0, 63)]
            },
            RNA)

        # test:
        # 1. multiple records in one file
        # 2. lowercase sequence
        # 3. DNA, RNA, Protein type
        # 4. variation of formats
        self.multi_fp = get_data_path('genbank_multi_records')
        self.multi_invs = (
            ('gsreildfk',
             {'ACCESSION': 'AAB29917',
              'COMMENT': 'Method: direct peptide sequencing.',
              'DBSOURCE': 'accession AAB29917.1',
              'DEFINITION': 'L-carnitine amidase {N-terminal}',
              'KEYWORDS': '.',
              'id': 'AAB29917',
              'LOCUS': {'date': '23-SEP-1994',
                        'division': 'BCT',
                        'locus_name': 'AAB29917',
                        'mol_type': None,
                        'shape': 'linear',
                        'size': 9,
                        'unit': 'aa'},
              'REFERENCE': [{'AUTHORS': 'Joeres,U. and Kula,M.R.',
                             'JOURNAL': 'AMB 40 (5), 606-610 (1994)',
                             'PUBMED': '7764422',
                             'REFERENCE': '1  (residues 1 to 9)',
                             'REMARK': 'from the original journal article.',
                             'TITLE': 'a microbial L-carnitine amidase'},
                            {'AUTHORS': 'Joeres,U. and Kula,M.R.',
                             'JOURNAL': 'AMB 40 (5), 606-610 (1994)',
                             'PUBMED': '7764422',
                             'REFERENCE': '1  (residues 1 to 9)',
                             'TITLE': 'a microbial L-carnitine amidase'}],
              'SOURCE': {'ORGANISM': 'Bacteria',
                         'taxonomy': 'Unclassified.'},
              'VERSION': 'AAB29917.1  GI:545426'},
             {
                              Feature(left_partial_=False,
                                      location='1..9',
                                      organism='"Bacteria"',
                                      rc_=False,
                                      right_partial_=False,
                                      type_='source'): [(0, 9)],
                              Feature(left_partial_=False,
                                      location='1..>9',
                                      product='"L-carnitine amidase"',
                                      rc_=False,
                                      right_partial_=True,
                                      type_='Protein'): [(0, 9)]
             },
             Protein),
            ('catgcaggc',
             {'ACCESSION': 'HQ018078',
              'DEFINITION': 'Uncultured Xylanimonas sp.16S, partial',
              'KEYWORDS': 'ENV.',
              'id': 'HQ018078',
              'LOCUS': {'date': '29-AUG-2010',
                        'division': 'ENV',
                        'locus_name': 'HQ018078',
                        'mol_type': 'DNA',
                        'shape': 'linear',
                        'size': 9,
                        'unit': 'bp'},
              'SOURCE': {'ORGANISM': 'uncultured Xylanimonas sp.',
                         'taxonomy': 'Bacteria; Actinobacteria; '
                         'Micrococcales; Promicromonosporaceae; '
                         'Xylanimonas; environmental samples.'},
              'VERSION': 'HQ018078.1  GI:304421728'},
             {
                              Feature(country='"Brazil: Parana, Paranavai"',
                                      environmental_sample='',
                                      left_partial_=False,
                                      location='1..9',
                                      rc_=False,
                                      right_partial_=False,
                                      type_='source'): [(0, 9)],
                              Feature(left_partial_=True,
                                      location='complement(<2..>8)',
                                      product='"16S ribosomal RNA"',
                                      rc_=True,
                                      right_partial_=True,
                                      type_='rRNA'): [(1, 8)]},
             DNA))


class ReaderTests(GenBankIOTests):
    def test_parse_reference(self):
        lines = '''
REFERENCE   1  (bases 1 to 154478)
  AUTHORS   Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.
  TITLE     Complete structure of the chloroplast genome of
            Arabidopsis thaliana
  JOURNAL   DNA Res. 6 (5), 283-290 (1999)
   PUBMED   10574454'''.split('\n')

        exp = {'AUTHORS': 'Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.',
               'JOURNAL': 'DNA Res. 6 (5), 283-290 (1999)',
               'PUBMED': '10574454',
               'REFERENCE': '1  (bases 1 to 154478)',
               'TITLE': ('Complete structure of the chloroplast genome of'
                         ' Arabidopsis thaliana')}
        self.assertEqual(_parse_reference(lines), exp)

    def test_parse_locus(self):
        for serialized, parsed in self.locus:
            self.assertEqual(_parse_locus(serialized), parsed)

    def test_parse_locus_invalid(self):
        lines = [
            # missing unit
            ['LOCUS       NC_005816               9609 '
             '    DNA     circular CON 07-FEB-2015'],
            # missing division
            ['LOCUS       SCU49845     5028 bp'
             '    DNA                    21-JUN-1999'],
            # wrong date format
            ['LOCUS       NP_001832                360 aa'
             '            linear   PRI 2001-12-18']]
        for line in lines:
            with self.assertRaisesRegex(GenBankFormatError,
                                        'Could not parse the LOCUS line:.*'):
                _parse_locus(line)

    def test_parse_section_default(self):
        lines = [
            ['FOO  blah blah',
             '     blah'],
            ['FOO=blah',
             '    blah'],
            ['FOO']]
        kwargs = [{'join_delimitor': '=', 'return_label': False},
                  {'label_delimitor': '=', 'join_delimitor': '',
                   'return_label': True},
                  {'label_delimitor': '=', 'join_delimitor': '=',
                   'return_label': True}]
        expects = ['blah blah=blah',
                   ('FOO', 'blahblah'),
                   ('FOO', '')]
        for i, j, k in zip(lines, kwargs, expects):
            self.assertEqual(k, _parse_section_default(i, **j))

    def test_parse_interval(self):
        length = 12

        examples = [
            '',
            '9',  # a single base in the presented sequence
            '3..8',
            '<3..8',
            '1..>8',
            'complement(3..8)',
            'complement(join(3..5,7..9))',
            'join(3..5,7..9)',
            'J00194.1:1..9',
            '1.9',
            '1^9']

        expects = [
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             []),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             [8]),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             [(2, 8)]),
            ({'right_partial_': False, 'left_partial_': True, 'rc_': False},
             [(2, 8)]),
            ({'right_partial_': True, 'left_partial_': False, 'rc_': False},
             [(0, 8)]),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             [(2, 8)]),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             [(2, 5), (6, 9)]),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             [(2, 5), (6, 9)]),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             []),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             []),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             [])]
        for example, expect in zip(examples, expects):
            parsed = _parse_interval(example, length)
            self.assertDictEqual(parsed[0], expect[0])
            npt.assert_equal(parsed[1], expect[1])

    def test_parse_interval_invalid(self):
        length = 12
        examples = [
            'abc',
            '3-8']
        for example in examples:
            with self.assertRaisesRegex(GenBankFormatError,
                                        'Could not parse location string: '
                                        '"%s"' % example):
                _parse_interval(example, length)

    def test_genbank_to_generator_single(self):
        # test single record and uppercase sequence
        for c in [Sequence, Protein]:
            obs = next(_genbank_to_generator(
                self.single_upper_fp, constructor=c))
            exp = c(self.single[0], metadata=self.single[1],
                    positional_metadata=self.single[2])
            self.assertEqual(exp, obs)

    def test_genbank_to_generator(self):
        for i, obs in enumerate(_genbank_to_generator(self.multi_fp)):
            seq, md, pmd, constructor = self.multi_invs[i]
            exp = constructor(seq, metadata=md, lowercase=True,
                              interval_metadata=pmd)
            self.assertEqual(exp, obs)

    def test_genbank_to_sequence(self):
        for i, exp in enumerate(self.multi_invs):
            obs = _genbank_to_sequence(self.multi_fp, seq_num=i+1)
            exp = Sequence(exp[0], metadata=exp[1], lowercase=True,
                           interval_metadata=exp[2])
            self.assertEqual(exp.interval_metadata,
                             obs.interval_metadata)
            self.assertEqual(exp, obs)

    def test_genbank_to_rna(self):
        self.maxDiff = None
        seq, md, pmd, constructor = self.single_rna
        obs = _genbank_to_rna(self.single_rna_fp)
        exp = constructor(seq, metadata=md,
                          lowercase=True, interval_metadata=pmd)
        self.assertEqual(exp.interval_metadata,
                         obs.interval_metadata)

        self.assertEqual(exp, obs)

    def test_genbank_to_dna(self):
        i = 1
        exp = self.multi_invs[i]
        obs = _genbank_to_dna(self.multi_fp, seq_num=i+1)
        exp = DNA(exp[0], metadata=exp[1], lowercase=True,
                  interval_metadata=exp[2])
        self.assertEqual(exp, obs)

    def test_genbank_to_protein(self):
        i = 0
        exp = self.multi_invs[i]
        obs = _genbank_to_protein(self.multi_fp, seq_num=i+1)
        exp = Protein(exp[0], metadata=exp[1],
                      lowercase=True, interval_metadata=exp[2])
        self.assertEqual(exp, obs)


class WriterTests(GenBankIOTests):
    def test_serialize_locus(self):
        for serialized, parsed in self.locus:
            self.assertEqual(
                _serialize_locus('LOCUS', parsed), serialized[0] + '\n')

    def test_generator_to_genbank(self):
        seq, md, pmd, constructor = self.single
        obj = constructor(seq, md, pmd)
        fh = io.StringIO()
        _generator_to_genbank([obj], fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_lower_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_sequence_to_genbank(self):
        fh = io.StringIO()
        for i, (seq, md, pmd, constructor) in enumerate(self.multi_invs):
            obj = Sequence(seq, md, interval_metadata=pmd, lowercase=True)

            _sequence_to_genbank(obj, fh)
        obs = fh.getvalue()
        fh.close()
        with io.open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_dna_protein_to_genbank(self):
        writers = [_protein_to_genbank,
                   _dna_to_genbank]
        fh = io.StringIO()
        for i, (seq, md, pmd, constructor) in enumerate(self.multi_invs):
            obj = constructor(seq, md, interval_metadata=pmd, lowercase=True)
            writers[i](obj, fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_rna_to_genbank(self):
        fh = io.StringIO()
        seq, md, pmd, constructor = self.single_rna
        obj = constructor(seq, md, interval_metadata=pmd, lowercase=True)
        _rna_to_genbank(obj, fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


class RoundtripTests(GenBankIOTests):
    def test_roundtrip_generator(self):
        fh = io.StringIO()
        _generator_to_genbank(_genbank_to_generator(self.multi_fp), fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.multi_fp) as fh:
            exp = fh.read()
        self.assertEqual(obs, exp)

    def test_roundtrip_rna(self):
        fh = io.StringIO()
        _rna_to_genbank(_genbank_to_rna(self.single_rna_fp), fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_dna(self):
        fh = io.StringIO()
        _dna_to_genbank(_genbank_to_dna(self.single_rna_fp), fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_protein(self):
        fh = io.StringIO()
        _protein_to_genbank(_genbank_to_protein(self.single_lower_fp), fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_lower_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_sequence(self):
        fh = io.StringIO()
        _sequence_to_genbank(_genbank_to_sequence(self.single_rna_fp), fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
