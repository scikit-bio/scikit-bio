from __future__ import absolute_import, division, print_function
from future.builtins import map, range, zip
import six
import io
import numpy as np
import pandas as pd
import numpy.testing as nptest
from unittest import TestCase, main
from datetime import datetime

from skbio import Protein, DNA, RNA, Sequence
from skbio.util import get_data_path
from skbio.io import GenbankFormatError
from skbio.io.format.genbank import (
    _genbank_sniffer,
    _genbank_to_generator, _genbank_to_biological_sequence,
    _genbank_to_dna, _genbank_to_rna, _genbank_to_protein,
    _parse_locus, _parse_reference,
    _parse_loc_str, _parse_section_default,
    _generator_to_genbank, _biological_sequence_to_genbank,
    _protein_to_genbank, _rna_to_genbank, _dna_to_genbank,
    _serialize_single_genbank,
    _serialize_locus)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'genbank_5_blanks_start_of_file',
            'genbank_single_record',
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


class GenbankIOTests(TestCase):
    # parent class to set up test data for the child class
    def setUp(self):
        # test locus line
        self.locus = (
            (['LOCUS       NC_005816   9609 bp   '
              'DNA   circular   CON   07-FEB-2015'],
             {'division': 'CON', 'mol_type': 'DNA', 'shape': 'circular',
              'locus_name': 'NC_005816', 'date': datetime(2015, 2, 7, 0, 0),
              'unit': 'bp', 'size': 9609}),
            (['LOCUS       SCU49845   5028 bp   '
              'DNA      PLN   21-JUN-1999'],
             {'division': 'PLN', 'mol_type': 'DNA', 'shape': None,
             'locus_name': 'SCU49845', 'date': datetime(1999, 6, 21, 0, 0),
              'unit': 'bp', 'size': 5028}),
            (['LOCUS       NP_001832   360 aa      '
              'linear   PRI   18-DEC-2001'],
             {'division': 'PRI', 'mol_type': None, 'shape': 'linear',
              'locus_name': 'NP_001832', 'date': datetime(2001, 12, 18, 0, 0),
              'unit': 'aa', 'size': 360}))

        self.single_fp = get_data_path('genbank_single_record')
        self.single = (
            'GSREILDFK',
            {'LOCUS': {'date': datetime(1994, 9, 23, 0, 0),
                       'division': 'BCT',
                       'locus_name': 'AAB29917',
                       'mol_type': None,
                       'shape': 'linear',
                       'size': 9,
                       'unit': 'aa'}},
            None,
            Protein)

        self.multi_fp = get_data_path('genbank_multi_records')
        self.multi = (
            ('gsreildfk',
             {'ACCESSION': 'AAB29917',
              'COMMENT': 'Method: direct peptide sequencing.',
              'DBSOURCE': 'accession AAB29917.1',
              'DEFINITION': 'L-carnitine amidase {N-terminal}',
              'FEATURES': [{'db_xref': 'taxon:2',
                            'index_': 0,
                            'left_partial_': False,
                            'location': '1..9',
                            'organism': 'Bacteria',
                            'rc_': False,
                            'right_partial_': False,
                            'type_': 'source'},
                           {'index_': 1,
                            'left_partial_': False,
                            'location': '1..>9',
                            'product': 'L-carnitine amidase',
                            'rc_': False,
                            'right_partial_': True,
                            'type_': 'Protein'}],
              'KEYWORDS': '.',
              'LOCUS': {'date': datetime(1994, 9, 23, 0, 0),
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
                             'TITLE': 'a microbial L-carnitine amidase'}],
              'SOURCE': {'ORGANISM': 'Bacteria',
                         'taxonomy': 'Unclassified.'},
              'VERSION': 'AAB29917.1  GI:545426'},
             pd.DataFrame({0: np.ones(9, dtype=bool),
                           1: np.ones(9, dtype=bool)}),
             Protein),

            ('gugaaacaaagcacuauugcacuggcugucuuaccguuacuguuuaccccugugacaaaagcc',
             # 'gtgaaacaaagcactattgcactggctgtcttaccgttactgtttacccctgtgacaaaagcc',
             {'ACCESSION': 'M14399',
              'COMMENT': 'Original source text: E.coli, cDNA to mRNA.',
              'DEFINITION': u"alkaline phosphatase signal mRNA, 5' end.",
              'FEATURES': [{'db_xref': 'taxon:562',
                            'index_': 0,
                            'left_partial_': False,
                            'location': '1..63',
                            'mol_type': 'mRNA',
                            'organism': 'Escherichia coli',
                            'rc_': False,
                            'right_partial_': False,
                            'type_': 'source'},
                           {'codon_start': '1',
                            'db_xref': ['GI:145230', 'taxon:562', 'taxon:561'],
                            'index_': 1,
                            'left_partial_': False,
                            'location': '1..>63',
                            'note': 'alkaline phosphatase signal peptide',
                            'protein_id': 'AAA23431.1',
                            'rc_': False,
                            'right_partial_': True,
                            'transl_table': '11',
                            'translation': 'MKQSTIALAVLPLLFTPVTKA',
                            'type_': 'CDS'}],
              'KEYWORDS': 'alkaline phosphatase; signal peptide.',
              'LOCUS': {'date': datetime(1993, 4, 26, 0, 0),
                        'division': 'BCT',
                        'locus_name': 'ECOALKP',
                        'mol_type': 'mRNA',
                        'shape': 'linear',
                        'size': 63,
                        'unit': 'bp'},
              'SOURCE': {'ORGANISM': 'Escherichia coli',
                         'taxonomy': 'Bacteria; Proteobacteria; '
                         'Gammaproteobacteria; Enterobacteriales; '
                         'Enterobacteriaceae; Escherichia.'},
              'VERSION': 'M14399.1  GI:145229'},
             pd.DataFrame({0: np.ones(63, dtype=bool),
                           1: np.ones(63, dtype=bool)}),
             RNA),

            ('catgcaggc',
             {'ACCESSION': 'HQ018078',
              'DEFINITION': 'Uncultured Xylanimonas sp.16S, partial',
              'FEATURES': [{'clone': 'ROP4R2E04',
                            'country': 'Brazil: Parana, Paranavai',
                            'db_xref': 'taxon:876087',
                            'environmental_sample': '',
                            'host': 'sugarcane',
                            'index_': 0,
                            'isolation_source': 'rhizosphere soil not '
                            'treated with nitrogen fertilizer',
                            'lat_lon': '23.0728 S 52.4650 W',
                            'left_partial_': False,
                            'location': '1..9',
                            'organism': 'uncultured Xylanimonas sp.',
                            'rc_': False,
                            'right_partial_': False,
                            'type_': 'source'},
                           {'index_': 1,
                            'left_partial_': True,
                            'location': 'complement(<2..>8)',
                            'product': '16S ribosomal RNA',
                            'rc_': True,
                            'right_partial_': True,
                            'type_': 'rRNA'}],
              'KEYWORDS': 'ENV.',
              'LOCUS': {'date': datetime(2010, 8, 29, 0, 0),
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
             pd.DataFrame({0: [True] * 9,
                           1: [False] + [True] * 7 + [False]}),
             DNA))


class ReaderTests(GenbankIOTests):
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
            with self.assertRaisesRegexp(
                    GenbankFormatError, 'Could not parse the LOCUS line:.*'):
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

    def test_parse_loc_str(self):
        length = 12

        examples = [
            '9',  # a single base in the presented sequence
            '3..8',
            '<3..8',
            '1..>8',
            'complement(3..8)',
            'complement(join(3..5,7..9))',
            'join(3..5,7..9)']

        expects = [
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': True, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': True, 'left_partial_': False, 'rc_': False},
             np.array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             np.array([0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0], dtype=bool))]
        for example, expect in zip(examples, expects):
            parsed = _parse_loc_str(example, length)
            self.assertDictEqual(parsed[0], expect[0])
            nptest.assert_equal(parsed[1], expect[1])

    def test_genbank_to_generator_single(self):
        # test single record and uppercase sequence
        fp = get_data_path('genbank_single_record')
        for c in [Sequence, Protein]:
            gb = next(_genbank_to_generator(fp, constructor=c))
            exp = c(self.single[0], self.single[1])
            self.assertEqual(exp, gb)

    def test_genbank_to_generator(self):
        fp = get_data_path('genbank_multi_records')
        for i, gb in enumerate(_genbank_to_generator(fp)):
            seq, md, pmd, constructor = self.multi[i]
            exp = constructor(seq, metadata=md, lowercase=True,
                              positional_metadata=pmd)
            self.assertEqual(exp, gb)

    def test_genbank_to_biological_sequence(self):
        fp = get_data_path('genbank_multi_records')
        for i, exp in enumerate(self.multi):
            gb = _genbank_to_biological_sequence(fp, seq_num=i+1)
            expect = Sequence(
                exp[0], metadata=exp[1], positional_metadata=exp[2])
            self.assertEqual(expect, gb)

    def test_genbank_to_dna(self):
        fp = get_data_path('genbank_multi_records')
        i = 2
        exp = self.multi[i]
        gb = _genbank_to_dna(fp, seq_num=i+1)
        expect = DNA(exp[0], metadata=exp[1], lowercase=True,
                     positional_metadata=exp[2])
        self.assertEqual(expect, gb)

    def test_genbank_to_rna(self):
        fp = get_data_path('genbank_multi_records')
        i = 1
        exp = self.multi[i]
        gb = _genbank_to_rna(fp, seq_num=i+1)
        expect = RNA(exp[0].replace('t', 'u'), metadata=exp[1],
                     lowercase=True, positional_metadata=exp[2])
        self.assertEqual(expect, gb)

    def test_genbank_to_protein(self):
        fp = get_data_path('genbank_multi_records')
        i = 0
        exp = self.multi[i]
        gb = _genbank_to_protein(fp, seq_num=i+1)
        expect = Protein(exp[0], metadata=exp[1],
                         lowercase=True, positional_metadata=exp[2])
        self.assertEqual(expect, gb)


class WriterTests(GenbankIOTests):
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

        with io.open(self.single_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_any_sequence_to_genbank(self):
        writers = {Protein: _protein_to_genbank,
                   RNA: _rna_to_genbank,
                   DNA: _dna_to_genbank}
        fh = io.StringIO()
        for seq, md, pmd, constructor in self.multi:
            obj = constructor(seq, md, pmd, lowercase=True)
            writers[constructor](obj, fh)
            obs = fh.getvalue()
        fh.close()

        with io.open(self.multi_fp) as fh:
            exp = fh.read()
        print(obs)
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
