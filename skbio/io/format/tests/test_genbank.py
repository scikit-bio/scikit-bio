# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import Protein, DNA, RNA, Sequence
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path
from skbio.io import GenBankFormatError
from skbio.io.format.genbank import (
    _genbank_sniffer,
    _genbank_to_generator, _genbank_to_sequence,
    _genbank_to_dna, _genbank_to_rna, _genbank_to_protein,
    _parse_locus, _parse_reference,
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
                       'unit': 'aa'}},
            None,
            Protein)

        self.single_rna_fp = get_data_path('genbank_single_record')
        imd = IntervalMetadata(63)
        imd.add([(0, 63)],
                [(False, False)],
                {'db_xref': '"taxon:562"',
                 'mol_type': '"mRNA"',
                 'organism': '"Escherichia coli"',
                 'type': 'source',
                 'strand': '+',
                 '__location': '1..63'})
        imd.add([(0, 63)],
                [(False, True)],
                {'phase': 0,
                 'db_xref': ['"taxon:562"', '"taxon:561"'],
                 '__location': '1..>63',
                 'strand': '+',
                 'note': '"alkaline phosphatase signal peptide"',
                 'protein_id': '"AAA23431.1"',
                 'transl_table': '11',
                 'translation': '"MKQSTIALAVLPLLFTPVTKA"',
                 'type': 'CDS'})
        self.single_rna = (
            'gugaaacaaagcacuauugcacuggcugucuuaccguuacuguuuaccccugugacaaaagcc',
            {'ACCESSION': 'M14399',
             'COMMENT': 'Original source text: E.coli, cDNA to mRNA.',
             'DEFINITION': "alkaline phosphatase signal mRNA, 5' end.",
             'KEYWORDS': 'alkaline phosphatase; signal peptide.',
             'LOCUS': {'date': '26-APR-1993',
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
             'VERSION': 'M14399.1'},
            imd,
            RNA)

        # test:
        # 1. multiple records in one file
        # 2. lowercase sequence
        # 3. DNA, RNA, Protein type
        # 4. variation of formats
        self.multi_fp = get_data_path('genbank_multi_records')
        imd_pro = IntervalMetadata(9)
        imd_pro.add([(0, 9)], [(False, False)],
                    {'organism': '"Bacteria"',
                     'type': 'source',
                     'strand': '+',
                     '__location': '1..9'},)
        imd_pro.add([(0, 9)], [(False, True)],
                    {'__location': '1..>9',
                     'product': '"L-carnitine amidase"',
                     'strand': '+',
                     'type': 'Protein'})
        imd_dna = IntervalMetadata(9)
        imd_dna.add([(0, 9)], [(False, False)],
                    {'country': '"Brazil: Parana, Paranavai"',
                     'type': 'source',
                     'strand': '+',
                     '__location': '1..9',
                     'environmental_sample': ''})
        imd_dna.add([(1, 8)], [(True, True)],
                    {'__location': 'complement(<2..>8)',
                     'product': '"16S ribosomal RNA"',
                     'strand': '-',
                     'type': 'rRNA'})

        self.multi = (
            ('gsreildfk',
             {'ACCESSION': 'AAB29917',
              'COMMENT': 'Method: direct peptide sequencing.',
              'DBSOURCE': 'accession AAB29917.1',
              'DEFINITION': 'L-carnitine amidase {N-terminal}',
              'KEYWORDS': '.',
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
             imd_pro,
             Protein),

            ('catgcaggc',
             {'ACCESSION': 'HQ018078',
              'DEFINITION': 'Uncultured Xylanimonas sp.16S, partial',
              'KEYWORDS': 'ENV.',
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
             imd_dna,
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
                                        r'Could not parse the LOCUS line:.*'):
                _parse_locus(line)

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
            seq, md, imd, constructor = self.multi[i]
            exp = constructor(seq, metadata=md, lowercase=True,
                              interval_metadata=imd)
            self.assertEqual(exp, obs)

    def test_genbank_to_sequence(self):
        for i, exp in enumerate(self.multi):
            obs = _genbank_to_sequence(self.multi_fp, seq_num=i+1)
            exp = Sequence(exp[0], metadata=exp[1], lowercase=True,
                           interval_metadata=exp[2])
            self.assertEqual(exp, obs)

    def test_genbank_to_rna(self):
        seq, md, imd, constructor = self.single_rna
        obs = _genbank_to_rna(self.single_rna_fp)
        exp = constructor(seq, metadata=md,
                          lowercase=True, interval_metadata=imd)

        self.assertEqual(exp, obs)

    def test_genbank_to_dna(self):
        i = 1
        exp = self.multi[i]
        obs = _genbank_to_dna(self.multi_fp, seq_num=i+1)
        exp = DNA(exp[0], metadata=exp[1], lowercase=True,
                  interval_metadata=exp[2])

        self.assertEqual(exp, obs)

    def test_genbank_to_protein(self):
        i = 0
        exp = self.multi[i]
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
        seq, md, imd, constructor = self.single
        obj = constructor(seq, md, interval_metadata=imd)
        with io.StringIO() as fh:
            _generator_to_genbank([obj], fh)
            obs = fh.getvalue()

        with open(self.single_lower_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_sequence_to_genbank(self):
        with io.StringIO() as fh:
            for i, (seq, md, imd, constructor) in enumerate(self.multi):
                obj = Sequence(seq, md, interval_metadata=imd, lowercase=True)
                _sequence_to_genbank(obj, fh)
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_dna_protein_to_genbank(self):
        writers = [_protein_to_genbank,
                   _dna_to_genbank]
        with io.StringIO() as fh:
            for i, (seq, md, imd, constructor) in enumerate(self.multi):
                obj = constructor(
                    seq, md, interval_metadata=imd, lowercase=True)
                writers[i](obj, fh)
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_rna_to_genbank(self):
        with io.StringIO() as fh:
            seq, md, imd, constructor = self.single_rna
            obj = constructor(seq, md, interval_metadata=imd, lowercase=True)
            _rna_to_genbank(obj, fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


class RoundtripTests(GenBankIOTests):
    def test_roundtrip_generator(self):
        with io.StringIO() as fh:
            _generator_to_genbank(_genbank_to_generator(self.multi_fp), fh)
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_rna(self):
        with io.StringIO() as fh:
            _rna_to_genbank(_genbank_to_rna(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_dna(self):
        with io.StringIO() as fh:
            _dna_to_genbank(_genbank_to_dna(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_protein(self):
        with io.StringIO() as fh:
            _protein_to_genbank(_genbank_to_protein(self.single_lower_fp), fh)
            obs = fh.getvalue()

        with open(self.single_lower_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_sequence(self):
        with io.StringIO() as fh:
            _sequence_to_genbank(_genbank_to_sequence(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
