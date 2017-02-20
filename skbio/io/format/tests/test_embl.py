# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import skbio.io

from skbio import DNA, RNA, Sequence, Protein
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path

# Module specific execption and functions
from skbio.io import EMBLFormatError

from skbio.io.format.embl import (
    _embl_sniffer, _parse_id, _parse_reference, _embl_to_generator,
    _get_embl_section, _embl_to_sequence, _embl_to_dna,
    _embl_to_rna, _embl_to_protein,
    _generator_to_embl, _sequence_to_embl,
    _protein_to_embl, _rna_to_embl, _dna_to_embl,
    _serialize_id)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'embl_single_record',
            'embl_multi_records']))

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only',
            'embl_uniprot_record']))

    def test_positives(self):
        for fp in self.positive_fps:
            self.assertEqual(_embl_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negative_fps:
            self.assertEqual(_embl_sniffer(fp), (False, {}))


# Boilerplate for EMBL IO tests
# TODO: implements all setUp needed
class EMBLIOTests(TestCase):
    def setUp(self):
        # to test ID line
        self.id = (
            # This is a derived record (non-coding, rRNA and spacer records)
            # (feature level record:
            # http://www.ebi.ac.uk/ena/browse/feature-level-products
            # TODO: a Uniprot record?
            (['ID   AB000684.1:1..275:rRNA; SV 1; linear; '
              'genomic DNA; STD; ENV; 275 BP.'],
             {'division': 'ENV', 'mol_type': 'genomic DNA', 'shape': 'linear',
              'accession': 'AB000684.1:1..275:rRNA', 'unit': 'bp',
              'size': 275, 'version': 1, 'class': 'STD'}),
            # A standard record
            (['ID   M14399; SV 1; linear; mRNA; STD; PRO; 63 BP.'],
             {'division': 'PRO', 'mol_type': 'mRNA', 'shape': 'linear',
              'accession': 'M14399', 'unit': 'bp',
              'size': 63, 'version': 1, 'class': 'STD'}))

        # define a single DNA record (with no interval metadata)
        # M14399; SV 1; linear; mRNA; STD; PRO; 63 BP.
        self.single = (
            'gtgaaacaaagcactattgcactggctgtcttaccgttactgtttacccctgtgacaaaagcc',
            {'LOCUS': {'accession': 'M14399',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp',
                       'version': 1}},
            None,
            DNA)

        # define a single protein record (uniprot)
        self.protein = (
            'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKG'
            'LIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAY'
            'NLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFK'
            'ALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWK'
            'FTPL',
            {'LOCUS': {'accession': '001R_FRG3G',
                       'status': 'Reviewed',
                       'size': 256,
                       'unit': 'aa'}},
            None,
            Protein)

        # define a single DNA record uppercase (filepath)
        self.single_upper_fp = get_data_path('embl_single_record_upper')

        # define a single RNA record lower
        self.single_lower_fp = get_data_path('embl_single_record_lower')

        # define a single RNA record file path
        self.single_rna_fp = get_data_path('embl_single_record')

        # define a interval metadata (see skbio.metadata.IntervalMetadata)
        imd = IntervalMetadata(63)

        # then add interval object to interval metadata. Add source
        imd.add([(0, 63)],
                [(False, False)],
                {'db_xref': '"taxon:562"',
                 'mol_type': '"mRNA"',
                 'organism': '"Escherichia coli"',
                 'type': 'source',
                 'strand': '+',
                 '__location': '1..63'})

        imd.add([(0, 63)],
                # the second True is beacause exact location is not known
                [(False, True)],
                {'phase': 0,
                 'db_xref': ['"GOA:P00634"',
                             '"InterPro:IPR001952"',
                             '"InterPro:IPR017849"',
                             '"InterPro:IPR017850"',
                             '"InterPro:IPR018299"',
                             '"PDB:1AJA"',
                             '"PDB:1AJB"',
                             '"PDB:1AJC"',
                             '"PDB:1AJD"',
                             '"PDB:1ALH"',
                             '"PDB:1ALI"',
                             '"PDB:1ALJ"',
                             '"PDB:1ALK"',
                             '"PDB:1ANI"',
                             '"PDB:1ANJ"',
                             '"PDB:1B8J"',
                             '"PDB:1ED8"',
                             '"PDB:1ED9"',
                             '"PDB:1ELX"',
                             '"PDB:1ELY"',
                             '"PDB:1ELZ"',
                             '"PDB:1EW8"',
                             '"PDB:1EW9"',
                             '"PDB:1HJK"',
                             '"PDB:1HQA"',
                             '"PDB:1KH4"',
                             '"PDB:1KH5"',
                             '"PDB:1KH7"',
                             '"PDB:1KH9"',
                             '"PDB:1KHJ"',
                             '"PDB:1KHK"',
                             '"PDB:1KHL"',
                             '"PDB:1KHN"',
                             '"PDB:1URA"',
                             '"PDB:1URB"',
                             '"PDB:1Y6V"',
                             '"PDB:1Y7A"',
                             '"PDB:2ANH"',
                             '"PDB:2G9Y"',
                             '"PDB:2GA3"',
                             '"PDB:2MLX"',
                             '"PDB:2MLY"',
                             '"PDB:2MLZ"',
                             '"PDB:3BDF"',
                             '"PDB:3BDG"',
                             '"PDB:3BDH"',
                             '"PDB:3CMR"',
                             '"PDB:3DPC"',
                             '"PDB:3DYC"',
                             '"PDB:3TG0"',
                             '"PDB:4KM4"',
                             '"PDB:4YR1"',
                             '"PDB:5C66"',
                             '"PDB:5GAD"',
                             '"PDB:5GAF"',
                             '"PDB:5GAG"',
                             '"PDB:5GAH"',
                             '"PDB:5JTL"',
                             '"PDB:5JTM"',
                             '"PDB:5JTN"',
                             '"PDB:5JTO"',
                             '"PDB:5JTP"',
                             '"UniProtKB/Swiss-Prot:P00634"'],
                 '__location': '1..>63',
                 'strand': '+',
                 'note': '"alkaline phosphatase signal peptide"',
                 'protein_id': '"AAA23431.1"',
                 'transl_table': '11',
                 'translation': '"MKQSTIALAVLPLLFTPVTKA"',
                 'type': 'CDS'})

        self.single_rna = (
            'gugaaacaaagcacuauugcacuggcugucuuaccguuacuguuuaccccugugacaaaagcc',
            {'LOCUS': {'accession': 'M14399',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp',
                       'version': 1},
             'ACCESSION': 'M14399;',  # accessions (could be more than one)
             'DATE': ["16-JUL-1988 (Rel. 16, Created)",
                      "02-SEP-1999 (Rel. 60, Last updated, Version 3)"],
             'DBSOURCE': 'MD5; c9b40131b8622946b5aafdf5473b3d43.',
             'DEFINITION': "E.coli alkaline phosphatase signal mRNA, 5' end.",
             'KEYWORDS': "alkaline phosphatase; signal peptide.",
             'SOURCE': {"ORGANISM": "Escherichia coli",
                        'taxonomy': "Bacteria; Proteobacteria; "
                        "Gammaproteobacteria; Enterobacterales; "
                        "Enterobacteriaceae; Escherichia."},
             'REFERENCE': [{'AUTHORS': 'Gray G.L., Baldridge J.S., '
                                       'McKeown K.S., Heyneker H.L., '
                                       'Chang C.N.;',
                            'JOURNAL': 'Gene 39(2-3):247-254(1985).',
                            'CROSS_REFERENCE': 'DOI; 10.1016/0378-1119(85)'
                                               '90319-1. PUBMED; 3912261.',
                            'REFERENCE': '1-63',
                            'TITLE': '"Periplasmic production of correctly '
                                     'processed human growth hormone in '
                                     'Escherichia coli: natural and bacterial '
                                     'signal sequences are '
                                     'interchangeable";'}], },
            imd,
            RNA)

        # define a multi record. File path
        self.multi_fp = get_data_path('embl_multi_records')

        # define interal metadata
        imd1 = IntervalMetadata(275)

        # then add interval object to interval metadata. Add source
        imd1.add([(0, 275)],
                 [(False, False)],
                 {'db_xref': '"taxon:77133"',
                  'mol_type': '"genomic DNA"',
                  'organism': '"uncultured bacterium"',
                  'type': 'source',
                  'strand': '+',
                  'environmental_sample': '',
                  'clone': '"2 of cluster I"',
                  '__location': '1..275'})

        imd1.add([(0, 275)],
                 [(True, True)],
                 {'product': '"16S rRNA"',
                  'type': 'rRNA',
                  'strand': '+',
                  '__location': "AB000684.1:<1..>275"})

        # define interal metadata
        imd2 = IntervalMetadata(275)

        # then add interval object to interval metadata. Add source
        imd2.add([(0, 275)],
                 [(False, False)],
                 {'db_xref': '"taxon:77133"',
                  'mol_type': '"genomic DNA"',
                  'organism': '"uncultured bacterium"',
                  'type': 'source',
                  'strand': '+',
                  'environmental_sample': '',
                  'clone': '"3 of cluster I"',
                  '__location': '1..275'})

        imd2.add([(0, 275)],
                 [(True, True)],
                 {'product': '"16S rRNA"',
                  'type': 'rRNA',
                  'strand': '+',
                  '__location': "AB000685.1:<1..>275"})

        # multi object
        self.multi = (
            ('CTGCCCTTAGTTGCCACTCTTCGGAGGGCACTCTAAGGGGACCGCCGGCGATAAGCCGAGGAA'
             'GGTGGGGATGACGTCAGGTCAGTATGCCCTTTATGCCCGGGGCTACACAGGCGCTACAGTGGC'
             'CAGGACAATGGGAAGCGACCCAGTAATGGGGAGCAAATCCCTAAACCTGGTCATGGTGCAGAT'
             'TGAGGGCTGAAACTCGCCCCTCATGAAGCCGGAATCGGTAGTAATGGCGGATCAGCTAAGCCG'
             'CCGTGAATACGTTCTCGGGCCTT',
             {'PARENT_ACCESSION': 'AB000684.1',
              'DATE': ['10-FEB-1997 (Rel. 50, Created)',
                       '29-JUL-2005 (Rel. 84, Last updated, Version 7)'],
              'DBSOURCE': 'MD5; 996348a90e49caf3f6155ad9478d6d90.',
              'DEFINITION': 'uncultured bacterium partial 16S rRNA',
              'KEYWORDS': 'ENV.',
              'LOCUS': {'accession': 'AB000684.1:1..275:rRNA',
                        'class': 'STD',
                        'division': 'ENV',
                        'mol_type': 'genomic DNA',
                        'shape': 'linear',
                        'size': 275,
                        'unit': 'bp',
                        'version': 1},
              'REFERENCE': [{'AUTHORS': 'Inagaki F., Hayashi S., Doi K., '
                                        'Motomura Y., Izawa E., Ogata S.;',
                             'JOURNAL': 'Submitted (24-JAN-1997) to the INSDC'
                                        '. Fumio Inagaki, Faculty of '
                                        'Agriculture, Kyushu University, '
                                        'Microbial Genetic Division, Institu'
                                        'te of Genetic Resourses; Higashi-ku'
                                        ' Hakozaki 6-10-1, Fukuoka-shi, '
                                        'Fukuoka, 812-81, Japan (E-mail'
                                        ':inagaki@agr.kyushu-u.ac.jp, Tel:'
                                        '+81-92-642-3059, Fax:+81-92-642-'
                                        '3059)',
                             'REFERENCE': '1-275',
                             'TITLE': ';'},
                            {'AUTHORS': 'Inagaki F., Hayashi S., Doi K., '
                                        'Motomura Y., Izawa E., Ogata S.;',
                             'JOURNAL': 'FEMS Microbiol. Ecol. 24:41-48'
                                        '(1997).',
                             'TITLE': '"Microbial participation in the for'
                                      'mation of siliceous deposits from '
                                      'geothermal water and analysis of the '
                                      'extremely thermophilic bacterial '
                                      'community";'}],
              'SOURCE': {'ORGANISM': 'uncultured bacterium',
                         'taxonomy': 'Bacteria; environmental samples.'}},
             imd1,
             DNA),
            ('CTGCCCTTAGTTGCCACCCTTCGGAGGGCACTCTAAGGGGACCGCCGGCGATAAGCCGAGGAA'
             'GGTGGGGATGACGTCAGGTCAGTATGCCCTTTATGCCCGGGGCTACACAGGCGCTACAGTGGC'
             'CAGGACAATGGGAAGCGACCCAGTAATGGGGAGCAAATCCCTAAACCTGGTCATGGTGCAGAT'
             'TGAGGGCTGAAACTCGCCCCTCATGAAGCCGGAAACGGTAGTAATGGCGGATCAGCTAAGCCG'
             'CCGTGAATACGTTCTCGGGCCTT',
             {'PARENT_ACCESSION': 'AB000685.1',
              'DATE': ['10-FEB-1997 (Rel. 50, Created)',
                       '29-JUL-2005 (Rel. 84, Last updated, Version 7)'],
              'DBSOURCE': 'MD5; 9999acf6687667a872da41b64d3c11f8.',
              'DEFINITION': 'uncultured bacterium partial 16S rRNA',
              'KEYWORDS': 'ENV.',
              'LOCUS': {'accession': 'AB000685.1:1..275:rRNA',
                        'class': 'STD',
                        'division': 'ENV',
                        'mol_type': 'genomic DNA',
                        'shape': 'linear',
                        'size': 275,
                        'unit': 'bp',
                        'version': 1},
              'REFERENCE': [{'AUTHORS': 'Inagaki F., Hayashi S., Doi K., '
                                        'Motomura Y., Izawa E., Ogata S.;',
                                        'JOURNAL': 'Submitted (24-JAN-1997) '
                                        'to the INSDC. Fumio Inagaki, Faculty'
                                        ' of Agriculture, Kyushu University, '
                                        'Microbial Genetic Division, '
                                        'Institute of Genetic Resourses; '
                                        'Higashi-ku Hakozaki 6-10-1, '
                                        'Fukuoka-shi, Fukuoka, 812-81, '
                                        'Japan (E-mail:inagaki@agr.kyushu-'
                                        'u.ac.jp, Tel:+81-92-642-3059, Fax:'
                                        '+81-92-642-3059)',
                             'REFERENCE': '1-275',
                             'TITLE': ';'},
                            {'AUTHORS': 'Inagaki F., Hayashi S., Doi K., '
                                        'Motomura Y., Izawa E., Ogata S.;',
                             'JOURNAL': 'FEMS Microbiol. Ecol. 24:41-'
                                        '48(1997).',
                             'TITLE': '"Microbial participation in the '
                                      'formation of siliceous deposits from '
                                      'geothermal water and analysis of the '
                                      'extremely thermophilic bacterial '
                                      'community";'}],
              'SOURCE': {'ORGANISM': 'uncultured bacterium',
                         'taxonomy': 'Bacteria; environmental samples.'}},
             imd2,
             DNA))


class ReaderTests(EMBLIOTests):
    """Implements test for reading EMBL data"""

    # TODO: implement test to deal with all EMBL methods
    def test_parse_id(self):
        """Parse ID record (first line of embl format)"""
        for serialized, parsed in self.id:
            self.assertEqual(_parse_id(serialized), parsed)

    def test_parse_id_invalid(self):
        lines = [
            # uniprot line (should this module handle it?)
            ['ID   001R_FRG3G              Reviewed;         256 AA.'],
            # missing unit
            ['ID   M14399; SV 1; linear; mRNA; STD; PRO; 63'],
            # missing division
            ['ID   M14399; SV 1; linear; mRNA; STD;      63 BP.']]

        for line in lines:
            with self.assertRaisesRegex(EMBLFormatError,
                                        'Could not parse the ID line:.*'):
                _parse_id(line)

    def test_no_protein_support(self):
        """Testing no protein support for embl"""
        # TODO: add protein support

        # a fake protein line.
        handle = io.StringIO('ID   M14399; SV 1; linear; mRNA; STD; '
                             'PRO; 63 AA.\n//\n')

        with self.assertRaisesRegex(EMBLFormatError,
                                    "There's no protein support for EMBL "
                                    "record"):
            # read a protein record
            Protein.read(handle)

        # return to 0
        handle.seek(0)

        with self.assertRaisesRegex(EMBLFormatError,
                                    "There's no protein support for EMBL "
                                    "record"):
            # read a generic record
            skbio.io.read(handle, format='embl')

    def test_parse_reference(self):
        lines = '''
RP   1-63
RX   DOI; 10.1016/0378-1119(85)90319-1.
RX   PUBMED; 3912261.
RA   Gray G.L., Baldridge J.S., McKeown K.S., Heyneker H.L., Chang C.N.;
RT   "Periplasmic production of correctly processed human growth hormone in
RT   Escherichia coli: natural and bacterial signal sequences are
RT   interchangeable";
RL   Gene 39(2-3):247-254(1985).'''.split('\n')

        exp = {'AUTHORS': 'Gray G.L., Baldridge J.S., '
                          'McKeown K.S., Heyneker H.L., Chang C.N.;',
               'JOURNAL': 'Gene 39(2-3):247-254(1985).',
               'CROSS_REFERENCE': 'DOI; 10.1016/0378-1119(85)90319-1. '
                                  'PUBMED; 3912261.',
               'REFERENCE': '1-63',
               'TITLE': '"Periplasmic production of correctly processed '
                        'human growth hormone in Escherichia coli: '
                        'natural and bacterial signal sequences are '
                        'interchangeable";'
               }

        # See all differences
        self.maxDiff = None
        self.assertEqual(_parse_reference(lines), exp)

    def test_embl_to_generator_single(self):
        # test single record and uppercase sequence
        for c in [Sequence, DNA]:
            obs = next(_embl_to_generator(
                self.single_upper_fp, constructor=c))
            exp = c(self.single[0], metadata=self.single[1],
                    positional_metadata=self.single[2], lowercase=True)
            self.assertEqual(exp, obs)

    def test_get_embl_section(self):
        """Verify to have a section for each embl ID"""

        with open(self.single_rna_fp) as fh:
            for line in fh:
                if line.startswith("//"):
                    continue

                # test that this function doesn't raise exceptions
                try:
                    _get_embl_section(line)

                except KeyError as err:
                    raise EMBLFormatError("Key {0} isn't defined in embl."
                                          "KEYS_2_SECTIONS".format(err))

    def test_embl_to_generator(self):
        for i, obs in enumerate(_embl_to_generator(self.multi_fp)):
            seq, md, imd, constructor = self.multi[i]
            exp = constructor(seq, metadata=md, lowercase=True,
                              interval_metadata=imd)
            self.assertEqual(exp, obs)

    def test_embl_to_sequence(self):
        for i, exp in enumerate(self.multi):
            obs = _embl_to_sequence(self.multi_fp, seq_num=i+1)
            exp = Sequence(exp[0], metadata=exp[1], lowercase=True,
                           interval_metadata=exp[2])
            self.assertEqual(exp, obs)

    def test_embl_to_rna(self):
        seq, md, imd, constructor = self.single_rna
        obs = _embl_to_rna(self.single_rna_fp)
        exp = constructor(seq, metadata=md,
                          lowercase=True, interval_metadata=imd)

        self.assertEqual(exp, obs)

    def test_embl_to_dna(self):
        i = 1
        exp = self.multi[i]
        obs = _embl_to_dna(self.multi_fp, seq_num=i+1)
        exp = DNA(exp[0], metadata=exp[1], lowercase=True,
                  interval_metadata=exp[2])

        self.assertEqual(exp, obs)

    def test_embl_to_protein(self):
        # TODO: add protein support
        i = 0
        # there is no support for protein at the moment
#        exp = self.multi[i]
#        obs = _embl_to_protein(self.multi_fp, seq_num=i+1)
#        exp = Protein(exp[0], metadata=exp[1],
#                      lowercase=True, interval_metadata=exp[2])
#        self.assertEqual(exp, obs)

        with self.assertRaisesRegex(EMBLFormatError,
                                    "There's no protein support for EMBL "
                                    "record"):
            # read a generic record
            _embl_to_protein(self.multi_fp, seq_num=i+1)


class WriterTests(EMBLIOTests):
    def test_serialize_id(self):
        for serialized, parsed in self.id:
            self.assertEqual(
                _serialize_id('ID', parsed), serialized[0] + '\n')

    def test_generator_to_embl(self):
        seq, md, imd, constructor = self.single
        obj = constructor(seq, md, interval_metadata=imd, lowercase=True)
        with io.StringIO() as fh:
            _generator_to_embl([obj], fh)
            obs = fh.getvalue()

        with open(self.single_lower_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_sequence_to_embl(self):
        with io.StringIO() as fh:
            for i, (seq, md, imd, constructor) in enumerate(self.multi):
                obj = Sequence(seq, md, interval_metadata=imd)
                _sequence_to_embl(obj, fh)
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_dna_to_embl(self):
        with io.StringIO() as fh:
            for i, (seq, md, imd, constructor) in enumerate(self.multi):
                obj = constructor(
                    seq, md, interval_metadata=imd, lowercase=True)
                _dna_to_embl(obj, fh)

            # read all records written
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_protein_to_embl(self):
        seq, md, imd, constructor = self.protein
        obj = constructor(seq, md, interval_metadata=imd)

        with io.StringIO() as fh:
            self.assertRaisesRegex(EMBLFormatError,
                                   "There's no protein support for EMBL "
                                   "record",
                                   _protein_to_embl, [obj], fh)

    def test_rna_to_embl(self):
        with io.StringIO() as fh:
            seq, md, imd, constructor = self.single_rna
            obj = constructor(seq, md, interval_metadata=imd, lowercase=True)
            _rna_to_embl(obj, fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


class RoundtripTests(EMBLIOTests):
    def test_roundtrip_generator(self):
        with io.StringIO() as fh:
            _generator_to_embl(_embl_to_generator(self.multi_fp), fh)
            obs = fh.getvalue()

        with open(self.multi_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_rna(self):
        with io.StringIO() as fh:
            _rna_to_embl(_embl_to_rna(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_roundtrip_dna(self):
        with io.StringIO() as fh:
            _dna_to_embl(_embl_to_dna(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    # TODO: test_roundtrip_protein

    def test_roundtrip_sequence(self):
        with io.StringIO() as fh:
            _sequence_to_embl(_embl_to_sequence(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
