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
    _serialize_id, _parse_assembly, _embl_parse_section_default,
    _serialize_dbsource)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'embl_single_record',
            'embl_multi_records']))

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only',
            'embl_uniprot_record',
            'embl_w_beginning_whitespace']))

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
              'locus_name': 'AB000684.1:1..275:rRNA', 'unit': 'bp',
              'size': 275, 'version': 1, 'class': 'STD', 'date': None}),
            # A standard record
            (['ID   M14399; SV 1; linear; mRNA; STD; PRO; 63 BP.'],
             {'division': 'PRO', 'mol_type': 'mRNA', 'shape': 'linear',
              'locus_name': 'M14399', 'unit': 'bp',
              'size': 63, 'version': 1, 'class': 'STD', 'date': None}))

        # define a single DNA record (with no interval metadata)
        # M14399; SV 1; linear; mRNA; STD; PRO; 63 BP.
        self.single = (
            'gtgaaacaaagcactattgcactggctgtcttaccgttactgtttacccctgtgacaaaagcc',
            {'LOCUS': {'locus_name': 'M14399',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp',
                       'version': 1,
                       'date': None}},
            None,
            DNA)

        # define a single protein record (uniprot)
        self.protein = (
            'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKG'
            'LIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAY'
            'NLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFK'
            'ALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWK'
            'FTPL',
            {'LOCUS': {'locus_name': '001R_FRG3G',
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

        # define a http://www.ebi.ac.uk/ena/browse/feature-level-products
        self.feature_level_fp = get_data_path("embl_feature_level_record")

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
            {'LOCUS': {'locus_name': 'M14399',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp',
                       'version': 1,
                       'date': '02-SEP-1999'},
             'ACCESSION': 'M14399;',  # accessions (could be more than one)
             'VERSION': 'M14399.1',  # a genbank like version
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
                            'REFERENCE': '1  (bases 1 to 63)',
                            'TITLE': '"Periplasmic production of correctly '
                                     'processed human growth hormone in '
                                     'Escherichia coli: natural and bacterial '
                                     'signal sequences are '
                                     'interchangeable";',
                            'PUBMED': '3912261'}],
             'CROSS_REFERENCE': ['DOI; 10.1016/0378-1119(85)'
                                 '90319-1. PUBMED; 3912261.']},
            imd,
            RNA)

        # define a multi record. File path
        self.multi_fp = get_data_path('embl_multi_records')

        # define interval metadata (as single metadata)
        imd1 = imd

        # define interal metadata for multi 2
        imd2 = IntervalMetadata(743)

        # then add interval object to interval metadata. Add source
        imd2.add([(0, 743)],
                 [(False, False)],
                 {'organism': '"Ruditapes philippinarum"',
                  'type': 'source',
                  '__location': '1..743',
                  'strand': '+',
                  'mol_type': '"mRNA"',
                  'db_xref': '"taxon:129788"'})

        imd2.add([(57, 444)],
                 [(False, False)],
                 {'translation': '"MPGGKAGKDSGKAKAKAVSRSARAGLQFPVGRIHRHLKNRT'
                                 'TSHG RVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRI'
                                 'TPRHLQLAIRGDEELDSLIKAT IAGGGVIPHIHKSLIGKKG'
                                 'GQQAK"',
                  'type': 'CDS',
                  '__location': '58..444',
                  'protein_id': '"APY18893.1"',
                  'strand': '+',
                  'phase': 0,
                  'product': '"histone"'})

        # multi object
        self.multi = (
            ('GTGAAACAAAGCACTATTGCACTGGCTGTCTTACCGTTACTGTTTACCCCTGTGACAAAAGCC',
             {'LOCUS': {'locus_name': 'M14399',
                        'class': 'STD',
                        'division': 'PRO',
                        'mol_type': 'mRNA',
                        'shape': 'linear',
                        'size': 63,
                        'unit': 'bp',
                        'version': 1,
                        'date': '02-SEP-1999'},
              'ACCESSION': 'M14399;',  # accessions (could be more than one)
              'VERSION': 'M14399.1',  # a genbank like version
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
                             'REFERENCE': '1  (bases 1 to 63)',
                             'TITLE': '"Periplasmic production of correctly '
                                      'processed human growth hormone in '
                                      'Escherichia coli: natural and '
                                      'bacterial signal sequences are '
                                      'interchangeable";',
                             'PUBMED': '3912261'}],
              'CROSS_REFERENCE': ['DOI; 10.1016/0378-1119(85)'
                                  '90319-1. PUBMED; 3912261.']},
             imd1,
             DNA),
            ('TGTGCACAGTCTACGCGTCATCTTGAAAGAAAGAACTACACTACTCCAAAAATAATCATGCC'
             'TGGTGGAAAAGCTGGTAAAGATTCCGGAAAGGCCAAGGCTAAGGCAGTGTCAAGGTCCGCAA'
             'GAGCTGGCTTACAGTTTCCAGTCGGACGTATTCACAGGCATTTGAAGAACAGAACCACTAGC'
             'CACGGTCGTGTTGGAGCTACAGCAGCCGTTTACAGTGCAGCAATCCTTGAATACCTGACCGC'
             'CGAAGTGCTTGAGTTGGCTGGAAACGCAAGTAAAGATCTCAAAGTAAAGAGAATCACCCCAC'
             'GTCACTTGCAGTTGGCAATCAGAGGAGATGAAGAGTTGGATTCCCTAATTAAAGCCACAATC'
             'GCTGGTGGTGGTGTTATTCCACATATCCACAAGTCACTTATTGGCAAGAAGGGAGGTCAGCA'
             'AGCCAAATAAATTGGACATACTCATTCATCAGGGAACAATGTGTAGTGAATGTGTTAAAAAG'
             'AACAATCTCATTGTGTAGCTCTTTAGTTTTATATGAATGTGTTAACATGGTCATTCACATCG'
             'TATGACTCATAGAATCATCTGTGTATCATTTCATCCTCTCATTTTATAGCTCCTCATTTTCC'
             'TTAGACTCATTAAAATTTTTATCTCGGAAAAATGTTTTTTCTACAATTTTAGCATTCATTTA'
             'TCTTCATCTTGCTTTTATGTTTAATAAAACGAACTTATAATACCAAAAAAAAAAAAAAAAA',
             {'ACCESSION': 'KX454487;',
              'VERSION': 'KX454487.1',
              'COMMENT': '##Assembly-Data-START##\nSequencing Technology '
                         ':: Sanger dideoxy sequencing\n##Assembly-Data-END##',
              'DATE': ['02-FEB-2017 (Rel. 131, Created)',
                       '02-FEB-2017 (Rel. 131, Last updated, Version 1)'],
              'DBSOURCE': 'MD5; cbc730cf7a8d694b50fb7dd6b993ae0d.',
              'DEFINITION': 'Ruditapes philippinarum histone mRNA, '
                            'complete cds.',
              'KEYWORDS': '.',
              'LOCUS': {'locus_name': 'KX454487',
                        'class': 'STD',
                        'division': 'INV',
                        'mol_type': 'mRNA',
                        'shape': 'linear',
                        'size': 743,
                        'unit': 'bp',
                        'version': 1,
                        'date': '02-FEB-2017'},
              'REFERENCE': [
                {'AUTHORS': 'Yang D., Zhao J., Wang Q.;',
                 'JOURNAL': 'Submitted (27-JUN-2016) to the INSDC. Key '
                            'Laboratory of Coastal Zone Environment Processes '
                            'and Ecological Remediation, Yantai Institute '
                            'of Coastal Zone Research (YIC), Chinese Academy '
                            'of Sciences (CAS), 17 Chunhui Road, Laishan '
                            'District, Yantai, Shandong 264003, China',
                 'REFERENCE': '1  (bases 1 to 743)',
                 'TITLE': ';'}],
              'CROSS_REFERENCE': [None],
              'SOURCE': {
                'ORGANISM': 'Ruditapes philippinarum',
                'taxonomy': 'Eukaryota; Metazoa; Lophotrochozoa; Mollusca; '
                            'Bivalvia; Heteroconchia; Euheterodonta; '
                            'Veneroida; Veneroidea; Veneridae; Ruditapes.'}},
             imd2,
             DNA))

        # define the feature level product obj
        self.feature_level = (
            'AAUUGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUC'
            'GAGCGGCAGCACAGAGGAACUUGUUCCUUGGGUGGCGAGCGGCGGACGGGUGAGUAAUGCCU'
            'AGGAAAUUGCCCUGAUGUGGGGGAUAACCAUUGGAAACGAUGGCUAAUACCGCAUGAUGCCU'
            'ACGGGCCAAAGAGGGGGACCUUCUGGCCUCUCGCGUCAGGAUAUGCCUAGGUGGGAUUAGCU'
            'AGUUGGUGAGGUAAUGGCUCACCAAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGAUCAGC'
            'CACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACA'
            'AUGGGCGCAAGCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUA'
            'CUUUCAGUCGUGAGGAAGGUGGUGUUGUUAAUAGCAGCAUCAUUUGACGUUAGCGACAGAAG'
            'AAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCGAGCGUUAAUCGGA'
            'AUUACUGGGCGUAAAGCGCAUGCAGGUGGUGGAUUAAGUCAGAUGUGAAAGCCCGGGGCUCA'
            'ACCUCGGAACCGCAUUUGAAACUGGUUCACUAGAGUACUGUAGAGGGGGGUAGAAUUUCAGG'
            'UGUAGCGGUGAAAUGCGUAGAGAUCUGAAGGAAUACCGGUGGCGAAGGCGGCCCCCUGGACA'
            'GAUACUGACACUCAGAUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCA'
            'CGCCGUAAACGAUGUCUACUUGGAGGUUGUGGCCUUGAGCCGUGGCUUUCGGAGCUAACGCG'
            'UUAAGUAGACCGCCUGGGGAGUACGGUCGCAAGAUUAAAACUCAAAUGAAUUGACGGGGGCC'
            'CGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAUGCAACGCGAAGAACCUUACCUACUCUUG'
            'ACAUCCAGAGAAGCCAGCGGAGACGCAGGUGUGCCUUCGGGAGCUCUGAGACAGGUGCUGCA'
            'UGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUA'
            'UCCUUGUUUGCCAGCGAGUCAUGUCGGGAACUCCAGGGAGACUGCCGGUGAUAAACCGGAGG'
            'AAGGUGGGGACGACGUCAAGUCAUCAUGGCCCUUACGAGUAGGGCUACACACGUGCUACAAU'
            'GGCGCAUACAGAGGGCAGCAAGCUAGCGAUAGUGAGCGAAUCCCAAAAAGUGCGUCGUAGUC'
            'CGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUAGAUCAGAAU'
            'GCUACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGCUG'
            'CAAAAGAAGUGGGUAGUUUAACCUUUCGGGGAGGACGCUCACCACUUUGUGGUUCAUGACUG'
            'GGGUGAAGUCGUAACAAGGUAGCGCUAGGGGAACCUGGCGCUGGAUCACCUCCUUA',
            {'DATE': ['02-JUN-2014 (Rel. 121, Created)',
                      '04-FEB-2016 (Rel. 127, Last updated, Version 5)'],
             'DBSOURCE': 'SILVA-LSU; LK021130. SILVA-SSU; LK021130. MD5; '
                         'afd116bf2c1a13acbf40d63d82f0218c. BioSample; '
                         'SAMEA3865288.',
             'DEFINITION': 'Vibrio anguillarum 16S rRNA',
             'KEYWORDS': '.',
             'LOCUS': {'locus_name': 'LK021130.1:74067..75610:rRNA',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'genomic DNA',
                       'shape': 'linear',
                       'size': 1544,
                       'unit': 'bp',
                       'version': 1,
                       'date': '04-FEB-2016'},
             'PARENT_ACCESSION': 'LK021130.1',
             'VERSION': 'LK021130.1',
             'PROJECT_IDENTIFIER': 'Project:PRJEB5701;',
             'REFERENCE': [
                {'AUTHORS': 'Holm K.;',
                 'JOURNAL': 'Submitted (26-MAR-2014) to the INSDC. '
                            'Norstruct, Dept of Chemistry, University of '
                            'Tromso, Science Park 3, NO-9037 Tromso, NORWAY.',
                 'TITLE': ';',
                 'REFERENCE': '1'},
                {'AUTHORS': 'Holm K.O., Nilsson K., Hjerde E., Willassen '
                            'N.P., Milton D.L.;',
                 'JOURNAL': 'Stand Genomic Sci. 10:60-60(2015).',
                 'TITLE': '"Complete genome sequence of Vibrio anguillarum '
                          'strain NB10, a virulent isolate from the Gulf '
                          'of Bothnia";',
                 'REFERENCE': '2',
                 'PUBMED': '26380645'}],
             'CROSS_REFERENCE': [
                None,
                'DOI; 10.1186/s40793-015-0060-7. PUBMED; 26380645.'],
             'SOURCE': {
                'ORGANISM': 'Vibrio anguillarum',
                'taxonomy': 'Bacteria; Proteobacteria; Gammaproteobacteria; '
                            'Vibrionales; Vibrionaceae; Vibrio.'}},
            None,
            RNA)

        # get the feature level file without FT
        self.feature_level_fp = get_data_path(
                "embl_feature_level_record_no_FT")

        # get a genbank file in order to to file conversion
        self.genbank_fp = get_data_path('genbank_single_record')

        # a embl constructed sequence file path
        self.embl_constructed_fp = get_data_path("embl_constructed")

        # a simple embl version to perform embl->gb->embl conversion
        self.single_rna_simple_fp = get_data_path(
                "embl_single_record_simple")


class ReaderTests(EMBLIOTests):
    """Implements test for reading EMBL data"""

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
                                        r'Could not parse the ID line:.*'):
                _parse_id(line)

    # current status of protein support is described in issue-1499
    # https://github.com/biocore/scikit-bio/issues/1499
    def test_no_protein_support(self):
        """Testing no protein support for embl"""
        # TODO: add protein support

        # a fake protein line.
        handle = io.StringIO('ID   M14399; SV 1; linear; mRNA; STD; '
                             'PRO; 63 AA.\n//\n')

        with self.assertRaisesRegex(EMBLFormatError,
                                    r"There's no protein support for EMBL "
                                    "record"):
            # read a protein record
            Protein.read(handle)

        # return to 0
        handle.seek(0)

        with self.assertRaisesRegex(EMBLFormatError,
                                    r"There's no protein support for EMBL "
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

        # DNA, Sequence and RNA data contain newlines
        lines = [line+"\n" for line in lines if line != '']

        exp = {'AUTHORS': 'Gray G.L., Baldridge J.S., '
                          'McKeown K.S., Heyneker H.L., Chang C.N.;',
               'JOURNAL': 'Gene 39(2-3):247-254(1985).',
               'CROSS_REFERENCE': 'DOI; 10.1016/0378-1119(85)90319-1. '
                                  'PUBMED; 3912261.',
               'REFERENCE': '(bases 1 to 63)',
               'TITLE': '"Periplasmic production of correctly processed '
                        'human growth hormone in Escherichia coli: '
                        'natural and bacterial signal sequences are '
                        'interchangeable";',
               'PUBMED': '3912261'
               }

        # read reference
        obs = _parse_reference(lines)

        # See all differences
        self.maxDiff = None
        self.assertEqual(obs, exp)

    def test_parse_assembly(self):
        lines = """
AH   LOCAL_SPAN     PRIMARY_IDENTIFIER     PRIMARY_SPAN     COMP
AS   1-426          AC004528.1             18665-19090
AS   427-526        AC001234.2             1-100            c
AS   527-1000       TI55475028             not_available
""".split('\n')

        # DNA, Sequence and RNA data contain newlines
        lines = [line+"\n" for line in lines if line != '']

        exp = [
            {
             'local_span': '1-426',
             'primary_identifier': 'AC004528.1',
             'primary_span': '18665-19090',
             'comp': ''
            },
            {
             'local_span': '427-526',
             'primary_identifier': 'AC001234.2',
             'primary_span': '1-100',
             'comp': 'c'
            },
            {
             'local_span': '527-1000',
             'primary_identifier': 'TI55475028',
             'primary_span': 'not_available',
             'comp': ''
            }
        ]

        # read reference
        obs = _parse_assembly(lines)

        # See all differences
        self.maxDiff = None
        self.assertEqual(obs, exp)

    def test_parse_bad_assembly(self):
        """test for a wrong assembly line"""

        lines = """
AH   LOCAL_SPAN     PRIMARY_IDENTIFIER     PRIMARY_SPAN     COMP
AS   1-426          AC004528.1
""".split("\n")

        # DNA, Sequence and RNA data contain newlines
        lines = [line+"\n" for line in lines if line != '']

        with self.assertRaisesRegex(EMBLFormatError,
                                    r"Can't parse assembly line"):
            # read a malformed assembly record
            _parse_assembly(lines)

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

    # current status of protein support is described in issue-1499
    # https://github.com/biocore/scikit-bio/issues/1499
    def test_embl_to_protein(self):
        # TODO: add protein support
        i = 0
        # there is no support for protein at the moment
        # when protein support will be added, this code must work
        # exp = self.multi[i]
        # obs = _embl_to_protein(self.multi_fp, seq_num=i+1)
        # exp = Protein(exp[0], metadata=exp[1],
        #               lowercase=True, interval_metadata=exp[2])
        # self.assertEqual(exp, obs)

        with self.assertRaisesRegex(EMBLFormatError,
                                    r"There's no protein support for EMBL "
                                    "record"):
            # read a generic record
            _embl_to_protein(self.multi_fp, seq_num=i+1)

    # deal with feature-level-products: ignore feature table
    def test_feature_level_products(self):
        seq, md, imd, constructor = self.feature_level
        obs = _embl_to_rna(self.feature_level_fp)
        exp = constructor(seq, metadata=md,
                          lowercase=True, interval_metadata=imd)

        self.assertEqual(obs, exp)

    # deal with constructed sequences: ignore interval metadata
    def test_constructed_sequences(self):
        with self.assertRaisesRegex(
                EMBLFormatError,
                r"There's no support for embl CON record"):

            _embl_to_dna(self.embl_constructed_fp)


class WriterTests(EMBLIOTests):
    def test_serialize_id(self):
        for serialized, parsed in self.id:
            self.assertEqual(
                _serialize_id('ID', parsed), serialized[0] + '\n')

    def test_serialize_dbsource(self):
        """Serialize a complex dbsource entry"""

        # test with a complex uniprot dbsource
        exp = """DR   EMBL; AY548484; AAT09660.1; -; Genomic_DNA.
DR   RefSeq; YP_031579.1; NC_005946.1.
DR   ProteinModelPortal; Q6GZX4; -.
DR   SwissPalm; Q6GZX4; -.
DR   GeneID; 2947773; -.
DR   KEGG; vg:2947773; -.
DR   Proteomes; UP000008770; Genome.
"""

        # split by lines
        lines = [line+"\n" for line in exp.split("\n") if line != '']

        # parse objects
        parsed = _embl_parse_section_default(lines)

        # now serialize them
        obs = _serialize_dbsource("DR", parsed)

        # test objects
        self.assertEqual(obs, exp)

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

    # TODO: add support for protein
    # current status of protein support is described in issue-1499
    # https://github.com/biocore/scikit-bio/issues/1499
    def test_protein_to_embl(self):
        seq, md, imd, constructor = self.protein
        obj = constructor(seq, md, interval_metadata=imd)

        with io.StringIO() as fh:
            self.assertRaisesRegex(EMBLFormatError,
                                   r"There's no protein support for EMBL "
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

    def test_rna_to_embl_flp(self):
        """Test writing feature level products"""

        with io.StringIO() as fh:
            seq, md, imd, constructor = self.feature_level
            obj = constructor(seq, md, interval_metadata=imd, lowercase=True)
            _rna_to_embl(obj, fh)
            obs = fh.getvalue()

        with open(self.feature_level_fp) as fh:
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
    # current status of protein support is described in issue-1499
    # https://github.com/biocore/scikit-bio/issues/1499

    def test_roundtrip_sequence(self):
        with io.StringIO() as fh:
            _sequence_to_embl(_embl_to_sequence(self.single_rna_fp), fh)
            obs = fh.getvalue()

        with open(self.single_rna_fp) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


class Convertertest(EMBLIOTests):
    def test_gb_to_embl(self):
        genbank = DNA.read(self.genbank_fp, format="genbank")

        with io.StringIO() as fh:
            DNA.write(genbank, format="embl", file=fh)

            # EMBL can't deal with genbank version (ie M14399.1  GI:145229)
            # read embl data and write to gb
            fh.seek(0)
            embl = DNA.read(fh, format="embl")

        with io.StringIO() as fh:
            DNA.write(embl, format="genbank", file=fh)

            # read gb data
            obs = fh.getvalue()

        with open(self.genbank_fp) as fh:
            exp = fh.read()

        self.assertEqual(exp, obs)

    def test_embl_to_gb(self):
        # EMBL records have more features than genbank, (ex more than one date,
        # embl class, DOI cross references) so I can't convert an embl to gb
        # and then to embl keeping all those data. But I can start from
        # genbank record

        # do embl file -> embl object -> gb file -> gb object ->
        # embl file. Ensure that first and last files are identical
        embl = DNA.read(self.single_rna_simple_fp, format="embl")

        # "write" genbank record in a embl file
        with io.StringIO() as fh:
            DNA.write(embl, format="genbank", file=fh)

            # read genbank file
            fh.seek(0)
            genbank = DNA.read(fh, format="genbank")

        # "write" genbank record in a embl file
        with io.StringIO() as fh:
            DNA.write(genbank, format="embl", file=fh)

            # read file object
            obs = fh.getvalue()

        # test objects
        with open(self.single_rna_simple_fp) as fh:
            exp = fh.read()

        self.assertEqual(exp, obs)


if __name__ == '__main__':
    main()
