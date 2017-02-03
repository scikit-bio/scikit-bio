# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import DNA, RNA, Sequence
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path

# Module specific execption and functions
from skbio.io import EMBLFormatError

from skbio.io.format.embl import (
    _embl_sniffer, _parse_id, _parse_reference, _embl_to_generator, 
    _get_embl_section)

# TODO: implement those methods
#    _genbank_to_generator, _genbank_to_sequence,
#    _genbank_to_dna, _genbank_to_rna, _genbank_to_protein,
#    _parse_locus, _parse_reference,
#    _generator_to_genbank, _sequence_to_genbank,
#    _protein_to_genbank, _rna_to_genbank, _dna_to_genbank,
#    _serialize_locus)

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
            # (feature level record: http://www.ebi.ac.uk/ena/browse/feature-level-products
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
            # TODO: a Uniprot record?
        
        # define a single DNA record (with no interval metadata)
        # M14399; SV 1; linear; mRNA; STD; PRO; 63 BP.
        self.single = (
            'gtgaaacaaagcactattgcactggctgtcttaccgttactgtttacccctgtgacaaaagcc',
            {'ID': {'accession': 'M14399',
                    'class': 'STD',
                    'division': 'PRO',
                    'mol_type': 'mRNA',
                    'shape': 'linear',
                    'size': 63,
                    'unit': 'bp',
                    'version': 1}},
            None,
            RNA)
        
        # define a single DNA record uppercase (filepath)
        self.single_upper_fp = get_data_path('embl_single_record_upper')
        
        # define a single RNA record file path
        self.single_rna_fp = get_data_path('embl_single_record')
        
        # define a interval metadata (see skbio.metadata.IntervalMetadata)
        imd = IntervalMetadata(63)
        
        # then add interval object to interval metadata. Add source
        imd.add([(0, 63)],
                [(False, False)],
                {'db_xref': 'taxon:562',
                 'mol_type': 'mRNA',
                 'organism': 'Escherichia coli',
                 'type': 'source',
                 'strand': '+',
                 '__location': '1..63'})
        
        imd.add([(0, 63)],
                [(False, True)], # the second True is beacause exact location is not known
                {'phase': 0,
                 'db_xref': ['GOA:P00634',
                             'InterPro:IPR001952',
                             'InterPro:IPR017849',
                             'InterPro:IPR017850',
                             'InterPro:IPR018299',
                             'PDB:1AJA',
                             'PDB:1AJB',
                             'PDB:1AJC',
                             'PDB:1AJD',
                             'PDB:1ALH',
                             'PDB:1ALI',
                             'PDB:1ALJ',
                             'PDB:1ALK',
                             'PDB:1ANI',
                             'PDB:1ANJ',
                             'PDB:1B8J',
                             'PDB:1ED8',
                             'PDB:1ED9',
                             'PDB:1ELX',
                             'PDB:1ELY',
                             'PDB:1ELZ',
                             'PDB:1EW8',
                             'PDB:1EW9',
                             'PDB:1HJK',
                             'PDB:1HQA',
                             'PDB:1KH4',
                             'PDB:1KH5',
                             'PDB:1KH7',
                             'PDB:1KH9',
                             'PDB:1KHJ',
                             'PDB:1KHK',
                             'PDB:1KHL',
                             'PDB:1KHN',
                             'PDB:1URA',
                             'PDB:1URB',
                             'PDB:1Y6V',
                             'PDB:1Y7A',
                             'PDB:2ANH',
                             'PDB:2G9Y',
                             'PDB:2GA3',
                             'PDB:2MLX',
                             'PDB:2MLY',
                             'PDB:2MLZ',
                             'PDB:3BDF',
                             'PDB:3BDG',
                             'PDB:3BDH',
                             'PDB:3CMR',
                             'PDB:3DPC',
                             'PDB:3DYC',
                             'PDB:3TG0',
                             'PDB:4KM4',
                             'PDB:4YR1',
                             'PDB:5C66',
                             'PDB:5GAD',
                             'PDB:5GAF',
                             'PDB:5GAG',
                             'PDB:5GAH',
                             'PDB:5JTL',
                             'PDB:5JTM',
                             'PDB:5JTN',
                             'PDB:5JTO',
                             'PDB:5JTP',
                             'UniProtKB/Swiss-Prot:P00634'],
                 '__location': '1..>63',
                 'strand': '+',
                 'note': 'alkaline phosphatase signal peptide',
                 'protein_id': 'AAA23431.1',
                 'transl_table': '11',
                 'translation': 'MKQSTIALAVLPLLFTPVTKA',
                 'type': 'CDS'})
        
        # define a full RNA record
        self.single_rna_fp = get_data_path('embl_single_record')
        
        self.single_rna = (
            'gtgaaacaaagcactattgcactggctgtcttaccgttactgtttacccctgtgacaaaagcc',
            {'LOCUS': {'accession': 'M14399',
                       'class': 'STD',
                       'division': 'PRO',
                       'mol_type': 'mRNA',
                       'shape': 'linear',
                       'size': 63,
                       'unit': 'bp',
                       'version': 1},
             'ACCESSION': 'M14399', # accessions (could be more than one)
             'DATE': ["16-JUL-1988 (Rel. 16, Created)", 
                      "02-SEP-1999 (Rel. 60, Last updated, Version 3)"],
             'DESCRIPTION': "E.coli alkaline phosphatase signal mRNA, 5' end.",
             'KEYWORDS': "alkaline phosphatase; signal peptide.",
             'ORGANISM': "Escherichia coli",
             'TAXONOMY': "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; "
                         "Enterobacteriaceae; Escherichia.",
             'REFERENCE': [{'AUTHORS': 'Gray G.L., Baldridge J.S., '
                                       'McKeown K.S., Heyneker H.L., Chang C.N.',
                            'JOURNAL': 'Gene 39(2-3):247-254(1985).',
                            'PUBMED': '3912261',
                            'REFERENCE': '1-63',
                            'TITLE': '"Periplasmic production of correctly processed '
                                     'human growth hormone in Escherichia coli: '
                                     'natural and bacterial signal sequences are '
                                     'interchangeable"'
                           }],
            },
            imd,
            RNA)

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
               'CROSS_REFERENCE': 'DOI; 10.1016/0378-1119(85)90319-1. PUBMED; 3912261.',
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
                    raise EMBLFormatError("Key {0} isn't defined in embl.KEYS_2_SECTIONS".format(err))
                    
if __name__ == '__main__':
    main()

