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

# TODO: add this type of exception
from skbio.io import EMBLFormatError

from skbio.io.format.embl import (
    _embl_sniffer, _parse_id)

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
            'whitespace_only']))

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

class ReaderTests(EMBLIOTests):
    """Implements test for reading EMBL data"""
    
    # TODO: implement test to deal with all EMBL methods
    def test_parse_id(self):
        """Parse ID record (first line of embl format)"""
        for serialized, parsed in self.id:
            self.assertEqual(_parse_id(serialized), parsed)

