#!/usr/bin/env python
"""Tests for FASTA sequence format writer.
"""
from unittest import TestCase, main
from skbio.format.sequences.fasta import (fasta_from_sequences,
                                          fasta_from_alignment)
from skbio.core.sequence import DNASequence, BiologicalSequence


class FastaTests(TestCase):

    """Tests for Fasta writer.
    """

    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ['AAAA', 'CCCC', 'gggg', 'uuuu']
        self.fasta_no_label = '>0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu'
        self.fasta_with_label =\
            '>1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU'
        self.fasta_with_label_lw2 =\
            '>1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU'
        self.alignment_dict = {'1st': 'AAAA', '2nd': 'CCCC', '3rd': 'GGGG',
                               '4th': 'UUUU'}
        self.sequence_objects_a = [DNASequence('ACTCGAGATC', 'seq1'),
                                   DNASequence('GGCCT', 'seq2')]
        self.sequence_objects_b = [BiologicalSequence('ACTCGAGATC', 'seq1'),
                                   BiologicalSequence('GGCCT', 'seq2')]

    def test_fasta_from_sequence_objects(self):
        """Check FASTA files are created correctly off of sequence objects"""
        self.assertEqual(fasta_from_sequences(self.sequence_objects_a),
                         FASTA_STRING)

        self.assertEqual(fasta_from_sequences(self.sequence_objects_b),
                         FASTA_STRING)

    def test_fasta_from_sequences(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_sequences(''), '')
        self.assertEqual(fasta_from_sequences(self.strings),
                         self.fasta_no_label)

    def test_fasta_from_alignment(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_alignment({}), '')
        self.assertEqual(fasta_from_alignment(self.alignment_dict),
                         self.fasta_with_label)
        self.assertEqual(fasta_from_alignment(self.alignment_dict,
                                              line_wrap=2),
                         self.fasta_with_label_lw2)

FASTA_STRING = '>seq1\nACTCGAGATC\n>seq2\nGGCCT'

if __name__ == "__main__":
    main()
