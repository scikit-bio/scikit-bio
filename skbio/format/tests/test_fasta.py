#!/usr/bin/env python
"""Tests for FASTA sequence format writer.
"""
from unittest import TestCase, main
from skbio.format.fasta import fasta_from_sequences, fasta_from_alignment


class FastaTests(TestCase):

    """Tests for Fasta writer.
    """

    def setUp(self):
        """Setup for Fasta tests."""
        self.strings = ['AAAA', 'CCCC', 'gggg', 'uuuu']
        self.labels = ['1st', '2nd', '3rd', '4th']
        self.infos = ["Dog", "Cat", "Mouse", "Rat"]
        self.fasta_no_label = '>0\nAAAA\n>1\nCCCC\n>2\ngggg\n>3\nuuuu'
        self.fasta_with_label =\
            '>1st\nAAAA\n>2nd\nCCCC\n>3rd\nGGGG\n>4th\nUUUU'
        self.fasta_with_label_lw2 =\
            '>1st\nAA\nAA\n>2nd\nCC\nCC\n>3rd\nGG\nGG\n>4th\nUU\nUU'
        self.alignment_dict = {'1st': 'AAAA', '2nd': 'CCCC', '3rd': 'GGGG',
                               '4th': 'UUUU'}
        self.fasta_with_label_species =\
            '>1st:Dog\nAAAA\n>2nd:Cat\nCCCC\n>3rd:Mouse\nGGGG\n>4th:Rat\nUUUU'


    def test_fasta_from_sequence(self):
        """should return correct fasta string."""
        self.assertEqual(fasta_from_sequences(''), '')
        self.assertEqual(fasta_from_sequences(self.strings),
                         self.fasta_no_label)

    def test_fasta_from_alignment_dict(self):
        """should return correct fasta string for a dictionary"""
        self.assertEqual(fasta_from_alignment({}), '')
        self.assertEqual(fasta_from_alignment(self.alignment_dict),
                         self.fasta_with_label)
        self.assertEqual(fasta_from_alignment(self.alignment_dict,
                                              line_wrap=2),
                         self.fasta_with_label_lw2)

if __name__ == "__main__":
    main()
