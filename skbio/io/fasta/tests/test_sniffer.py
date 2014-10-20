# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

from skbio.io.fasta._fasta import _fasta_sniffer
from skbio.io.fasta._fastaqual import _qual_sniffer
from skbio.util import get_data_path


class SnifferTests(TestCase):
    def setUp(self):
        self.fasta_positive_fps = map(get_data_path, [
            'fasta_3_seqs_defaults',
            'fasta_max_width_1',
            'fasta_single_bio_seq_non_defaults',
            'fasta_single_prot_seq_non_defaults',
            'fasta_3_seqs_non_defaults',
            'fasta_max_width_5',
            'fasta_single_dna_seq_defaults',
            'fasta_single_rna_seq_defaults',
            'fasta_description_newline_replacement_empty_str',
            'fasta_multi_seq',
            'fasta_single_dna_seq_non_defaults',
            'fasta_single_rna_seq_non_defaults',
            'fasta_description_newline_replacement_multi_char',
            'fasta_prot_seqs_odd_labels',
            'fasta_single_nuc_seq_defaults',
            'fasta_single_seq',
            'fasta_id_whitespace_replacement_empty_str',
            'fasta_sequence_collection_different_type',
            'fasta_single_nuc_seq_non_defaults',
            'fasta_id_whitespace_replacement_multi_char',
            'fasta_single_bio_seq_defaults',
            'fasta_single_prot_seq_defaults',
            'fasta_10_seqs',
            'fasta_invalid_after_10_seqs'
        ])

        self.fasta_negative_fps = map(get_data_path, [
            'empty',
            'whitespace_only',
            'fasta_invalid_missing_header',
            'fasta_invalid_blank_line',
            'fasta_invalid_whitespace_only_line',
            'fasta_invalid_missing_seq_data_first',
            'fasta_invalid_missing_seq_data_middle',
            'fasta_invalid_missing_seq_data_last',
            'fasta_invalid_legacy_format',
            'fasta_id_whitespace_replacement_none',
            'fasta_description_newline_replacement_none'
        ])

        self.qual_positive_fps = map(get_data_path, [
            'qual_3_seqs_defaults',
        ])

    def test_fasta_positives(self):
        for fp in self.fasta_positive_fps:
            self.assertEqual(_fasta_sniffer(fp), (True, {}))

    def test_fasta_negatives(self):
        for fp in self.fasta_negative_fps + self.qual_positive_fps:
            self.assertEqual(_fasta_sniffer(fp), (False, {}))

    def test_qual_positives(self):
        for fp in self.qual_positive_fps:
            self.assertEqual(_qual_sniffer(fp), (True, {}))

    def test_qual_negatives(self):
        # none of the positive or negative fasta files should be identified as
        # qual
        for fp in self.fasta_positive_fps + self.fasta_negative_fps:
            self.assertEqual(_qual_sniffer(fp), (False, {}))


if __name__ == '__main__':
    main()
