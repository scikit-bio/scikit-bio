# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import unittest

from skbio.io import PhylipFormatError
from skbio.io.format.phylip import (
    _tabular_msa_to_phylip, _phylip_to_tabular_msa, _phylip_sniffer)
from skbio import TabularMSA, DNA, RNA
from skbio.util import get_data_path


class TestSniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [get_data_path(e) for e in [
            'phylip_dna_3_seqs',
            'phylip_single_seq_long',
            'phylip_single_seq_short',
            'phylip_two_chunks',
            'phylip_variable_length_ids',
            'phylip_varied_whitespace_in_seqs',
            'phylip_whitespace_in_header_1',
            'phylip_whitespace_in_header_2',
            'phylip_whitespace_in_header_3',
        ]]

        # negative tests for sniffer don't include
        # phylip_invalid_empty_line_between_seqs, phylip_invalid_too_few_seqs,
        # phylip_invalid_too_many_seqs - because sniffer only reads first seq
        self.negatives = [get_data_path(e) for e in [
            'empty',
            'whitespace_only',
            'phylip_invalid_empty_line_after_header',
            'phylip_invalid_empty_line_before_header',
            'phylip_invalid_header_too_long',
            'phylip_invalid_header_too_short',
            'phylip_invalid_no_header',
            'phylip_invalid_seq_too_long',
            'phylip_invalid_seq_too_short',
            'phylip_invalid_zero_seq_len',
            'phylip_invalid_zero_seqs',
        ]]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_phylip_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_phylip_sniffer(fp), (False, {}))


class TestReaders(unittest.TestCase):
    def setUp(self):
        self.valid_configurations = [
            ([get_data_path('phylip_dna_3_seqs')],
             [('..ACC-GTTGG..', 'd1'), ('TTACCGGT-GGCC', 'd2'),
              ('.-ACC-GTTGC--', 'd3')]
             ),

            ([get_data_path('phylip_single_seq_long')],
             [('..ACC-GTTGG..AATGC.C----', 'foo')]
             ),

            ([get_data_path('phylip_single_seq_short')],
             [('-', '')]
             ),

            ([get_data_path('phylip_two_chunks'),
              get_data_path('phylip_varied_whitespace_in_seqs'),
              get_data_path('phylip_whitespace_in_header_1'),
              get_data_path('phylip_whitespace_in_header_2'),
              get_data_path('phylip_whitespace_in_header_3'),
              ],
             [('..ACC-GTTGG..AATGC.C', 'foo'), ('TTACCGGT-GGCCTA-GCAT', 'bar')]
             ),

            ([get_data_path('phylip_variable_length_ids')],
             [('.-ACGT', ''), ('TGCA-.', 'a'), ('.ACGT-', 'bb'),
              ('TGCA-.', '1'), ('AAAAAA', 'abcdefghij'),
              ('GGGGGG', 'ab def42ij')]
             ),

        ]

        self.positive_fps = list(map(get_data_path, [
            'phylip_dna_3_seqs',
            'phylip_single_seq_long',
            'phylip_single_seq_short',
            'phylip_two_chunks',
            'phylip_variable_length_ids',
            'phylip_varied_whitespace_in_seqs',
            'phylip_whitespace_in_header_1',
            'phylip_whitespace_in_header_2',
            'phylip_whitespace_in_header_3',
        ]))

        self.invalid_files = [(get_data_path(e[0]), e[1], e[2]) for e in [
            ('empty', PhylipFormatError,
             r'This file is empty.'),

            ('whitespace_only', PhylipFormatError,
             r'Found non-header line .*: ""'),

            ('phylip_invalid_empty_line_after_header', PhylipFormatError,
             r'Empty lines are not allowed.'),

            ('phylip_invalid_empty_line_before_header', PhylipFormatError,
             r'Found non-header line .*: ""'),

            ('phylip_invalid_empty_line_between_seqs', PhylipFormatError,
             r'Empty lines are not allowed.'),

            ('phylip_invalid_header_too_long', PhylipFormatError,
             r'Found non-header line .*: "2 20 extra_text"'),

            ('phylip_invalid_header_too_short', PhylipFormatError,
             r'Found non-header line .*: " 20"'),

            ('phylip_invalid_no_header', PhylipFormatError,
             r'Found non-header line .*: "foo .*"'),

            ('phylip_invalid_seq_too_long', PhylipFormatError,
             r'The length of sequence foo is not 20 as specified .*.'),

            ('phylip_invalid_seq_too_short', PhylipFormatError,
             r'The length of sequence foo is not 20 as specified .*.'),

            ('phylip_invalid_too_few_seqs', PhylipFormatError,
             r'The number of sequences is not .* as specified .*.'),

            ('phylip_invalid_too_many_seqs', PhylipFormatError,
             r'The number of sequences is not .* as specified in the header.'),

            ('phylip_invalid_zero_seq_len', PhylipFormatError,
             r'The number of sequences and the length must be positive.'),

            ('phylip_invalid_zero_seqs', PhylipFormatError,
             r'The number of sequences and the length must be positive.'),
        ]]

    def test_phylip_to_tabular_msa_invalid_files(self):
        for fp, error_type, error_msg_regex in self.invalid_files:
            with self.assertRaisesRegex(error_type, error_msg_regex):
                _phylip_to_tabular_msa(fp, constructor=DNA)

    def test_phylip_to_tabular_msa_no_constructor(self):
        with self.assertRaisesRegex(ValueError, '`constructor`'):
            _phylip_to_tabular_msa(get_data_path('phylip_dna_3_seqs'))

    def test_phylip_to_tabular_msa_valid_files(self):
        for valid_files, components in self.valid_configurations:
            for valid in valid_files:
                observed = _phylip_to_tabular_msa(valid, constructor=DNA)

                expected_seqs = []
                expected_index = []
                for seq, ID in components:
                    expected_seqs.append(DNA(seq))
                    expected_index.append(ID)
                expected = TabularMSA(expected_seqs, index=expected_index)

                self.assertEqual(observed, expected)


class TestWriters(unittest.TestCase):
    def setUp(self):
        # ids all same length, seqs longer than 10 chars
        dna_3_seqs = TabularMSA([
            DNA('..ACC-GTTGG..', metadata={'id': "d1"}),
            DNA('TTACCGGT-GGCC', metadata={'id': "d2"}),
            DNA('.-ACC-GTTGC--', metadata={'id': "d3"})], minter='id')

        # id lengths from 0 to 10, with mixes of numbers, characters, and
        # spaces. sequences are shorter than 10 chars
        variable_length_ids = TabularMSA([
            DNA('.-ACGT', metadata={'id': ''}),
            DNA('TGCA-.', metadata={'id': 'a'}),
            DNA('.ACGT-', metadata={'id': 'bb'}),
            DNA('TGCA-.', metadata={'id': '1'}),
            DNA('AAAAAA', metadata={'id': 'abcdefghij'}),
            DNA('GGGGGG', metadata={'id': 'ab def42ij'})], minter='id')

        # sequences with 20 chars = exactly two chunks of size 10
        two_chunks = TabularMSA([
            DNA('..ACC-GTTGG..AATGC.C', metadata={'id': 'foo'}),
            DNA('TTACCGGT-GGCCTA-GCAT', metadata={'id': 'bar'})], minter='id')

        # single sequence with more than two chunks
        single_seq_long = TabularMSA([
            DNA('..ACC-GTTGG..AATGC.C----', metadata={'id': 'foo'})],
            minter='id')

        # single sequence with only a single character (minimal writeable
        # alignment)
        single_seq_short = TabularMSA([DNA('-', metadata={'id': ''})],
                                      minter='id')

        # alignments that can be written in phylip format
        self.objs = [dna_3_seqs, variable_length_ids, two_chunks,
                     single_seq_long, single_seq_short]
        self.fps = map(get_data_path,
                       ['phylip_dna_3_seqs', 'phylip_variable_length_ids',
                        'phylip_two_chunks', 'phylip_single_seq_long',
                        'phylip_single_seq_short'])

        # alignments that cannot be written in phylip format, paired with their
        # expected error message regexps
        self.invalid_objs = [
            # no seqs
            (TabularMSA([]), 'one sequence'),

            # no positions
            (TabularMSA([DNA('', metadata={'id': "d1"}),
                         DNA('', metadata={'id': "d2"})]), 'one position'),

            # ids too long
            (TabularMSA([RNA('ACGU', metadata={'id': "foo"}),
                         RNA('UGCA', metadata={'id': "alongsequenceid"})],
                        minter='id'),
             '10.*alongsequenceid')
        ]

    def test_write(self):
        for fp, obj in zip(self.fps, self.objs):
            fh = io.StringIO()
            _tabular_msa_to_phylip(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with io.open(fp) as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_write_invalid_alignment(self):
        for invalid_obj, error_msg_regexp in self.invalid_objs:
            fh = io.StringIO()
            with self.assertRaisesRegex(PhylipFormatError, error_msg_regexp):
                _tabular_msa_to_phylip(invalid_obj, fh)

            # ensure nothing was written to the file before the error was
            # thrown
            obs = fh.getvalue()
            fh.close()
            self.assertEqual(obs, '')


if __name__ == '__main__':
    unittest.main()
