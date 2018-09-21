# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import numpy as np

from skbio.util import get_data_path, assert_data_frame_almost_equal
from skbio.io import BLAST7FormatError
from skbio.io.format.blast7 import _blast7_to_data_frame, _blast7_sniffer


class TestBLAST7Sniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [get_data_path(e) for e in [
            'blast7_default_single_line',
            'blast7_default_multi_line',
            'blast7_custom_minimal',
            'blast7_custom_single_line',
            'blast7_custom_multi_line',
            'blast7_custom_mixed_nans',
            'blast7_invalid_differing_fields',
            'blast7_invalid_no_data',
            'blast7_invalid_too_many_columns',
            'legacy9_and_blast7_default',
            'legacy9_invalid_too_many_columns',
            'legacy9_mixed_nans',
            'legacy9_multi_line',
            'legacy9_single_line']]

        self.negatives = [get_data_path(e) for e in [
            'blast7_invalid_gibberish',
            'blast7_invalid_for_sniffer',
            'blast7_invalid_for_sniffer_2',
            'empty']]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_blast7_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_blast7_sniffer(fp), (False, {}))


class TestBlast7Reader(unittest.TestCase):
    def test_default_valid_single_line(self):
        fp = get_data_path('blast7_default_single_line')
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1', 'subject2', 100.00, 8.0, 0.0, 0.0, 1.0,
                             8.0, 3.0, 10.0, 9e-05, 16.9]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

        fp = get_data_path('legacy9_single_line')
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1', 'subject1', 90.00, 7.0, 1.0, 0.0, 0.0,
                             8.0, 4.0, 10.0, 1e-05, 15.5]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

    def test_default_valid_multi_line(self):
        fp = get_data_path('blast7_default_multi_line')
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1', 'subject2', 70.00, 5.0, 0.0, 0.0, 7.0,
                             60.0, 3.0, 100.0, 9e-05, 10.5],
                            ['query1', 'subject2', 30.00, 8.0, 0.0, 0.0, 6.0,
                             15.0, 1.0, 100.0, 0.053, 12.0],
                            ['query1', 'subject2', 90.00, 2.0, 0.0, 0.0, 9.0,
                             35.0, 2.0, 100.0, 0.002, 8.3]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

        fp = get_data_path('legacy9_multi_line')
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1', 'subject1', 90.00, 7.0, 1.0, 0.0, 0.0,
                             8.0, 4.0, 10.0, 1e-05, 15.5],
                            ['query1', 'subject1', 70.00, 8.0, 0.0, 1.0, 0.0,
                             9.0, 5.0, 7.0, 0.231, 7.8],
                            ['query1', 'subject1', 90.00, 5.0, 1.0, 1.0, 0.0,
                             0.0, 2.0, 10.0, 0.022, 13.0]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

    def test_default_valid_mixed_output(self):
        fp = get_data_path('legacy9_and_blast7_default')
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query2', 'subject2', 100.00, 8.0, 0.0, 1.0, 0.0,
                             9.0, 3.0, 10.0, 2e-05, 9.8],
                            ['query2', 'subject1', 70.00, 9.0, 1.0, 0.0, 1.0,
                             8.0, 4.0, 9.0, 0.025, 11.7]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

    def test_custom_valid_minimal(self):
        fp = get_data_path("blast7_custom_minimal")
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1']], columns=['qseqid'])
        assert_data_frame_almost_equal(df, exp)

    def test_custom_valid_single_line(self):
        fp = get_data_path("blast7_custom_single_line")
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([['query1', 100.00, 100.00, 8.0, 0.0, 16.9, 8.0,
                             'PAAWWWWW']],
                           columns=['qseqid', 'ppos', 'pident', 'length',
                                    'sgi', 'bitscore', 'qend', 'qseq'])
        assert_data_frame_almost_equal(df, exp)

    def test_custom_valid_multi_line(self):
        fp = get_data_path("blast7_custom_multi_line")
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([[1.0, 8.0, 3.0, 10.0, 8.0, 0.0, 1.0, 'query1',
                             'subject2'],
                            [2.0, 5.0, 2.0, 15.0, 8.0, 0.0, 2.0, 'query1',
                             'subject2'],
                            [1.0, 6.0, 2.0, 12.0, 8.0, 0.0, 1.0, 'query1',
                             'subject2']],
                           columns=['qstart', 'qend', 'sstart', 'send',
                                    'nident', 'mismatch', 'sframe',
                                    'qaccver', 'saccver'])
        assert_data_frame_almost_equal(df, exp)

    def test_custom_valid_mixed_nans(self):
        fp = get_data_path("blast7_custom_mixed_nans")
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([[0.0, np.nan, 8.0, 13.0, 1.0, 1.0, np.nan,
                             'subject2'],
                            [np.nan, 0.0, 8.0, np.nan, 1.0, 1.0, 'query1',
                            np.nan]],
                           columns=['qgi', 'sgi', 'qlen', 'slen', 'qframe',
                                    'sframe', 'qseqid', 'sseqid'])
        assert_data_frame_almost_equal(df, exp)

    def test_legacy9_valid_mixed_nans(self):
        fp = get_data_path("legacy9_mixed_nans")
        df = _blast7_to_data_frame(fp)
        exp = pd.DataFrame([[np.nan, 'subject1', np.nan, 7.0, 1.0, 0.0, np.nan,
                             8.0, 4.0, 10.0, np.nan, 15.5],
                            ['query2', 'subject1', 90.00, 8.0, np.nan, 0.0,
                             0.0, 8.0, np.nan, 9.0, 1e-05, np.nan]],
                           columns=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
        assert_data_frame_almost_equal(df, exp)

    def test_differing_fields_error(self):
        fp = get_data_path("blast7_invalid_differing_fields")
        with self.assertRaisesRegex(
                BLAST7FormatError,
                r"Fields \[.*'qseqid', .*'sseqid', .*'qstart'\]"
                r" do.*\[.*'qseqid', .*'sseqid', .*'score'\]"):
            _blast7_to_data_frame(fp)
        fp = get_data_path("legacy9_invalid_differing_fields")
        with self.assertRaisesRegex(
                BLAST7FormatError,
                r"Fields \[.*'qseqid', .*'sseqid', .*'qstart'\]"
                r" do.*\[.*'qseqid', .*'sseqid', "
                r".*'sallseqid'\]"):
            _blast7_to_data_frame(fp)

    def test_no_data_error(self):
        fp = get_data_path("blast7_invalid_gibberish")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"File contains no"):
            _blast7_to_data_frame(fp)
        fp = get_data_path("blast7_invalid_no_data")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"File contains no"):
            _blast7_to_data_frame(fp)
        fp = get_data_path("empty")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"File contains no"):
            _blast7_to_data_frame(fp)

    def test_wrong_amount_of_columns_error(self):
        fp = get_data_path("blast7_invalid_too_many_columns")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"Number of fields.*\(2\)"):
            _blast7_to_data_frame(fp)
        fp = get_data_path("legacy9_invalid_too_many_columns")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"Number of fields.*\(12\)"):
            _blast7_to_data_frame(fp)

    def test_unrecognized_field_error(self):
        fp = get_data_path("blast7_invalid_unrecognized_field")
        with self.assertRaisesRegex(BLAST7FormatError,
                                    r"Unrecognized field \(.*'sallid'\)"):
            _blast7_to_data_frame(fp)


if __name__ == '__main__':
    unittest.main()
