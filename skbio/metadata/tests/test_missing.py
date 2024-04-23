# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import pandas.testing as pdt
import numpy as np


from skbio.metadata.missing import series_encode_missing, series_extract_missing


class RoundTripMixin:
    def check_roundtrip(self, real_value, dtype):
        notna_exp = [real_value]
        series = pd.Series(notna_exp + self.missing_terms)

        encoded = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded)

        self.assertEqual(encoded.dtype, dtype)
        # the non-null side of the series
        self.assertEqual(list(encoded[encoded.notna()]), notna_exp)
        # the null end (but in the original vocabulary)
        pdt.assert_series_equal(missing, series[1:].astype(object))

    def test_roundtrip_float(self):
        self.check_roundtrip(0.05, float)

    def test_roundtrip_string(self):
        self.check_roundtrip('hello', object)

    def test_roundtrip_int(self):
        self.check_roundtrip(42, float)

    def test_roundtrip_bool(self):
        self.check_roundtrip(True, object)

    def test_roundtrip_all_missing_object(self):
        expected = [None, float('nan')] + self.missing_terms
        series = pd.Series(expected, dtype=object)

        encoded = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded)

        self.assertEqual(encoded.dtype, object)
        pdt.assert_series_equal(missing, series.astype(object))


class TestISNDC(RoundTripMixin, unittest.TestCase):
    def setUp(self):
        self.enum = 'INSDC:missing'
        self.missing_terms = ['not applicable', 'missing', 'not collected',
                              'not provided', 'restricted access']


class TestOmitted(RoundTripMixin, unittest.TestCase):
    def setUp(self):
        self.enum = 'blank'
        self.missing_terms = [None, float('nan')]

    # test_roundtrip_all_missing_float is not possible with other schemes
    def test_roundtrip_all_missing_float(self):
        expected = [None, float('nan')] + self.missing_terms
        series = pd.Series(expected, dtype=float)

        encoded = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded)

        self.assertEqual(encoded.dtype, float)
        pdt.assert_series_equal(missing, series.astype(object))


class TestError(RoundTripMixin, unittest.TestCase):
    def setUp(self):
        self.enum = 'no-missing'
        self.missing_terms = []

    # no missing values, so bool and int are not object and float
    def test_roundtrip_bool(self):
        self.check_roundtrip(True, bool)

    def test_roundtrip_int(self):
        self.check_roundtrip(42, np.int64)

    def test_roundtrip_all_missing_object(self):
        with self.assertRaisesRegex(ValueError, 'Missing values.*name=None'):
            super().test_roundtrip_all_missing_object()

if __name__ == '__main__':
    unittest.main()
