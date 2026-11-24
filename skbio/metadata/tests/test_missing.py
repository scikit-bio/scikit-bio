# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from packaging.version import Version

import pandas as pd
import pandas.testing as pdt
import numpy as np

from skbio.metadata.missing import series_encode_missing, series_extract_missing

PANDAS_3 = Version(pd.__version__) >= Version("3.0.0")

def _is_string_dtype(dtype):
    """Check if dtype is a string-like dtype (handles Pandas 2.x and 3.x)."""
    if dtype == object:
        return True
    dtype_str = str(dtype)
    return dtype_str == 'str' or dtype_str.startswith('string') or 'StringDtype' in dtype_str


class RoundTripMixin:
    def check_roundtrip(self, real_value, dtype):
        notna_exp = [real_value]
        series = pd.Series(notna_exp + self.missing_terms)

        encoded, mask = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded, mask)

        # For string types, Pandas 3.x uses StringDtype instead of object
        if dtype == object:
            self.assertTrue(_is_string_dtype(encoded.dtype),
                f"Expected string-like dtype, got {encoded.dtype}")
        else:
            self.assertEqual(encoded.dtype, dtype)
        # the non-null side of the series
        self.assertEqual(list(encoded[encoded.notna()]), notna_exp)
        # the null end (but in the original vocabulary)
        # Compare values loosely to handle None vs nan differences in Pandas 3.x
        pdt.assert_series_equal(missing, series[1:], check_dtype=False)

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

        encoded, mask = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded, mask)

        self.assertTrue(_is_string_dtype(encoded.dtype),
            f"Expected string-like dtype, got {encoded.dtype}")
        # Compare values loosely to handle None vs nan differences in Pandas 3.x
        pdt.assert_series_equal(missing, series.astype(object), check_dtype=False)


class TestISNDC(RoundTripMixin, unittest.TestCase):
    def setUp(self):
        self.enum = 'INSDC:missing'
        self.missing_terms = ['not applicable', 'missing', 'not collected',
                              'not provided', 'restricted access']

    def test_roundtrip_float(self):
        self.check_roundtrip(0.05, object)

    def test_roundtrip_int(self):
        self.check_roundtrip(42, object)


class TestOmitted(RoundTripMixin, unittest.TestCase):
    def setUp(self):
        self.enum = 'blank'
        self.missing_terms = [None, float('nan')]

    # test_roundtrip_all_missing_float is not possible with other schemes
    def test_roundtrip_all_missing_float(self):
        expected = [None, float('nan')] + self.missing_terms
        series = pd.Series(expected, dtype=float)

        encoded, mask = series_encode_missing(series, self.enum)
        missing = series_extract_missing(encoded, mask)

        self.assertEqual(encoded.dtype, float)
        pdt.assert_series_equal(missing, series)#.astype(object))


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
