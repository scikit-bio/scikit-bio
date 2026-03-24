# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

import numpy as np

from skbio.util import get_data_path
from skbio.io.format.phylip_dm import (
    _parse_header,
    _parse_line,
    _phylip_dm_sniffer,
    _phylip_dm_to_distance_matrix,
    _matrix_to_phylip_dm,
)
from skbio.io import PhylipDMFormatError
from skbio.stats.distance import DistanceMatrix


class TestHelpers(unittest.TestCase):
    def test_parse_header(self):
        self.assertEqual(_parse_header("5\n"), 5)
        self.assertEqual(_parse_header("   5\n"), 5)
        self.assertEqual(_parse_header("5  \n"), 5)

        msg = "Header line must be a single integer."
        for header in ("  3 7\n", "1.25\n", "hello\n"):
            with self.assertRaises(PhylipDMFormatError) as cm:
                _parse_header(header)
            self.assertEqual(str(cm.exception), msg)

        msg = "The number of objects must be positive."
        for n in ("-5\n", "0\n"):
            with self.assertRaises(PhylipDMFormatError) as cm:
                _parse_header(f"{n}\n")
            self.assertEqual(str(cm.exception), msg)

    def test_parse_line(self):
        for line in ("", "\n", "   \n", "\t\r\n"):
            with self.assertRaises(PhylipDMFormatError) as cm:
                _parse_line(line, 0, 1, False, False, float)
            self.assertEqual(str(cm.exception), "Empty lines are not allowed.")

        with self.assertRaises(PhylipDMFormatError) as cm:
            _parse_line("            1.23", 0, 1, False, True, float)
        self.assertEqual(str(cm.exception), "Empty IDs are not allowed.")

        with self.assertRaises(PhylipDMFormatError) as cm:
            _parse_line("sample\t0.25\there\t0.17", 7, 3, True, False, float)
        # The following test doesn't work in Python 3.10. Skipping...
        # self.assertEqual(str(cm.exception), (
        #     "Non-numeric cell values encountered in line 9."))


class TestSniffer(unittest.TestCase):
    def test_positives(self):
        # A typical file in lower layout with 14 objects. Header line has preceding
        # whitespaces. ID column has a width of 12, whereas all IDs themselves are
        # within 10 characters, leaving at least 2 spaces separating IDs from values.
        # This is often seen in real-world PHYLIP files. Yet, it is considered as
        # "relaxed" in scikit-bio. Columns are separated by whitespaces but one column
        # separation was intentionally filled with tabs for testing purpose.
        fp = get_data_path("phylip_dm_valid_lt.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="lower", strict=False))

        # square, relaxed, 5 objects
        fp = get_data_path("phylip_dm_valid_sq.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="square", strict=False))

        # lower, relaxed, 4 objects; IDs are incremental integers
        fp = get_data_path("phylip_dm_simple_lt.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="lower", strict=False))

        # square, relaxed, 4 objects, as above
        fp = get_data_path("phylip_dm_simple_sq.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="square", strict=False))

        # square, strict, 4 objects; first and last IDs are exactly 10 characters
        # (therefore touching values without any separator)
        fp = get_data_path("phylip_dm_good_simple_strict_square.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="square", strict=True))

        # square, strict, 4 objects; first and last IDs are exactly 10 characters
        # (therefore touching values without any separator)
        fp = get_data_path("phylip_dm_good_simple_strict_square.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="square", strict=True))

        # same as above, but only last ID is exactly 10 characters
        fp = get_data_path("phylip_dm_good_simple_strict_square_reader.dist")
        obs = _phylip_dm_sniffer(fp)
        self.assertTrue(obs[0])
        self.assertDictEqual(obs[1], dict(layout="square", strict=True))

    def test_negatives(self):
        negative_files = [
            "empty",
            "whitespace_only",
            "phylip_dm_invalid_empty_line_after_header.dist",
            "phylip_dm_invalid_empty_line_before_header.dist",
            "phylip_dm_sq_invalid_empty_line_after_header.dist",
            "phylip_dm_sq_invalid_empty_line_before_header.dist",
            "phylip_dm_invalid_header_too_long.dist",
            "phylip_dm_invalid_no_header.dist",
            "phylip_dm_invalid_no_dists.dist",
            "phylip_dm_invalid_zero_header.dist",
            "phylip_dm_invalid_too_many_columns.dist",
            "phylip_dm_invalid_too_few_columns_sq.dist",
            "phylip_dm_invalid_wrong_number_dists_lt.dist",
        ]
        for file in negative_files:
            self.assertEqual(_phylip_dm_sniffer(get_data_path(file)), (False, {}))


class TestReaders(unittest.TestCase):
    def test_read_positives(self):
        # a small, square matrix
        fp = get_data_path("phylip_dm_valid_sq.dist")
        kwargs = dict(layout="square", strict=False)
        exp = np.array([
            [0.0, 1.0, 2.0, 3.0, 3.0],
            [1.0, 0.0, 2.0, 3.0, 3.0],
            [2.0, 2.0, 0.0, 3.0, 3.0],
            [3.0, 3.0, 3.0, 0.0, 1.0],
            [3.0, 3.0, 3.0, 1.0, 0.0],
        ])
        obs = _phylip_dm_to_distance_matrix(fp, **kwargs)
        self.assertTrue((obs.data == exp).all())

        # a large, lower matrix
        fp = get_data_path("phylip_dm_valid_lt.dist")
        kwargs = dict(layout="lower", strict=False)
        exp = np.array([
            [0.0, 1.7043, 2.0235, 2.1378, 1.5232, 1.8261, 1.9182, 2.0039, 1.9431,
                1.9663, 2.0593, 1.6664, 1.732, 1.7101,],
            [1.7043, 0.0, 1.1901, 1.3287, 1.2423, 1.2508, 1.2536, 1.3066, 1.2827,
                1.3296, 1.2005, 1.346, 1.3757, 1.3956,],
            [2.0235, 1.1901, 0.0, 1.2905, 1.3199, 1.3887, 1.4658, 1.4826, 1.4502,
                1.8708, 1.5356, 1.4577, 1.7803, 1.6661,],
            [2.1378, 1.3287, 1.2905, 0.0, 1.7878, 1.3137, 1.3788, 1.3826, 1.4543,
                1.6683, 1.6606, 1.5935, 1.7119, 1.7599,],
            [1.5232, 1.2423, 1.3199, 1.7878, 0.0, 1.0642, 1.1124, 0.9832, 1.0629,
                0.9228, 1.0681, 0.9127, 1.0635, 1.0557,],
            [1.8261, 1.2508, 1.3887, 1.3137, 1.0642, 0.0, 0.1022, 0.2061, 0.3895,
                0.8035, 0.7239, 0.7278, 0.7899, 0.6933,],
            [1.9182, 1.2536, 1.4658, 1.3788, 1.1124, 0.1022, 0.0, 0.2681, 0.393,
                0.7109, 0.729, 0.7412, 0.8742, 0.7118,],
            [2.0039, 1.3066, 1.4826, 1.3826, 0.9832, 0.2061, 0.2681, 0.0, 0.3665,
                0.8132, 0.7894, 0.8763, 0.8868, 0.7589,],
            [1.9431, 1.2827, 1.4502, 1.4543, 1.0629, 0.3895, 0.393, 0.3665, 0.0,
                0.7858, 0.714, 0.7966, 0.8288, 0.8542,],
            [1.9663, 1.3296, 1.8708, 1.6683, 0.9228, 0.8035, 0.7109, 0.8132, 0.7858,
                0.0, 0.7095, 0.5959, 0.6213, 0.5612,],
            [2.0593, 1.2005, 1.5356, 1.6606, 1.0681, 0.7239, 0.729, 0.7894, 0.714,
                0.7095, 0.0, 0.4604, 0.5065, 0.47,],
            [1.6664, 1.346, 1.4577, 1.5935, 0.9127, 0.7278, 0.7412, 0.8763, 0.7966,
                0.5959, 0.4604, 0.0, 0.3502, 0.3097,],
            [1.732, 1.3757, 1.7803, 1.7119, 1.0635, 0.7899, 0.8742, 0.8868, 0.8288,
                0.6213, 0.5065, 0.3502, 0.0, 0.2712,],
            [1.7101, 1.3956, 1.6661, 1.7599, 1.0557, 0.6933, 0.7118, 0.7589, 0.8542,
                0.5612, 0.47, 0.3097, 0.2712, 0.0,],
        ])
        obs = _phylip_dm_to_distance_matrix(fp, **kwargs)
        self.assertTrue((obs.data == exp).all())

    def test_read_dtype(self):
        fn = _phylip_dm_to_distance_matrix
        fp = get_data_path("phylip_dm_valid_sq.dist")
        kwargs = dict(layout="square", strict=False)
        exp = np.array([
            [0.0, 1.0, 2.0, 3.0, 3.0],
            [1.0, 0.0, 2.0, 3.0, 3.0],
            [2.0, 2.0, 0.0, 3.0, 3.0],
            [3.0, 3.0, 3.0, 0.0, 1.0],
            [3.0, 3.0, 3.0, 1.0, 0.0],
        ])

        # float64 (default)
        for dtype in (None, float, np.float64, "float64"):
            obs = _phylip_dm_to_distance_matrix(fp, dtype=dtype, **kwargs)
            self.assertTrue((obs.data == exp).all())
            self.assertEqual(obs.dtype, np.float64)

        # float32
        exp = exp.astype(np.float32)
        for dtype in (np.float32, "float32"):
            obs = _phylip_dm_to_distance_matrix(fp, dtype=dtype, **kwargs)
            self.assertTrue((obs.data == exp).all())
            self.assertEqual(obs.dtype, np.float32)

        # invalid data type
        for dtype in ("xyz", int, "float16", "uint8", np.float16):
            self.assertRaises(
                TypeError, _phylip_dm_to_distance_matrix, fp, dtype=dtype, **kwargs
            )

    def test_error_empty_file(self):
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(get_data_path("empty"))
        self.assertEqual(str(cm.exception), "This file is empty.")

    def test_error_wrong_row_number(self):
        fp = get_data_path("phylip_dm_invalid_too_many_rows.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="square")
        self.assertEqual(str(cm.exception), (
            "The number of rows in the matrix body exceeds 5 as specified in the "
            "header."
        ))

        fp = get_data_path("phylip_dm_invalid_too_few_rows.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="square")
        self.assertEqual(str(cm.exception), (
            "The number of rows in the matrix body (4) does not match the number of "
            "objects specified in the header (5)."
        ))

    def test_error_wrong_column_number(self):
        fp = get_data_path("phylip_dm_invalid_too_many_columns.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="square")
        self.assertEqual(str(cm.exception), (
            "The number of distance values (6) in line 2 does not match the "
            "expectation (5, as specified in the header)."
        ))

        fp = get_data_path("phylip_dm_invalid_too_few_columns_sq.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="square")
        self.assertEqual(str(cm.exception), (
            "The number of distance values (4) in line 2 does not match the "
            "expectation (5, as specified in the header)."
        ))

        fp = get_data_path("phylip_dm_invalid_wrong_number_dists_lt.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="lower")
        self.assertEqual(str(cm.exception), (
            "The number of distance values (1) in line 2 does not match the "
            "expectation (0, which is line number - 2)."
        ))

    def test_error_matrix_data_parsed_as_id(self):
        # A special case where a strict object ID is exactly 10 characters and it
        # touches the values. If parsed in relaxed mode, it will raise.
        fp = get_data_path("phylip_dm_good_simple_strict_square_reader.dist")
        with self.assertRaises(PhylipDMFormatError) as cm:
            _phylip_dm_to_distance_matrix(fp, layout="square", strict=False)
        self.assertEqual(str(cm.exception), (
            "The number of distance values (3) in line 5 does not match the "
            "expectation (4, as specified in the header)."
        ))


class TestWriters(unittest.TestCase):
    def setUp(self):
        self.dm_condensed = DistanceMatrix([1, 2, 3, 4, 5, 6], condensed=True)
        self.dm_square = DistanceMatrix(
            [[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]]
        )
        self.exp_lower = (
            "4\n0\n1\t1.0\n2\t2.0\t4.0\n3\t3.0\t5.0\t6.0\n"
        )
        self.exp_square = (
            "4\n0\t0.0\t1.0\t2.0\t3.0\n1\t1.0\t0.0\t4.0\t5.0\n"
            "2\t2.0\t4.0\t0.0\t6.0\n3\t3.0\t5.0\t6.0\t0.0\n"
        )

    def test_condensed_to_lower(self):
        """Test writing condensed DistanceMatrix to lower triangular layout."""
        fh = io.StringIO()
        self.dm_condensed.write(fh, format="phylip_dm", layout="lower")
        result = fh.getvalue()
        fh.close()
        self.assertEqual(result, self.exp_lower)

    def test_condensed_to_square(self):
        """Test writing condensed DistanceMatrix to square layout."""
        fh = io.StringIO()
        self.dm_condensed.write(fh, format="phylip_dm", layout="square")
        result = fh.getvalue()
        fh.close()
        self.assertEqual(result, self.exp_square)

    def test_square_to_lower(self):
        """Test writing square DistanceMatrix to lower triangular layout."""
        fh = io.StringIO()
        self.dm_square.write(fh, format="phylip_dm", layout="lower")
        result = fh.getvalue()
        fh.close()
        self.assertEqual(result, self.exp_lower)

    def test_square_to_square(self):
        """Test writing square DistanceMatrix to square layout."""
        fh = io.StringIO()
        self.dm_square.write(fh, format="phylip_dm", layout="square")
        result = fh.getvalue()
        fh.close()
        self.assertEqual(result, self.exp_square)

    def test_not_enough_objects(self):
        fh = io.StringIO()
        msg = (
            "DistanceMatrix can only be written in PHYLIP format if there are at "
            "least two samples in the matrix."
        )
        with self.assertRaises(PhylipDMFormatError) as e:
            _matrix_to_phylip_dm(DistanceMatrix([]), fh, "\t", layout="lower")
        self.assertEqual(str(e.exception), msg)
        fh.close()

    def test_invalid_layout(self):
        """Test writing DistanceMatrix to invalid layout."""
        fh = io.StringIO()
        msg = "Upper triangular layout is currently not supported."
        with self.assertRaises(PhylipDMFormatError) as e:
            self.dm_square.write(fh, format="phylip_dm", layout="upper")
        self.assertEqual(str(e.exception), msg)
        msg = "'xyz' is not a supported matrix layout."
        with self.assertRaises(PhylipDMFormatError) as e:
            self.dm_square.write(fh, format="phylip_dm", layout="xyz")
        self.assertEqual(str(e.exception), msg)
        fh.close()


if __name__ == "__main__":
    unittest.main()
