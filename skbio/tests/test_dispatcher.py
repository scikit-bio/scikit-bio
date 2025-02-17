# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt

try:
    import polars as pl
    import polars.testing as plt
except ImportError:
    has_polars = False
else:
    has_polars = True

from skbio._dispatcher import create_table, create_table_1d
from skbio._config import get_option, set_option
from skbio.util._testing import assert_data_frame_almost_equal


class TestPandas(unittest.TestCase):
    def setUp(self):
        set_option("tabular_backend", "pandas")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def test_create_table_no_backend(self):
        obs = create_table(data=self.data, columns=self.columns, index=self.index)
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = create_table(self.data)
        exp = pd.DataFrame(self.data)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_numpy_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index)
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = create_table_1d(self.data_1d)
        exp = pd.Series(self.data_1d)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_numpy_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="numpy")
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table_1d(self.data, backend="nonsense")


class TestNumpy(unittest.TestCase):
    def setUp(self):
        set_option("tabular_backend", "numpy")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_option("tabular_backend", "pandas")

    def test_create_table_no_backend(self):
        obs = create_table(data=self.data, columns=self.columns, index=self.index)
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = create_table(self.data)
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_pandas_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index)
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = create_table_1d(self.data_1d)
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="numpy")
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_pandas_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table_1d(self.data, backend="nonsense")


@unittest.skipIf(not has_polars, "Polars is not available for unit tests.")
class TestPolars(unittest.TestCase):
    def setUp(self):
        set_option("tabular_backend", "polars")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_option("tabular_backend", "pandas")

    def test_create_table_no_backend(self):
        obs = create_table(data=self.data, columns=self.columns, index=self.index)
        exp = pl.DataFrame(self.data, schema=self.columns)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = create_table(self.data)
        exp = pl.DataFrame(self.data)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="polars"
        )
        exp = pl.DataFrame(self.data, schema=self.columns)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_numpy_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_pandas_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index)
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = create_table_1d(self.data_1d)
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="polars")
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_pandas_backend(self):
        obs = create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_numpy_backend(self):
        obs = create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            create_table_1d(self.data, backend="nonsense")


if __name__ == "__main__":
    unittest.main()
