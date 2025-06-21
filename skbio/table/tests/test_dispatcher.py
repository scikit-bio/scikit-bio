# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, skipIf

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt

from skbio._config import set_config
from skbio.table import Table
from skbio.util import get_package
from skbio.util._testing import assert_data_frame_almost_equal

from skbio.table._dispatcher import _create_table, _create_table_1d, _ingest_table


pl = get_package("polars", raise_error=False)
# polars.testing isn't imported along with polars. Therefore one needs to import it
# separately.
if pl is not None:
    plt = get_package("polars.testing")
adt = get_package("anndata", raise_error=False)


class TestPandas(TestCase):
    def setUp(self):
        set_config("output", "pandas")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def test_create_table_no_backend(self):
        obs = _create_table(data=self.data, columns=self.columns, index=self.index)
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = _create_table(self.data)
        exp = pd.DataFrame(self.data)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_numpy_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index)
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = _create_table_1d(self.data_1d)
        exp = pd.Series(self.data_1d)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_numpy_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="numpy")
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table_1d(self.data, backend="nonsense")


class TestNumpy(TestCase):
    def setUp(self):
        set_config("output", "numpy")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_config("output", "pandas")

    def test_create_table_no_backend(self):
        obs = _create_table(data=self.data, columns=self.columns, index=self.index)
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = _create_table(self.data)
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_pandas_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index)
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = _create_table_1d(self.data_1d)
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="numpy")
        exp = np.array(self.data_1d)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_pandas_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table_1d(self.data, backend="nonsense")


@skipIf(pl is None, "Polars is not available for unit tests.")
class TestPolars(TestCase):
    def setUp(self):
        set_config("output", "polars")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_config("output", "pandas")

    def test_create_table_no_backend(self):
        obs = _create_table(data=self.data, columns=self.columns, index=self.index)
        exp = pl.DataFrame(self.data, schema=self.columns)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_no_optionals(self):
        obs = _create_table(self.data)
        exp = pl.DataFrame(self.data)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_same_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="polars"
        )
        exp = pl.DataFrame(self.data, schema=self.columns)
        plt.assert_frame_equal(obs, exp)

    def test_create_table_numpy_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_pandas_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="pandas"
        )
        exp = pd.DataFrame(self.data, columns=self.columns, index=self.index)
        assert_data_frame_almost_equal(obs, exp)

    def test_create_table_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table(self.data, backend="nonsense")

    def test_create_table_1d_no_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index)
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_no_optionals(self):
        obs = _create_table_1d(self.data_1d)
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_same_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="polars")
        exp = pl.Series(self.data_1d)
        plt.assert_series_equal(obs, exp)

    def test_create_table_1d_pandas_backend(self):
        obs = _create_table_1d(self.data_1d, index=self.index, backend="pandas")
        exp = pd.Series(self.data_1d, index=self.index)
        pdt.assert_series_equal(obs, exp)

    def test_create_table_1d_numpy_backend(self):
        obs = _create_table(
            data=self.data, columns=self.columns, index=self.index, backend="numpy"
        )
        exp = np.array(self.data)
        npt.assert_array_equal(obs, exp)

    def test_create_table_1d_bad_backend(self):
        with self.assertRaisesRegex(ValueError, "Unsupported backend: 'nonsense'"):
            _create_table_1d(self.data, backend="nonsense")


class TestBIOM(TestCase):
    def setUp(self):
        self.data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.samples = ["A", "B", "C"]
        self.features = ["f1", "f2", "f3"]

    def test_biom_input(self):
        tbl = Table(
            self.data.T, observation_ids=self.features, sample_ids=self.samples
        )
        data, row_ids, col_ids = _ingest_table(tbl)
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, self.samples)
        self.assertEqual(col_ids, self.features)

    def test_biom_input_pass_ids(self):
        tbl = Table(
            self.data.T, observation_ids=self.features, sample_ids=self.samples
        )
        data, row_ids, col_ids = _ingest_table(
            tbl, sample_ids=self.samples, feature_ids=self.features
        )
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, self.samples)
        self.assertEqual(col_ids, self.features)


@skipIf(adt is None, "Anndata is not available for unit tests.")
class TestAnndata(TestCase):
    def setUp(self):
        self.data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.samples = pd.DataFrame(index=["A", "B", "C"])
        self.features = pd.DataFrame(index=["f1", "f2", "f3"])

    def test_anndata_input(self):
        tbl = adt.AnnData(self.data, obs=self.samples, var=self.features)
        data, row_ids, col_ids = _ingest_table(tbl)
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, list(self.samples.index))
        self.assertEqual(col_ids, list(self.features.index))

    def test_anndata_input_pass_ids(self):
        tbl = adt.AnnData(self.data, obs=self.samples, var=self.features)
        data, row_ids, col_ids = _ingest_table(
            tbl, sample_ids=list(self.samples.index), feature_ids=list(self.features.index)
        )
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, list(self.samples.index))
        self.assertEqual(col_ids, list(self.features.index))


if __name__ == "__main__":
    main()
