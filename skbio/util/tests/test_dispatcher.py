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

try:
    import polars as pl
    import polars.testing as plt
except (ImportError, ModuleNotFoundError):
    has_polars = False
else:
    has_polars = True

try:
    import anndata as adt
except (ImportError, ModuleNotFoundError):
    has_anndata = False
else:
    has_anndata = True

from skbio.table import Table
from skbio.util.config._dispatcher import create_table, create_table_1d, ingest_array
from skbio.util import set_config
from skbio.util._testing import assert_data_frame_almost_equal


class TestPandas(TestCase):
    def setUp(self):
        set_config("output", "pandas")
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


@skipIf(not has_polars, "Polars is not available for unit tests.")
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


class TestBIOM(TestCase):
    def setUp(self):
        self.data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.samples = ["A", "B", "C"]
        self.features = ["f1", "f2", "f3"]

    def test_biom_input(self):
        # need to transpose to ensure proper orientation of data
        tbl = Table(
            self.data.T, observation_ids=self.features, sample_ids=self.samples
        ).transpose()
        with self.assertWarnsRegex(
            UserWarning,
            "BIOM format uses samples as columns and features as rows. Most "
            "scikit-bio functions expect samples as rows and features as columns. "
            "Please ensure your input is in the correct orientation.\n",
        ):
            data, row_ids, col_ids = ingest_array(tbl)
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, self.samples)
        self.assertEqual(col_ids, self.features)

    def test_biom_input_pass_ids(self):
        # need to transpose to ensure proper orientation of data
        tbl = Table(
            self.data.T, observation_ids=self.features, sample_ids=self.samples
        ).transpose()
        with self.assertWarnsRegex(
            UserWarning,
            "BIOM format uses samples as columns and features as rows. Most "
            "scikit-bio functions expect samples as rows and features as columns. "
            "Please ensure your input is in the correct orientation.\n",
        ):
            data, row_ids, col_ids = ingest_array(
                tbl, row_ids=self.samples, col_ids=self.features
            )
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, self.samples)
        self.assertEqual(col_ids, self.features)


@skipIf(not has_anndata, "Anndata is not available for unit tests.")
class TestAnndata(TestCase):
    def setUp(self):
        self.data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.samples = pd.DataFrame(index=["A", "B", "C"])
        self.features = pd.DataFrame(index=["f1", "f2", "f3"])

    def test_anndata_input(self):
        tbl = adt.AnnData(self.data, obs=self.samples, var=self.features)
        data, row_ids, col_ids = ingest_array(tbl)
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, list(self.samples.index))
        self.assertEqual(col_ids, list(self.features.index))

    def test_anndata_input_pass_ids(self):
        tbl = adt.AnnData(self.data, obs=self.samples, var=self.features)
        data, row_ids, col_ids = ingest_array(
            tbl, row_ids=list(self.samples.index), col_ids=list(self.features.index)
        )
        npt.assert_array_equal(data, self.data)
        self.assertEqual(row_ids, list(self.samples.index))
        self.assertEqual(col_ids, list(self.features.index))


if __name__ == "__main__":
    main()
