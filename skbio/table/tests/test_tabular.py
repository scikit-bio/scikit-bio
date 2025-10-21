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

from skbio.table._tabular import _create_table, _create_table_1d, _ingest_table


# import optional dependencies
pl = get_package("polars", raise_error=False)
# polars.testing isn't imported along with polars. Therefore one needs to import it
# separately.
if pl is not None:
    plt = get_package("polars.testing")
adt = get_package("anndata", raise_error=False)
torch = get_package("torch", raise_error=False)
jnp = get_package("jax.numpy", raise_error=False)


class TestIngest(TestCase):
    def setUp(self):
        self.data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        self.data_1d = self.data[0]
        self.samples = ["A", "B", "C"]
        self.features = ["f1", "f2", "f3"]

    def test_ingest_numpy(self):
        table = self.data
        obs = _ingest_table(table)
        self.assertIs(obs[0], table)
        self.assertIsNone(obs[1])
        self.assertIsNone(obs[2])

        obs = _ingest_table(table, self.samples, self.features)
        self.assertIs(obs[0], table)
        self.assertIs(obs[1], self.samples)
        self.assertIs(obs[2], self.features)

    def test_ingest_pandas(self):
        table = pd.DataFrame(self.data, index=self.samples, columns=self.features)
        obs = _ingest_table(table)
        self.assertIsInstance(obs[0], np.ndarray)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(obs[1].tolist(), self.samples)
        self.assertListEqual(obs[2].tolist(), self.features)

        # no samples / features
        table = pd.DataFrame(self.data)
        obs = _ingest_table(table)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(obs[1].tolist(), list(range(obs[0].shape[0])))
        self.assertListEqual(obs[2].tolist(), list(range(obs[0].shape[1])))

        # override samples / features
        obs = _ingest_table(table, self.samples, self.features)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(obs[1], self.samples)
        self.assertListEqual(obs[2], self.features)

    @skipIf(pl is None, "Polars is not available for unit tests.")
    def test_ingest_polars(self):
        table = pl.DataFrame(self.data, schema=self.features)
        obs = _ingest_table(table, sample_ids=self.samples)
        self.assertIsInstance(obs[0], np.ndarray)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(obs[1], self.samples)
        self.assertListEqual(obs[2].names(), self.features)

        # no samples, override features
        table = pl.DataFrame(self.data)
        obs = _ingest_table(table, feature_ids=self.features)
        npt.assert_array_equal(obs[0], self.data)
        self.assertIsNone(obs[1])
        self.assertListEqual(obs[2], self.features)

    def test_ingest_biom(self):
        table = Table(
            self.data.T, observation_ids=self.features, sample_ids=self.samples
        )
        obs = _ingest_table(table)
        self.assertIsInstance(obs[0], np.ndarray)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(obs[1], self.samples)
        self.assertListEqual(obs[2], self.features)

    @skipIf(adt is None, "Anndata is not available for unit tests.")
    def test_ingest_anndata(self):
        table = adt.AnnData(
            self.data,
            obs=pd.DataFrame(index=self.samples),
            var=pd.DataFrame(index=self.features),
        )
        obs = _ingest_table(table)
        self.assertIsInstance(obs[0], np.ndarray)
        npt.assert_array_equal(obs[0], self.data)
        self.assertListEqual(list(obs[1]), self.samples)
        self.assertListEqual(list(obs[2]), self.features)

    def test_ingest_sequence(self):
        table = self.data.tolist()
        obs = _ingest_table(table)
        npt.assert_array_equal(obs[0], self.data)
        self.assertIsNone(obs[1])
        self.assertIsNone(obs[2])

        obs = _ingest_table(tuple(table))
        npt.assert_array_equal(obs[0], self.data)

        obs = _ingest_table(tuple(tuple(x) for x in table))
        npt.assert_array_equal(obs[0], self.data)

    @skipIf(torch is None, "PyTorch is not available for unit tests.")
    def test_ingest_torch(self):
        table = torch.tensor(self.data)
        obs = _ingest_table(table)
        npt.assert_array_equal(obs[0], self.data)
        self.assertIsNone(obs[1])
        self.assertIsNone(obs[2])

    @skipIf(jnp is None, "JAX is not available for unit tests.")
    def test_ingest_jax(self):
        table = jnp.asarray(self.data)
        obs = _ingest_table(table)
        npt.assert_array_equal(obs[0], self.data)
        self.assertIsNone(obs[1])
        self.assertIsNone(obs[2])

    def test_ingest_error(self):
        msg = "'int' is not a supported table format."
        with self.assertRaises(TypeError) as cm:
            _ingest_table(123)
        self.assertEqual(str(cm.exception), msg)

        msg = "'str' is not a supported table format."
        with self.assertRaises(TypeError) as cm:
            _ingest_table("hello")
        self.assertEqual(str(cm.exception), msg)

        msg = "'type' is not a supported table format."
        with self.assertRaises(TypeError) as cm:
            _ingest_table(TypeError)
        self.assertEqual(str(cm.exception), msg)

        msg = "Input table has less than 2 dimensions."
        with self.assertRaises(ValueError) as cm:
            _ingest_table(np.array([1, 2, 3]), expand=False)
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            _ingest_table([1, 2, 3], expand=False)
        self.assertEqual(str(cm.exception), msg)

        msg = "Input table has 3 samples whereas 2 sample IDs were provided."
        with self.assertRaises(ValueError) as cm:
            _ingest_table(self.data, sample_ids=["A", "B"])
        self.assertEqual(str(cm.exception), msg)

        msg = "Input table has 3 features whereas 2 feature IDs were provided."
        with self.assertRaises(ValueError) as cm:
            _ingest_table(self.data, feature_ids=["f2", "f3"])
        self.assertEqual(str(cm.exception), msg)


class TestPandas(TestCase):
    def setUp(self):
        set_config("table_output", "pandas")
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
        set_config("table_output", "numpy")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_config("table_output", "pandas")

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
        set_config("table_output", "polars")
        self.data = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        self.data_1d = self.data[0]
        self.index = ["A", "B", "C"]
        self.columns = ["f1", "f2", "f3"]

    def tearDown(self):
        set_config("table_output", "pandas")

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
        self.assertListEqual(list(row_ids), list(self.samples.index))
        self.assertListEqual(list(col_ids), list(self.features.index))

    def test_anndata_input_pass_ids(self):
        tbl = adt.AnnData(self.data, obs=self.samples, var=self.features)
        data, row_ids, col_ids = _ingest_table(
            tbl, sample_ids=list(self.samples.index), feature_ids=list(self.features.index)
        )
        npt.assert_array_equal(data, self.data)
        self.assertListEqual(list(row_ids), list(self.samples.index))
        self.assertListEqual(list(col_ids), list(self.features.index))


if __name__ == "__main__":
    main()
